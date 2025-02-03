/**
 * Runs on page load, sets up the sketchers,
 * switches to the active editor and sets up the click event to switch in the future,
 * and if save is enabled, calls functions to reload, lock, and set up the autosave if needed
 */
$(async function () {
  /*
    Need to move the functions that are used in search.js into utils to only have to load this script when making the actual reaction page
     */
  showSketcherLoadingCircle();
  await setupNewKetcherSketcher();

  // sleep used to allow sketchers to load scripts and make js Objects
  await sleep(1000);
  // sketcher changes when user clicks radio button
  await switchActiveEditor();

  $('input[name="sketcher-select"]').click(function () {
    switchActiveEditor();
  });
  if (checkIfSaveEnabled()) {
    //
    reloadReaction();
    sketcherLockHandler();
    setupSketcherAutosave();
  }
  hideSketcherLoadingCircle();
  // Use Promise.race to give setupNewMarvinSketcher() a maximum of 2 seconds
  // try {
  //   await Promise.race([
  //     setupNewMarvinSketcher(),
  //     new Promise((_, reject) => {
  //       setTimeout(() => {
  //         reject(new Error("Timeout after 2 seconds"));
  //       }, 2000);
  //     }),
  //   ]);
  // } catch (error) {
  //   abortMarvinSketcherCreation();
  //   console.error(error);
  //   // Handle errors or take appropriate action
  // }
});

/**
 * Event listener for the retrosynthesis-export button click
 * Exports the SMILES of either a single structure or the product to the retrosynthesis page
 * @return {Promise<void>}
 */
async function exportDataToRetrosynthesis() {
  let smiles = await exportSmilesFromActiveEditor();
  if (smiles.includes(">>")) {
    smiles = smiles.split(">>")[1];
  }
  window.open(`/retrosynthesis/${smiles}`, "_blank");
}

/**
 * Event listener for the "Example" button click
 * Calls a function to load example smiles into the selected sketcher.
 */
function loadExampleSmiles() {
  let selectedSketcher = $('input[name="sketcher-select"]:checked').attr("id");
  if (selectedSketcher === "marvin-select") {
    marvinExampleSmiles();
  } else if (selectedSketcher === "ketcher-select") {
    ketcherExampleSmiles();
  }
}

/**
 * Event listener for the "Submit" button click
 * Sends an ajax request of all the data required to make the reaction table.
 * Then handles the response from the ajax request and either makes the reaction table,
 * prompts for compound data or alerts error
 */
async function createReactionTable() {
  let smiles;
  try {
    smiles = await exportSmilesFromActiveEditor();
  } catch (error) {
    alert("An error occurred:\n\n" + error.message);
  }
  let rxn = await exportRXNFromActiveEditor();
  $("#js-reaction-rxn").val(rxn);
  await exportImageFromActiveEditor();
  sketcherDataLossHandler();
  let [reactants, products] =
    reactionSmilesToReactantsAndProductsSmiles(smiles);
  $(".loading-bar").css("display", "block");
  let smilesNew = removeReagentsFromSmiles(smiles);
  $("#js-reaction-smiles").val(smilesNew);
  let workgroup = getVal($("#js-active-workgroup"));
  let workbook = getVal($("#js-active-workbook"));
  let reaction_id = getVal($("#js-reaction-id"));
  let demo_mode = getVal($("#js-demo"));
  let tutorial = getVal($("#js-tutorial"));
  let polymer_mode = $('input[id="polymer-mode-select"]').prop("checked");
  let polymer_indices = identifyPolymers(await exportRXNFromActiveEditor());
  smiles = replaceSmilesSymbols(smiles);

  // Asynchronous request to _process in routes.py
  $.ajax({
    type: "GET",
    contentType: "application/json;charset-utf-08",
    dataType: "json",
    url:
      "/_process?reactants=" +
      reactants +
      "&products=" +
      products +
      "&reactionSmiles=" +
      smiles +
      "&demo=" +
      demo_mode +
      "&tutorial=" +
      tutorial +
      "&workgroup=" +
      workgroup +
      "&workbook=" +
      workbook +
      "&reaction_id=" +
      reaction_id +
      "&polymer=" +
      polymer_mode +
      "&polymerIndices=" +
      polymer_indices,
  }).done(function (data) {
    $(".loading-bar").css("display", "none");
    if (data.novelCompound) {
      // one of the compounds is not in database - must be added as novel compound
      $("#reaction-table-div").html(data.reactionTable).show();
      sketcherNovelCompoundInputValidate();
    } else if (data.error) {
      $("#errorAlert").text(data.error).show();
      $("#reaction-table-div").hide();
    } else if (data.reactionTable === "Demo") {
      alert(
        "One or more reactants or products are not in the database. Adding novel compounds is not available in Demo mode.",
      );
    } else {
      generateReactionTable(data);
    }
  });
  // when ajax is done
  document.body.className = "";
}

/**
 * Adds the reaction table rendered template from the response to the page HTML.
 * Calls the required functions to enable reaction table autofill/styling
 * @param {Object} data - The AJAX response object
 */
function generateReactionTable(data) {
  $("#reaction-table-div").html(data.reactionTable).show(); // Sends data to the reaction table
  initialiseReactionTable();
  $("#action-summary").show();
  $("#js-reaction-name2").val($("#js-reaction-name").val());
  // initial mass and equivalents to highlight
  let limitingReactantTableNumber = String(
    Number($("input[name='reactant-limiting']:checked").val()),
  );
  let colorRoundedReactantMassID =
    "#js-reactant-rounded-mass" + limitingReactantTableNumber;
  $("#js-reactant-equivalent" + limitingReactantTableNumber).val(1);
  autoChangeRequiredStyling2(colorRoundedReactantMassID);
  let numberOfReactants = Number($("#js-number-of-reactants").val());
  for (i = 1; i <= numberOfReactants; i++) {
    autoChangeRequiredStyling2("#js-reactant-equivalent" + i);
    autoChangeRequiredStyling2("#js-reactant-physical-form" + i);
  }
  let numberOfProducts = Number($("#js-number-of-products").val());
  for (i = 1; i <= numberOfProducts; i++) {
    autoChangeRequiredStyling2("#js-product-physical-form" + i);
  }
}
/**
 * Alerts user that continuing with their current input may overwrite existing data
 */
function sketcherDataLossHandler() {
  const $reactionDiv = $("#reaction-table-div");
  const reactionDivChildren = $reactionDiv.children().length;
  const firstChildNodeID = $reactionDiv.children().first().id;
  let tutorial_mode = $("#js-tutorial").val();
  // dont ask if in tutorial mode or if novel compound submission form is present
  if (tutorial_mode !== "yes") {
    if (
      reactionDivChildren !== 0 &&
      firstChildNodeID !== "js-novel-compound-input-form"
    ) {
      const text =
        "Please note that any reaction data already inputted will be lost. Do you wish to continue?";
      if (confirm(text) === false) {
        return;
      }
    }
  }
}

/**
 * Calls functions to reload the reaction sketcher, and reaction table if data for it is found in the json object
 */
async function reloadReaction() {
  // disable select sketcher buttons until SMILES are loaded into sketcher
  let $sketcherSelectRadioButtons = $('input[name="sketcher-select"]');
  $sketcherSelectRadioButtons.prop("disabled", true);

  const response = await getPolymerMode();
  if (response === true) {
    await reloadSketcherFromRXN();
  } else {
    await reloadSketcherFromSmiles();
  }
  $sketcherSelectRadioButtons.prop("disabled", false);
  $("#js-load-status").val("loaded");
  // show compound if reloading a reaction where the reaction table has previously been loaded and description is in json
  let js_reaction_table_data = JSON.parse(
    document.getElementById("js-reaction-table-data").value,
  );
  if (Object.keys(js_reaction_table_data).includes("reaction_description")) {
    await createReactionTable();
  }
}

/**
 * Loads the reloaded reaction smiles to the active sketcher
 */
async function reloadSketcherFromSmiles() {
  let $reactionSmiles = $("#js-reaction-smiles");
  let reloadedReaction = $reactionSmiles.val();
  if (reloadedReaction) {
    let selectedSketcher = $('input[name="sketcher-select"]:checked').attr(
      "id",
    );
    if (selectedSketcher === "marvin-select") {
      marvin.importStructure("cxsmiles", reloadedReaction);
      await sleep(500);
    } else if (selectedSketcher === "ketcher-select") {
      // wait for ketcher to load and if it has loaded then reload the scheme
      for (let i = 0; i < 5; i++) {
        let ketcher = getKetcher();
        if (ketcher !== undefined) {
          ketcher.setMolecule(reloadedReaction);
          break;
        } else {
          await sleep(250);
        }
      }
      // sleep used to allow ketcher to load smiles before proceeding
      await sleep(2000);
    }
  }
}

/**
 * Loads the reloaded RXN file to the active sketcher
 */
async function reloadSketcherFromRXN() {
  let $reactionRXN = $("#js-reaction-rxn");
  let reloadedReaction = $reactionRXN.val();
  if (reloadedReaction) {
    let selectedSketcher = $('input[name="sketcher-select"]:checked').attr(
      "id",
    );
    if (selectedSketcher === "marvin-select") {
      marvin.importStructure("rxn", reloadedReaction);
      await sleep(500);
    } else if (selectedSketcher === "ketcher-select") {
      let ketcher = getKetcher();
      ketcher.setMolecule(reloadedReaction);
      // sleep used to allow ketcher to load smiles before proceeding
      await sleep(2000);
    }
  }
}

/**
 * Locks the reaction if the current user is not the creator or if the reaction is saved as complete.
 */
function sketcherLockHandler() {
  if (ifCurrentUserIsNotCreator() || $("#js-complete").val() === "complete") {
    $("#marvin-sketcher").css("pointer-events", "none");
    $("#ketcher-sketcher").css("pointer-events", "none");
    $("#demo-button").prop("disabled", true);
    $("#action-button-submit").prop("disabled", true);
  }
}

/**
 * Sets up the autosave functionality for Ketcher and Marvin
 */
function setupSketcherAutosave() {
  if (getKetcher()) {
    setupKetcherAutosave();
  }
  if (marvin) {
    setupMarvinAutosave();
  }
}
