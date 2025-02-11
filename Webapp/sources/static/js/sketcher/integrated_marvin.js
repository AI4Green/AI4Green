let marvin;
/**
 * Assigns the global variable 'marvin' to the marvin sketcher object and prepares the sketcher ready for use
 */
async function setupNewMarvinSketcher() {
  let marvinKey = $("#marvin-js-key").val();
  marvin = await newMarvinSketcher("#marvin-sketcher", marvinKey);
  if (marvin) {
    marvin.setDisplaySettings({ toolbars: "reaction" });
    $("#marvin-sketcher").hide();
  }
}

/**
 * Loads the required scripts and returns the sketcher or aborts if scripts fail to load
 * @param {string} IDSelector - The ID of the new DOM element in the HTML
 * @param {string} marvinKey - The license key tied to the URL to access Marvin JS services
 * @return {Object} MarvinJS Chemical Editor object with access to the MarvinJS methods
 */
async function newMarvinSketcher(IDSelector, marvinKey) {
  let marvinKey1 = await getMarvinKey();

  try {
    // Use Promise.race to set a timeout for the entire operation
    const timeoutPromise = new Promise((_, reject) => {
      setTimeout(() => reject(new Error("Timeout")), 5000);
    });

    // Use Promise.all to wait for both scripts to load or for the timeout
    await Promise.all([
      $.getScript("https://marvinjs.chemicalize.com/v1/client.js"),
      $.getScript(
        `https://marvinjs.chemicalize.com/v1/${marvinKey1}/client-settings.js`,
      ),
      timeoutPromise, // Timeout promise
    ]);
  } catch (error) {
    // Check if the error is due to a timeout
    if (error.message === "Timeout") {
      console.log("Operation timed out");
    } else {
      console.log(error.status);
      abortMarvinSketcherCreation();
    }
    return null;
  }
  $("#marvin-select").prop("disabled", false);
  return ChemicalizeMarvinJs.createEditor(IDSelector);
}

/**
 * Gets the marvin key from the application backend
 * @return {Promise<string>} - The marvin JS license key for accessing their API services
 */
async function getMarvinKey() {
  let marvinJsKey;
  await $.ajax({
    method: "POST",
    url: "/get_marvinjs_key",
  }).then(function (response) {
    marvinJsKey = response.marvinjs_key;
  });
  return marvinJsKey;
}

/**
 * Shows the Marvin Sketcher and hides Ketcher and exports SMILES from ketcher to Marvin
 * @return {Promise<void>}
 */
async function switchToMarvinSketcher() {
  $("#ketcher-sketcher").hide();
  if (marvin) {
    await exportReactionFromKetcherToMarvin();
    $("#marvin-sketcher").show();
  } else {
    alert("Marvin JS is temporarily unavailable");
  }
}

/**
 * Switches back to Ketcher, disables further switching and hides the loading circle again
 */
function abortMarvinSketcherCreation() {
  // runs when
  flashMarvinDownMessage();
  $("#ketcher-select").prop("checked", true);
  $("#marvin-select").prop("disabled", true);
}

function flashMarvinDownMessage() {
  let $userMessageElement = $("#reaction-saved-indicator");
  $userMessageElement.text("Marvin JS unavailable");
  $userMessageElement
    .removeClass()
    .addClass("reaction-save-success")
    .fadeIn("fast");
  setTimeout(fade_save_message, 3000);
}

/**
 * Autosave when the sketcher is edited
 */
function setupMarvinAutosave() {
  marvin.on("molchange", function () {
    // sketcher arg to true to save just smiles
    setTimeout(autoSaveCheck(null, true), 100);
  });
}

/**
 * Exports structures from ketcher to marvin via SMILES (or RXN in polymer mode)
 * @return {Promise<void>}
 */
async function exportReactionFromKetcherToMarvin() {
  let reaction;
  if (getKetcher()) {
    const polymerMode = await getPolymerMode();
    if (polymerMode === true) {
      reaction = await exportRXNFromKetcher();
      if (reaction) {
        marvin.importStructure("rxn", reaction);
      }
    } else {
      reaction = await exportSmilesFromKetcher();
      if (reaction) {
        marvin.importStructure("cxsmiles", reaction);
      }
    }
  }
}

/**
 * Exports the current contents of Marvin JS as SMILES
 * @return {Promise<string>} Reaction SMILES string where reagents are optional and later discarded
 * in format: "reactant1.reactantX>reagent1.reagentX>product1.productX"
 */
function exportSmilesFromMarvin() {
  return marvin
    .exportStructure("cxsmiles", { extra: "f" })
    .then(function (smiles) {
      return smiles;
    });
}

/**
 * Exports the current contents of Marvin JS as RXN
 * @return {Promise<string>} RXN file
 */
function exportRXNFromMarvin() {
  return marvin
    .exportStructure("rxn", { multiplesgroup: false })
    .then(function (reaction) {
      return reaction;
    });
}

/**
 * Exports the current contents of Marvin JS as MOL
 * @return {Promise<string>} MOL file
 */
function exportMOLFromMarvin() {
  return marvin
    .exportStructure("mol", { multiplesgroup: false })
    .then(function (reaction) {
      return reaction;
    });
}

/**
 * Enters the example reaction into the Marvin editor
 */
function marvinExampleSmiles() {
  let smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1";

  getPolymerMode().then((polymerMode) => {
    if (polymerMode === true) {
      marvin.importStructure("rxn", getExamplePolymer());
    } else {
      marvin.importStructure("cxsmiles", smiles);
    }
  });
}

/**
 * Makes an image of the current reaction scheme and saves to a hidden HTML input
 */
async function exportMarvinImage() {
  let settings = { width: 600, height: 400, multiplesgroup: true };
  let reactionSchemeImage;
  await marvin.exportStructure("jpeg", settings).then(function (source) {
    reactionSchemeImage = source;
  });
  return reactionSchemeImage;
}
