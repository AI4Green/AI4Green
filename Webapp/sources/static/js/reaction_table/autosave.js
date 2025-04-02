// run observer if not in tutorial mode
if ($("#js-tutorial").val() === "no" && $("#js-demo").val() !== "demo") {
  observer();
}

function observer() {
  // this line detects any changes user makes to any input field and saves 0.5 seconds after user stops focus on input
  $(document).on("change", ":input:not(#polymer-mode-select)", function (e) {
    setTimeout(autoSaveCheck(e), 500);
  });
  // on press remove solvent save
  $(document).on("click", ".js-remove-solvent", function (e) {
    setTimeout(autoSaveCheck(e), 500);
  });
  // on press remove reagent save
  $(document).on("click", ".js-remove-reagent", function (e) {
    setTimeout(autoSaveCheck(e), 500);
  });
}
/**
 * checks whether to autosave and if so what type of save to perform or returns blank if no save should be done
 *
 * @param {HTMLElement} [e=null] - the changed element that triggered the autosave function
 * @param {boolean} [sketcher=false] - true if the changed element is a chemical sketcher
 */
function autoSaveCheck(e = null, sketcher = false) {
  if (!checkIfSaveEnabled()) {
    return;
  }
  let listOfClasses = [];
  if (e !== null) {
    listOfClasses = e.currentTarget.classList.value;
  }
  let load_status = $("#js-load-status").val();
  let complete_status = $("#js-complete").val();
  // only autosave if: not loading, autosave field, not complete, and editor is creator
  if (
    !(
      load_status === "loading" ||
      listOfClasses.includes("no_autosave") ||
      complete_status === "complete" ||
      ifCurrentUserIsNotCreator()
    )
  ) {
    // use sketcher variable to either do the sketcher autosave function or the normal autosave function
    if (sketcher === false) {
      autoSave();
    } else if (sketcher === true) {
      // dont autosave sketcher if reaction table has loaded
      let reactionDiv = document.getElementById("reaction-table-div");
      sketcherAutoSave();
      // }
    }
  }
}

function autoSave() {
  // normal autosave
  postReactionData();
}

function completeReaction() {
  // from pressing the complete reaction button, it sends the "complete" argument to the save PostReactionData function
  let complete = "complete";
  postReactionData(complete);
}

function ifCurrentUserIsNotCreator() {
  let currentUser = $("#js-email").val();
  let creator = $("#js-creator-email").val();
  return creator !== currentUser;
}

async function sketcherAutoSave() {
  // autosave for when the sketcher has been updated
  let smiles;
  let rxn;
  let selectedSketcher = $('input[name="sketcher-select"]:checked').attr("id");
  if (selectedSketcher === "marvin-select") {
    smiles = await exportSmilesFromMarvin();
    rxn = await exportRXNFromMarvin(); // needed for polymer mode
  } else if (selectedSketcher === "ketcher-select") {
    smiles = await exportSmilesFromKetcher();
    rxn = await exportRXNFromKetcher();
  }
  if (smiles === ">>" || smiles === "" || smiles === undefined) {
    // catches if autosave occurs during a sketcher crash
    return;
  }
  // change here for reagent support
  let smilesNew = removeReagentsFromSmiles(smiles);
  $("#js-reaction-smiles").val(smilesNew);
  $("#js-reaction-rxn").val(rxn);
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let reactionID = $("#js-reaction-id").val();
  let polymerMode = $("#polymer-mode-select").is(":checked");
  let userEmail = "{{ current_user.email }}";
  updateReactionTable(smiles, workgroup, workbook, reactionID);
  $.ajax({
    url: "/_autosave_sketcher",
    type: "post",
    data: {
      workgroup: workgroup,
      workbook: workbook,
      reactionID: reactionID,
      userEmail: userEmail,
      reactionSmiles: smilesNew,
      reactionRXN: rxn,
      polymerMode: polymerMode,
    },
    dataType: "json",
    success: function () {
      flashUserSaveMessage("SKETCHER");
    },
    error: function () {
      flashUserErrorSavingMessage();
    },
  });
  updateReactionTable(smiles, workgroup, workbook, reactionID);
}

function updateReactionTable(smiles, workgroup, workbook, reactionID) {
  fetch("/autoupdate_reaction_table", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
    body: JSON.stringify({
      reaction_smiles: smiles,
      workgroup: workgroup,
      workbook: workbook,
      reaction_id: reactionID,
    }),
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      generateReactionTable(item);
    });
}

function postReactionData(complete = "not complete") {
  // this function posts data to the backend for saving and then updates the save indicator
  let {
    polymerIndices,
    polymerisationType,
    reactionID,
    workgroup,
    workbook,
    reactionSmiles,
    reactionRXN,
    userEmail,
    reactionName,
    reactionDescription,
    reactionImage,
    reactantPrimaryKeys,
    productPrimaryKeys,
    reagentPrimaryKeys,
    solventPrimaryKeys,
    numberOfReactants,
    limitingReactantTableNumber,
    reactantMasses,
    reactantMassesRaw,
    reactantAmounts,
    reactantAmountsRaw,
    reactantVolumes,
    reactantVolumesRaw,
    reactantEquivalents,
    reactantPhysicalForms,
    reactantDensities,
    reactantConcentrations,
    reagentNames,
    reagentDensities,
    reagentConcentrations,
    reagentMolecularWeights,
    reagentAmounts,
    reagentAmountsRaw,
    reagentEquivalents,
    reagentPhysicalForms,
    reagentHazards,
    reagentMasses,
    reagentMassesRaw,
    reagentVolumes,
    reagentVolumesRaw,
    reagentSmiles,
    solventPhysicalForms,
    solventNames,
    solventConcentrations,
    solventVolumes,
    solventHazards,
    productPhysicalForms,
    productAmounts,
    productAmountsRaw,
    productMasses,
    productMassesRaw,
    mainProductTableNumber,
    amountUnits,
    massUnits,
    volumeUnits,
    solventVolumeUnits,
    productAmountUnits,
    productMassUnits,
    // summary table elements
    realProductMass,
    unreactedReactantMass,
    polymerMn,
    polymerMw,
    polymerDispersity,
    polymerMassMethod,
    polymerMassCalibration,
    polymerTg,
    polymerTm,
    polymerTc,
    polymerThermalMethod,
    polymerThermalCalibration,
    reactionTemperature,
    elementSustainability,
    batchFlow,
    isolationMethod,
    catalystUsed,
    catalystRecovered,
    otherHazardTextArea,
    customProtocol1,
    customProtocol2,
    selectedRadioButtons,
    researcher,
    supervisor,
    // extra elements for data export
    reactantNames,
    reactantMolecularWeights,
    reactantHazards,
    reactantPhysicalFormsText,
    reagentPhysicalFormsText,
    solventPhysicalFormsText,
    productNames,
    productMolecularWeights,
    productHazards,
    productPhysicalFormsText,
    productIntendedDPs,
    massEfficiency,
    conversion,
    selectivity,
    toExport,
  } = getFieldData();
  let summary_to_print_el = document.getElementById("section-to-print");
  let summary_to_print = "no summary data";
  if (summary_to_print_el !== null) {
    summary_to_print = summary_to_print_el.outerHTML;
  }

  let $reactionSaveIndicator = $("#reaction-saved-indicator");
  $.ajax({
    url: "/_autosave",
    type: "post",
    data: {
      polymerIndices: polymerIndices,
      polymerisationType: polymerisationType,
      workgroup: workgroup,
      workbook: workbook,
      reactionID: reactionID,
      reactionSmiles: reactionSmiles,
      reactionRXN: reactionRXN,
      userEmail: userEmail,
      reactionName: reactionName,
      reactionImage: reactionImage,
      reactionDescription: reactionDescription,
      reactantPrimaryKeys: reactantPrimaryKeys,
      productPrimaryKeys: productPrimaryKeys,
      reagentPrimaryKeys: reagentPrimaryKeys,
      solventPrimaryKeys: solventPrimaryKeys,
      numberOfReactants: numberOfReactants,
      limitingReactantTableNumber: limitingReactantTableNumber,
      reactantMasses: reactantMasses,
      reactantMassesRaw: reactantMassesRaw,
      reactantAmounts: reactantAmounts,
      reactantAmountsRaw: reactantAmountsRaw,
      reactantVolumes: reactantVolumes,
      reactantVolumesRaw: reactantVolumesRaw,
      reactantEquivalents: reactantEquivalents,
      reactantPhysicalForms: reactantPhysicalForms,
      reactantDensities: reactantDensities,
      reactantConcentrations: reactantConcentrations,
      reagentNames: reagentNames,
      reagentDensities: reagentDensities,
      reagentConcentrations: reagentConcentrations,
      reagentMolecularWeights: reagentMolecularWeights,
      reagentAmounts: reagentAmounts,
      reagentAmountsRaw: reagentAmountsRaw,
      reagentEquivalents: reagentEquivalents,
      reagentPhysicalForms: reagentPhysicalForms,
      reagentHazards: reagentHazards,
      reagentMasses: reagentMasses,
      reagentMassesRaw: reagentMassesRaw,
      reagentVolumes: reagentVolumes,
      reagentVolumesRaw: reagentVolumesRaw,
      reagentSmiles: reagentSmiles,
      solventPhysicalForms: solventPhysicalForms,
      solventNames: solventNames,
      solventConcentrations: solventConcentrations,
      solventVolumes: solventVolumes,
      solventHazards: solventHazards,
      productPhysicalForms: productPhysicalForms,
      productAmounts: productAmounts,
      productAmountsRaw: productAmountsRaw,
      productMasses: productMasses,
      productMassesRaw: productMassesRaw,
      mainProductTableNumber: mainProductTableNumber,
      amountUnits: amountUnits,
      massUnits: massUnits,
      volumeUnits: volumeUnits,
      solventVolumeUnits: solventVolumeUnits,
      productAmountUnits: productAmountUnits,
      productMassUnits: productMassUnits,
      realProductMass: realProductMass,
      unreactedReactantMass: unreactedReactantMass,
      polymerMn: polymerMn,
      polymerMw: polymerMw,
      polymerDispersity: polymerDispersity,
      polymerMassMethod: polymerMassMethod,
      polymerMassCalibration: polymerMassCalibration,
      polymerTg: polymerTg,
      polymerTm: polymerTm,
      polymerTc: polymerTc,
      polymerThermalMethod: polymerThermalMethod,
      polymerThermalCalibration: polymerThermalCalibration,
      reactionTemperature: reactionTemperature,
      elementSustainability: elementSustainability,
      batchFlow: batchFlow,
      isolationMethod: isolationMethod,
      catalystUsed: catalystUsed,
      catalystRecovered: catalystRecovered,
      otherHazardTextArea: otherHazardTextArea,
      customProtocol1: customProtocol1,
      customProtocol2: customProtocol2,
      selectedRadioButtons: selectedRadioButtons,
      researcher: researcher,
      supervisor: supervisor,
      complete: complete,
      reactantNames: reactantNames,
      reactantMolecularWeights: reactantMolecularWeights,
      reactantHazards: reactantHazards,
      reactantPhysicalFormsText: reactantPhysicalFormsText,
      reagentPhysicalFormsText: reagentPhysicalFormsText,
      solventPhysicalFormsText: solventPhysicalFormsText,
      productNames: productNames,
      productMolecularWeights: productMolecularWeights,
      productHazards: productHazards,
      productPhysicalFormsText: productPhysicalFormsText,
      productIntendedDPs: productIntendedDPs,
      summary_to_print: summary_to_print,
      massEfficiency: massEfficiency,
      conversion: conversion,
      selectivity: selectivity,
      toExport: JSON.stringify(toExport),
    },
    dataType: "json",
    success: async function (response) {
      if (complete === "complete") {
        if (
          response.feedback ===
          "Please enter unreacted and product mass to mark reaction as complete!"
        ) {
          $reactionSaveIndicator.text(
            "Enter unreacted and product mass to mark reaction as complete",
          );
          $reactionSaveIndicator
            .removeClass()
            .addClass("reaction-save-failure")
            .fadeIn("fast");
          document
            .querySelector("#scrollToUnreacted")
            .scrollIntoView({ behavior: "smooth" });
          // if the user is trying to complete reaction, highlight red empty fields that are compulsory.
          autoChangeRequiredStyling2("#js-unreacted-reactant-mass");
          autoChangeRequiredStyling2("#js-real-product-mass");
        } else if (
          response.feedback ===
          "Please fill in all required data to mark reaction as complete!"
        ) {
          $reactionSaveIndicator.text(
            "Please fill in all required data to mark reaction as complete",
          );
          $reactionSaveIndicator
            .removeClass()
            .addClass("reaction-save-failure")
            .fadeIn("fast");
          document
            .querySelector("#js-reaction-description")
            .scrollIntoView({ behavior: "smooth" });
        } else if (response.feedback === "Reaction locked") {
          $reactionSaveIndicator.text("Reaction Changes Saved & Locked");
          $reactionSaveIndicator
            .removeClass()
            .addClass("reaction-save-success")
            .fadeIn("fast");
          autoChangeRequiredStyling2("#js-real-product-mass");
          autoChangeRequiredStyling2("#js-unreacted-reactant-mass");
          $("#js-real-product-mass")
            .removeClass("add-highlight-unfilled-cell")
            .addClass("remove-highlight-summary-cell");
          $("#js-unreacted-reactant-mass")
            .removeClass("add-highlight-unfilled-cell")
            .addClass("remove-highlight-summary-cell");
          $("#js-complete").val("complete");
          controlLockedReactionFunctionality();
          showLoadingOverlay("Locking Reaction");
          displayOverlayWhilstMakingPDF("locked");
        }
      } else {
        // if not locking reaction save as normal
        flashUserSaveMessage("AUTOSAVE");
      }
    },
    error: function () {
      flashUserErrorSavingMessage();
    },
  });
}

/**
 * Disables page contents and buttons outside of page contents div
 * Enables specific buttons for locked reactions.
 * Upload, view, delete files. Print summary. Reaction notes.
 */
function controlLockedReactionFunctionality() {
  // disabling functionality
  $("#page-contents :input").prop("disabled", true);
  $("#complete-reaction-button").prop("disabled", true);
  // restoring specific functionality
  $("#print-pdf").prop("disabled", false);
  $("#reaction-note-button").show().prop("disabled", false);
  $("#new-reaction-note-modal").find("*").prop("disabled", false);
  $("#file-list").find("*").prop("disabled", false);
  // only the creator can upload or delete files or add reaction note comments
  if (ifCurrentUserIsNotCreator()) {
    $(".delete-file-button").prop("disabled", true);
    $("#file-upload-div").find("*").prop("disabled", true);
    $("#reaction-note-button").hide();
  } else {
    $("#file-upload-div").find("*").prop("disabled", false);
  }
}

function getFieldData() {
  let polymerisationType = $("#js-polymerisation-type").val();
  if (!polymerisationType) {
    polymerisationType = "";
  }
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let reactionID = $("#js-reaction-id").val();
  let userEmail = "{{ current_user.email }}";
  // variables in the reaction table stage
  // reactant data
  let reactionSmiles = $("#js-reaction-smiles").val();
  let reactionRXN = $("#js-reaction-rxn").val();
  let polymerIndices = JSON.stringify(identifyPolymers(reactionRXN));
  let reactionName = $("#js-reaction-name").val();
  let reactionDescription = $("#js-reaction-description").val();
  let reactionImage = $("#image").attr("src");
  if (!reactionImage) {
    reactionImage = "";
  }
  let limitingReactantTableNumber = $(
    "#js-limiting-reactant-table-number",
  ).val();
  let numberOfReactants = Number($("#js-number-of-reactants").val());
  let reactantMassID, reactantEquivalentID, reactantPhysicalFormID;
  let reactantNames = "";
  let reactantPrimaryKeys = "";
  let reactantMasses = "";
  let reactantMassesRaw = "";
  let reactantAmounts = "";
  let reactantAmountsRaw = "";
  let reactantVolumes = "";
  let reactantVolumesRaw = "";
  let reactantEquivalents = "";
  let reactantPhysicalForms = "";
  let reactantDensities = "";
  let reactantConcentrations = "";
  let reactantMolecularWeights = "";
  let reactantHazards = "";
  let reactantPhysicalFormsText = "";

  for (let i = 1; i <= numberOfReactants; i++) {
    let reactantNamesID = "#js-reactant" + i;
    reactantNames += $(reactantNamesID).val() + ";";
    let reactantPrimaryKeyID = "#js-reactant-primary-key" + i;
    reactantPrimaryKeys += $(reactantPrimaryKeyID).val() + ";";
    reactantMassID = "#js-reactant-rounded-mass" + i;
    reactantMasses += $(reactantMassID).val() + ";";
    let reactantMassRawID = "#js-reactant-mass" + i;
    reactantMassesRaw += $(reactantMassRawID).val() + ";";
    let reactantAmountID = "#js-reactant-rounded-amount" + i;
    reactantAmounts += $(reactantAmountID).val() + ";";
    let reactantAmountRawID = "#js-reactant-amount" + i;
    reactantAmountsRaw += $(reactantAmountRawID).val() + ";";
    let reactantVolumeRawID = "#js-reactant-volume" + i;
    reactantVolumesRaw = $(reactantVolumeRawID).val() + ";";
    let reactantVolumeID = "#js-reactant-rounded-volume" + i;
    reactantVolumes += $(reactantVolumeID).val() + ";";
    let reactantDensityID = "#js-reactant-density" + i;
    reactantDensities += $(reactantDensityID).val() + ";";
    let reactantConcentrationID = "#js-reactant-concentration" + i;
    reactantConcentrations += $(reactantConcentrationID).val() + ";";
    reactantEquivalentID = "#js-reactant-equivalent" + i;
    reactantEquivalents += $(reactantEquivalentID).val() + ";";
    reactantPhysicalFormID = "#js-reactant-physical-form" + i;
    let reactantPhysicalForm = $(reactantPhysicalFormID).prop("selectedIndex");
    reactantPhysicalForms += reactantPhysicalForm + ";";
    let reactantMolecularWeightsID = "#js-reactant-molecular-weight" + i;
    reactantMolecularWeights += $(reactantMolecularWeightsID).val() + ";";
    let reactantHazardsID = "#js-reactant-hazards" + i;
    reactantHazards += $(reactantHazardsID).val() + ";";
    reactantPhysicalFormsText += $(reactantPhysicalFormID).val() + ";";
  }
  // reagent data
  let reagentPrimaryKeys = "";
  let numberOfReagents = Number($("#js-number-of-reagents").val());
  let reagentEquivalentID;
  let reagentEquivalents = "";
  let reagentPhysicalForms = "";
  let reagentNames = "";
  let reagentMolecularWeights = "";
  let reagentDensities = "";
  let reagentConcentrations = "";
  let reagentAmounts = "";
  let reagentAmountsRaw = "";
  let reagentHazards = "";
  let reagentMasses = "";
  let reagentMassesRaw = "";
  let reagentVolumes = "";
  let reagentVolumesRaw = "";
  let reagentPhysicalFormsText = "";
  let reagentSmiles = "";
  for (let i = 1; i <= numberOfReagents; i++) {
    let reagentPrimaryKeyID = "#js-reagent-primary-key" + i;
    reagentPrimaryKeys += $(reagentPrimaryKeyID).val() + ";";
    reagentEquivalentID = "#js-reagent-equivalent" + i;
    reagentEquivalents += $(reagentEquivalentID).val() + ";";
    let reagentPhysicalFormID = "#js-reagent-physical-form" + i;
    let reagentPhysicalForm = $(reagentPhysicalFormID).prop("selectedIndex");
    reagentPhysicalForms += reagentPhysicalForm + ";";
    let reagentNameID = "#js-reagent" + i;
    reagentNames += $(reagentNameID).val() + ";";
    let reagentMolecularWeightID = "#js-reagent-molecular-weight" + i;
    reagentMolecularWeights += $(reagentMolecularWeightID).val() + ";";
    let reagentDensityID = "#js-reagent-density" + i;
    reagentDensities += $(reagentDensityID).val() + ";";
    let reagentConcentrationID = "#js-reagent-concentration" + i;
    reagentConcentrations += $(reagentConcentrationID).val() + ";";
    let reagentAmountID = "#js-reagent-rounded-amount" + i;
    reagentAmounts += $(reagentAmountID).val() + ";";
    let reagentAmountRawID = "#js-reagent-amount" + i;
    reagentAmountsRaw += $(reagentAmountRawID).val() + ";";
    let reagentHazardID = "#js-reagent-hazards" + i;
    reagentHazards += $(reagentHazardID).val() + ";";
    let reagentMassesID = "#js-reagent-rounded-mass" + i;
    reagentMasses += $(reagentMassesID).val() + ";";
    let reagentMassesRawID = "#js-reagent-mass" + i;
    reagentMassesRaw += $(reagentMassesRawID).val() + ";";
    let reagentVolumeID = "#js-reagent-rounded-volume" + i;
    reagentVolumes += $(reagentVolumeID).val() + ";";
    let reagentVolumeRawID = "#js-reagent-volume" + i;
    reagentVolumesRaw += $(reagentVolumeRawID).val() + ";";
    reagentPhysicalFormsText += $(reagentPhysicalFormID).val() + ";";
    let reagentSmilesID = "#js-reagent-smiles" + i;
    reagentSmiles += $(reagentSmilesID).val() + ";";
  }
  // solvent data
  let solventPrimaryKeys = "";
  let numberOfSolvents = Number($("#js-number-of-solvents").val());
  let solventVolumeID, solventPhysicalFormID;
  let solventNames = "";
  let solventVolumes = "";
  let solventConcentrations = "";
  let solventPhysicalForms = "";
  let solventHazards = "";
  let solventPhysicalFormsText = "";
  for (let i = 1; i <= numberOfSolvents; i++) {
    let solventPrimaryKeyID = "#js-solvent-primary-key" + i;
    solventPrimaryKeys += $(solventPrimaryKeyID).val() + ";";
    let solventNameID = "#js-solvent" + i;
    solventNames += $(solventNameID).val() + ";";
    solventVolumeID = "#js-solvent-volume" + i;
    solventVolumes += $(solventVolumeID).val() + ";";
    let solventConcentrationID = "#js-solvent-rounded-concentration" + i;
    solventConcentrations += $(solventConcentrationID).val() + ";";
    solventPhysicalFormID = "#js-solvent-physical-form" + i;
    let solventPhysicalForm = $(solventPhysicalFormID).prop("selectedIndex");
    solventPhysicalForms += solventPhysicalForm + ";";
    let solventHazardID = "#js-solvent-hazards" + i;
    solventHazards += $(solventHazardID).val() + ";";
    solventPhysicalFormsText += $(solventPhysicalFormID).val() + ";";
  }
  // product data
  let productPrimaryKeys = "";
  let productNames = "";
  let productMolecularWeights = "";
  let productHazards = "";
  let mainProductTableNumber = $("#js-main-product-table-number").val();
  let numberOfProducts = Number($("#js-number-of-products").val());
  let productPhysicalFormID;
  let productPhysicalForms = "";
  let productPhysicalFormsText = "";
  let productIntendedDPs = "";
  let productMassesRaw = "";
  let productAmounts = "";
  let productAmountsRaw = "";
  let productMasses = "";
  for (let i = 1; i <= numberOfProducts; i++) {
    let productPrimaryKeyID = "#js-product-primary-key" + i;
    productPrimaryKeys += $(productPrimaryKeyID).val() + ";";
    let productNameID = "#js-product" + i;
    productNames += $(productNameID).val() + ";";
    let productMolecularWeightID = "#js-product-molecular-weight" + i;
    productMolecularWeights += $(productMolecularWeightID).val() + ";";
    let productHazardsID = "#js-product-hazard" + i;
    productHazards += $(productHazardsID).val() + ";";
    productPhysicalFormID = "#js-product-physical-form" + i;
    let productPhysicalForm = $(productPhysicalFormID).prop("selectedIndex");
    productPhysicalForms += productPhysicalForm + ";";
    productPhysicalFormsText += $(productPhysicalFormID).val() + ";";
    let productIntendedDPID = "#js-product-intended-dp" + i;
    productIntendedDPs += $(productIntendedDPID).val() + ";";
    let productMassFormID = "#js-product-rounded-mass" + i;
    productMasses += $(productMassFormID).val() + ";";
    let productMassRawFormID = "#js-product-mass" + i;
    productMassesRaw += $(productMassRawFormID).val() + ";";
    let productAmountFormID = "#js-product-rounded-amount" + i;
    productAmounts += $(productAmountFormID).val() + ";";
    let productAmountRawID = "#js-product-amount" + i;
    productAmountsRaw += $(productAmountRawID).val() + ";";
  }
  // unit data
  let amountUnits = unitData("#js-amount-unit", "mmol");
  let volumeUnits = unitData("#js-volume-unit", "mL");
  let massUnits = unitData("#js-mass-unit", "mg");
  let solventVolumeUnits = unitData("#js-solvent-volume-unit", "mL");
  let productAmountUnits = unitData("#js-product-amount-unit", "mmol");
  let productMassUnits = unitData("#js-product-mass-unit", "mg");
  // summary table fields
  // yield data
  let summaryTableData = JSON.parse($("#js-summary-table-data").val());
  let realProductMass = summaryInput("#js-real-product-mass", "");
  let unreactedReactantMass = summaryInput("#js-unreacted-reactant-mass", "");
  let polymerMn = summaryInput("#js-polymer-mn", "");
  let polymerMw = summaryInput("#js-polymer-mw", "");
  let polymerDispersity = summaryInput("#js-polymer-dispersity", "");
  let polymerMassMethod = summaryInput("#js-polymer-mass-method", "-select-");
  let polymerMassCalibration = summaryInput("#js-polymer-mass-calibration", "");
  let polymerTg = summaryInput("#js-polymer-tg", "");
  let polymerTm = summaryInput("#js-polymer-tm", "");
  let polymerTc = summaryInput("#js-polymer-tc", "");
  let polymerThermalMethod = summaryInput(
    "#js-polymer-thermal-method",
    "-select-",
  );
  let polymerThermalCalibration = summaryInput(
    "#js-polymer-thermal-calibration",
    "",
  );
  let reactionTemperature = summaryInput("#js-temperature", "");
  let elementSustainability = summaryDropDown(
    "#js-elements",
    summaryTableData,
    "element_sustainability",
  );
  let batchFlow = summaryInput("#js-batch-flow", "-select-");
  let isolationMethod = summaryDropDown(
    "#js-isolation",
    summaryTableData,
    "isolation_method",
  );
  let conversion = summaryInput("#js-conversion", "");
  let selectivity = summaryInput("#js-selectivity", "");
  let catalystUsed = summaryInput("#js-catalyst", "-select-");
  let catalystRecovered = summaryInput("#js-recovery", "-select-");
  let otherHazardTextArea = summaryInput("#other-risks-textbox", "");
  let customProtocol1 = summaryInput("#field1-text", "");
  let customProtocol2 = summaryInput("#field2-text", "");
  let researcher = summaryInput("#js-researcher", "");
  let supervisor = summaryInput("#js-supervisor", "");
  // standard protocol radio buttons
  let standard_protocol_radio_buttons = [
    "#cyanide",
    "#hplc",
    "#massSpec",
    "#pyrophorics",
    "#microwave",
    "#diazotisation",
    "#hydrogenation",
    "#peptideSynthesis",
    "#ozone",
    "#freeRadicals",
    "#liquidAmmonia",
    "#peroxides",
    "#sealedTube",
    "#field1",
    "#field2",
    "#nonHalogenatedSolvent",
    "#halogenatedSolvent",
    "#specialistContainer",
    "#sinkWithExcessWater",
    "#standardSpillResponse",
    "#otherSpillResponse",
    "#slight",
    "#serious",
    "#major",
    "#lowLikelihood",
    "#possible",
    "#frequentOccur",
    "#individual",
    "#localLabs",
    "#buildingWide",
  ];
  let selectedRadioButtons = "";
  standard_protocol_radio_buttons.filter(function (currentElement) {
    if ($(currentElement).is(":checked")) {
      currentElement = currentElement.substring(1);
      selectedRadioButtons += currentElement + ";";
    }
  });
  let mass_efficiency = summaryInput("#js-me", "");
  let toExport = [];
  $(".to-export").each(function () {
    toExport.push({ key: this.name, value: $(this).val() });
  });
  $(".to-export-risk").each(function () {
    if ($(this).is(":checked")) {
      toExport.push({ key: this.name, value: $(this).val() });
    }
  });
  $(".to-export-checkbox").each(function () {
    if ($(this).is(":checked")) {
      toExport.push({ key: this.name, value: "Yes" });
    }
  });
  $(".to-export-hazard-matrix").each(function () {
    if (this.innerText !== "") {
      toExport.push({ key: this.id, value: this.innerText });
    }
  });
  const summary_sustainable_metrics = {
    "#js-temperature": "Temperature Sustainability",
    "#js-elements": "Elements Sustainability",
    "#js-batch-flow": "Batch or Flow Sustainability",
    "#js-isolation": "Isolation Sustainability",
    "#js-catalyst": "Catalyst Sustainability",
    "#js-recovery": "Recovery Sustainability",
    "#js-ae-cell": "Atom Economy Sustainability",
    "#js-me": "Mass Efficiency Sustainability",
    "#js-yield": "Yield Sustainability",
    "#js-conversion": "Conversion Sustainability",
    "#js-selectivity": "Selectivity Sustainability",
  };
  for (const property in summary_sustainable_metrics) {
    const colour = $(`${property}`).attr("class");
    if (colour !== "to-export") {
      toExport.push({
        key: `${summary_sustainable_metrics[property]}`,
        value: colour,
      });
    }
  }
  return {
    polymerIndices: polymerIndices,
    polymerisationType: polymerisationType,
    reactionID: reactionID,
    workgroup: workgroup,
    workbook: workbook,
    reactionSmiles: reactionSmiles,
    reactionRXN: reactionRXN,
    userEmail: userEmail,
    reactionName: reactionName,
    reactionDescription: reactionDescription,
    reactionImage: reactionImage,
    reactantPrimaryKeys: reactantPrimaryKeys,
    productPrimaryKeys: productPrimaryKeys,
    reagentPrimaryKeys: reagentPrimaryKeys,
    solventPrimaryKeys: solventPrimaryKeys,
    numberOfReactants: numberOfReactants,
    limitingReactantTableNumber: limitingReactantTableNumber,
    reactantMasses: reactantMasses,
    reactantMassesRaw: reactantMassesRaw,
    reactantAmounts: reactantAmounts,
    reactantAmountsRaw: reactantAmountsRaw,
    reactantVolumes: reactantVolumes,
    reactantVolumesRaw: reactantVolumesRaw,
    reactantEquivalents: reactantEquivalents,
    reactantPhysicalForms: reactantPhysicalForms,
    reactantDensities: reactantDensities,
    reactantConcentrations: reactantConcentrations,
    reagentNames: reagentNames,
    reagentDensities: reagentDensities,
    reagentConcentrations: reagentConcentrations,
    reagentMolecularWeights: reagentMolecularWeights,
    reagentAmounts: reagentAmounts,
    reagentAmountsRaw: reagentAmountsRaw,
    reagentEquivalents: reagentEquivalents,
    reagentPhysicalForms: reagentPhysicalForms,
    reagentHazards: reagentHazards,
    reagentMasses: reagentMasses,
    reagentMassesRaw: reagentMassesRaw,
    reagentVolumes: reagentVolumes,
    reagentVolumesRaw: reagentVolumesRaw,
    reagentSmiles: reagentSmiles,
    solventPhysicalForms: solventPhysicalForms,
    solventNames: solventNames,
    solventConcentrations: solventConcentrations,
    solventVolumes: solventVolumes,
    solventHazards: solventHazards,
    productPhysicalForms: productPhysicalForms,
    productAmounts: productAmounts,
    productAmountsRaw: productAmountsRaw,
    productMasses: productMasses,
    productMassesRaw: productMassesRaw,
    mainProductTableNumber: mainProductTableNumber,
    amountUnits: amountUnits,
    massUnits: massUnits,
    volumeUnits: volumeUnits,
    solventVolumeUnits: solventVolumeUnits,
    productAmountUnits: productAmountUnits,
    productMassUnits: productMassUnits,
    realProductMass: realProductMass,
    unreactedReactantMass: unreactedReactantMass,
    polymerMn: polymerMn,
    polymerMw: polymerMw,
    polymerDispersity: polymerDispersity,
    polymerMassMethod: polymerMassMethod,
    polymerMassCalibration: polymerMassCalibration,
    polymerTg: polymerTg,
    polymerTm: polymerTm,
    polymerTc: polymerTc,
    polymerThermalMethod: polymerThermalMethod,
    polymerThermalCalibration: polymerThermalCalibration,
    reactionTemperature: reactionTemperature,
    elementSustainability: elementSustainability,
    batchFlow: batchFlow,
    isolationMethod: isolationMethod,
    catalystUsed: catalystUsed,
    catalystRecovered: catalystRecovered,
    otherHazardTextArea: otherHazardTextArea,
    customProtocol1: customProtocol1,
    customProtocol2: customProtocol2,
    selectedRadioButtons: selectedRadioButtons,
    researcher: researcher,
    supervisor: supervisor,
    reactantNames: reactantNames,
    reactantMolecularWeights: reactantMolecularWeights,
    reactantHazards: reactantHazards,
    reactantPhysicalFormsText: reactantPhysicalFormsText,
    reagentPhysicalFormsText: reagentPhysicalFormsText,
    solventPhysicalFormsText: solventPhysicalFormsText,
    productNames: productNames,
    productMolecularWeights: productMolecularWeights,
    productHazards: productHazards,
    productPhysicalFormsText: productPhysicalFormsText,
    productIntendedDPs: productIntendedDPs,
    massEfficiency: mass_efficiency,
    conversion: conversion,
    selectivity: selectivity,
    toExport: toExport,
  };
}

function unitData(unit_id, default_unit) {
  // returns unit data but changes to a default value if null
  let savedUnit = $(unit_id).val();
  if (savedUnit === null) {
    savedUnit = default_unit;
  }
  return savedUnit;
}

function summaryInput(field_id, default_value) {
  // if undefined then assign the default value
  let summaryInputField = $(field_id).val();
  if (summaryInputField === undefined) {
    summaryInputField = default_value;
  }
  return summaryInputField;
}

function summaryDropDown(field_id, summary_json, json_key) {
  // if dropdown filled - use that value. if dropdown not present but saved value is - use that. if neither, save as undefined
  let summaryDropDownField = $(field_id).prop("selectedIndex");
  if (
    summaryDropDownField === undefined &&
    summary_json[json_key] === "undefined"
  ) {
    summaryDropDownField = "undefined";
  } else if (
    summaryDropDownField === undefined &&
    summary_json[json_key] !== "undefined"
  ) {
    summaryDropDownField = summary_json[json_key];
  }
  return summaryDropDownField;
}

function fade_save_message() {
  // fades out the save message
  $("#reaction-saved-indicator").fadeOut("slow");
}

function flashUserSaveMessage(message) {
  let $reactionSaveIndicator = $("#reaction-saved-indicator");
  $reactionSaveIndicator.text(message);
  $reactionSaveIndicator
    .removeClass()
    .addClass("reaction-save-success")
    .fadeIn("fast");
  setTimeout(fade_save_message, 1000);
}

function flashUserErrorSavingMessage() {
  // shows error message to user
  let $reactionSaveIndicator = $("#reaction-saved-indicator");
  $reactionSaveIndicator.text("Failed To Save Reaction Changes");
  $reactionSaveIndicator
    .removeClass()
    .addClass("reaction-save-failure")
    .fadeIn("fast");
}
