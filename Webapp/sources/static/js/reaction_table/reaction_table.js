// Functions for the reaction table
let reactionTable, amountFactor, massFactor, volumeFactor;
function initialiseReactionTable() {
  // called from integrated marvin js
  // runs all the functions required to initialise the reaction table page
  [reactionTable, amountFactor, massFactor, volumeFactor] =
    commonGlobalVariables();
  // gets the page ready by checking first reactant as limiting and adjusting if in demo mode
  $('input[name="reactant-limiting"]').first().prop("checked", true);
  hideButtonsInDemoMode();
  // reload before autofill functions
  // The next functions make changes based on the limiting reactant/main product
  updateTableAfterLimitingReactantChange();
  updateMainProduct();

  // set up event listeners to change table amounts/volumes/masses
  autoFillLimitingReactant();
  autofillReactantFields2();
  autofillProductFields2();
  initialiseReagentListeners();
  initialiseSolventListeners();
  updateStyling();
  setColours();
  $("#js-load-status").on("change", function () {
    autoSaveCheck(null);
  });
}

/**
 * Hides specific buttons when the application is in demo mode.
 *
 * This function checks if the element with the ID "js-demo" has a value of "demo".
 * If it does, it hides the elements with the following IDs:
 * - "reaction-name-description"
 * - "js-add-new-reagent-by-table"
 * - "js-add-new-solvent-by-table"
 */
function hideButtonsInDemoMode() {
  if (getVal($("#js-demo")) === "demo") {
    document.getElementById("reaction-name-description").style.display = "none";
    document.getElementById("js-add-new-reagent-by-table").style.display =
      "none";
    document.getElementById("js-add-new-solvent-by-table").style.display =
      "none";
  }
}

/**
 * Event listener for limiting reactant change.
 * Updates the reaction table and related fields when the limiting reactant selection changes.
 *
 * This function listens for clicks on elements with the class "js-reactant-limiting".
 * When triggered, it:
 * - Updates the limiting reactant table number in the HTML.
 * - Calls `updateReactantsAfterLimitingReactantChange()` to refresh reactant data.
 * - Sends an AJAX request to autofill reactant fields and updates required styling.
 * - Calls functions to autofill solvent and product fields.
 */
function updateTableAfterLimitingReactantChange() {
  $(".js-reactant-limiting").click(function () {
    // update the limiting reactant table number in the html
    let limitingReactantTableNumber = getLimitingReactantTableNumber();
    $("#js-limiting-reactant-table-number").val(limitingReactantTableNumber);
    updateReactantsAfterLimitingReactantChange();
    // when the limiting reactant is changed, the autofill functions need to be reran
    // styling updating after limiting reactant is changed
    $.ajax({
      url: autofillReactantFields2(),
      success: function () {
        updateRequiredStylingLimited();
      },
    });
    autofillSolventFields2();
    autofillProductFields2();
  });
}

/**
 * Updates reactant fields when the limiting reactant changes.
 *
 * This function iterates through all reactants in the reaction table and:
 * - Recalculates equivalents relative to the new limiting reactant.
 * - Updates the "readonly" and "editable" states of equivalent and mass fields accordingly.
 *
 * If a reactant is marked as limiting, its equivalents are recalculated.
 * Otherwise:
 * - The equivalent field becomes editable.
 * - The rounded mass field becomes read-only.
 */
function updateReactantsAfterLimitingReactantChange() {
  for (let i = 1; i < reactionTable.numberOfReactants + 1; i++) {
    if ($("#js-reactant-limiting" + i).is(":checked")) {
      // whichever reactant is limiting - all equivs are recalculated relative to this
      updateEquivalentsRelativeToNewLimitingReactant(i);
    } else {
      $("#js-reactant-equivalent" + i)
        .attr("readonly", false)
        .removeClass("readonly-cell")
        .addClass("editable-cell");
      $("#js-reactant-rounded-mass" + i)
        .attr("readonly", true)
        .removeClass("editable-cell")
        .addClass("readonly-cell");
    }
  }
}

/**
 * Updates the equivalents of all reactants and reagents when a new limiting reactant is selected.
 *
 * @param {number} x - The index of the newly selected limiting reactant.
 *
 * This function:
 * - Calculates a change factor to normalize the new limiting reactant's equivalents to 1.
 * - Scales the equivalents of all other reactants accordingly.
 * - Updates reagent equivalents based on the same scaling factor.
 * - Adjusts the readonly and editable states of relevant input fields:
 *   - The limiting reactant's equivalent is set to 1 and marked as readonly.
 *   - The limiting reactant's mass field becomes editable.
 */
function updateEquivalentsRelativeToNewLimitingReactant(x) {
  // x is the limiting reactant integer and i is the integer for each reactant
  // updates the classes of reactants - controls readonly and style properties when limiting reactant changes
  // Gets the change factor required to change the new limiting reactant equivalents to 1 and update other equivalents
  let $limitingReactantEquivalent = $("#js-reactant-equivalent" + x);
  let oldEquiv = getVal($limitingReactantEquivalent);
  let changeFactor = oldEquiv / 1;
  // Scales the equivalents of all reactants by this change factor
  for (let i = 1; i < reactionTable.numberOfReactants + 1; i++) {
    let $reactantEquivalent = $("#js-reactant-equivalent" + i);
    let newEquivalentValue = getVal($reactantEquivalent) / changeFactor;
    $reactantEquivalent.val(newEquivalentValue);
  }
  // Scales and updates equivalents for reagents
  let reagentNumber = Number(getVal($("#js-number-of-reagents")));
  for (let i = 1; i < reagentNumber + 1; i++) {
    let $reagentEquivalent = $("#js-reagent-equivalent" + i);
    let newEquivalentValue = getVal($reagentEquivalent) / changeFactor;
    $reagentEquivalent.val(newEquivalentValue);
  }
  // updates the limiting reactant equivalent and mass fields
  $limitingReactantEquivalent
    .val(1)
    .attr("readonly", true)
    .removeClass("editable-cell")
    .addClass("readonly-cell");
  $("#js-reactant-rounded-mass" + x)
    .attr("readonly", false)
    .removeClass("readonly-cell")
    .addClass("editable-cell");
}

/**
 * Adds an event listener to update the hidden variable storing the main product's table number upon main product
 * selection change.
 *
 * This function:
 * - Listens for clicks on elements with the class "js-main-product".
 * - Retrieves the selected main product's table number.
 * - Adjusts the table number by subtracting the total number of reactants, reagents, and solvents.
 * - Stores the updated value in the hidden input field "js-main-product-table-number".
 */
function updateMainProduct() {
  $(".js-main-product").click(function () {
    // updates the hidden variable main product table number
    let mainProductTableNumber = getNum(
      $("input[name='js-main-product']:checked"),
    );
    let numberOfReagents = getNum($("#js-number-of-reagents"));
    let numberOfSolvents = getNum($("#js-number-of-solvents"));
    $("#js-main-product-table-number").val(
      mainProductTableNumber -
        reactionTable.numberOfReactants -
        numberOfReagents -
        numberOfSolvents,
    );
  });
}

/**
 * Sets up event listeners for parameters affecting the limiting reactant calculations.
 *
 * This function defines an array of parameters (`changedParameters`) that impact molar calculations
 * when modified. For each parameter, it triggers:
 * - `autofillLimitingReactantAmount()`
 * - `autofillRoundedLimitingReactantAmount()`
 * - `autofillLimitingReactantMass()`
 *
 * These functions ensure that changes in mass or unit selection automatically update the relevant
 * limiting reactant values.
 */
function autoFillLimitingReactant() {
  let changedParameters = [
    ".js-reactant-rounded-masses",
    "#js-amount-unit",
    "#js-mass-unit",
  ];
  for (const param of changedParameters) {
    autofillLimitingReactantAmount(param);
    autofillRoundedLimitingReactantAmount(param);
    autofillLimitingReactantMass(param);
  }
}

/**
 * Adds an event listener to recalculate and autofill the limiting reactant amount.
 *
 * @param {string} changedParameter - A jQuery selector (`#id` or `.class`) that triggers recalculation
 *                                    when its value changes.
 *
 * This function:
 * - Listens for `input` or `change` events on the specified parameter.
 * - Retrieves the limiting reactant's table number.
 * - Gets the reactant's mass and molecular weight.
 * - Converts the mass using the appropriate unit factor (`massFactor`).
 * - Converts the molecular weight using the appropriate unit factor (`amountFactor`).
 * - Calculates the reactant amount and updates the corresponding input field.
 */
function autofillLimitingReactantAmount(changedParameter) {
  $(changedParameter).on("input change", function () {
    const limitingReactantTableNumber = getLimitingReactantTableNumber();
    let reactantMass = getVal(
      $("#js-reactant-rounded-mass" + limitingReactantTableNumber),
    );
    let reactantMolecularWeight = getVal(
      $("#js-reactant-molecular-weight" + limitingReactantTableNumber),
    );
    if (!reactantMolecularWeight) {
      let polymerMnSelector = $(
        "#js-reactant-mn" + limitingReactantTableNumber,
      );
      let polymerMn =
        getVal(polymerMnSelector) || polymerMnSelector.attr("placeholder");
      if (polymerMn) {
        // if polymer
        reactantMolecularWeight = polymerMn;
      }
    }
    const reactantMassUnit = getVal($("#js-mass-unit"));
    const reactantAmountUnit = getVal($("#js-amount-unit"));
    let reactantAmount =
      (reactantMass * massFactor[reactantMassUnit]) /
      (reactantMolecularWeight * amountFactor[reactantAmountUnit]);
    $("#js-reactant-amount" + limitingReactantTableNumber).val(reactantAmount);
  });
}

/**
 * Adds an event listener to round and autofill the limiting reactant amount when mass or units change.
 *
 * @param {string} changedParameter - A jQuery selector (`#id` or `.class`) that triggers recalculation
 *                                    when its value changes.
 *
 * This function:
 * - Listens for `input` or `change` events on the specified parameter.
 * - Retrieves the limiting reactant's table number.
 * - Gets the current reactant amount.
 * - Rounds the reactant amount using `roundedNumber()`.
 * - Updates the corresponding rounded amount input field.
 */
function autofillRoundedLimitingReactantAmount(changedParameter) {
  $(changedParameter).on("input change", function () {
    const limitingReactantTableNumber = getLimitingReactantTableNumber();
    let reactantAmount = getNum(
      $("#js-reactant-amount" + limitingReactantTableNumber),
    );
    let roundedReactantAmount = roundedNumber(reactantAmount);
    $("#js-reactant-rounded-amount" + limitingReactantTableNumber).val(
      roundedReactantAmount,
    );
  });
}

/**
 * Adds an event listener to recalculate and autofill the hidden input for the limiting reactant mass.
 *
 * @param {string} changedParameter - A jQuery selector (`#id` or `.class`) that triggers recalculation
 *                                    when its value changes.
 *
 * This function:
 * - Listens for `input` or `change` events on the specified parameter.
 * - Retrieves the limiting reactant's table number.
 * - Gets the updated rounded mass of the limiting reactant.
 * - Updates the corresponding hidden input field for the reactant mass.
 */
function autofillLimitingReactantMass(changedParameter) {
  $(changedParameter).on("input change", function () {
    let limitingReactantTableNumber = getLimitingReactantTableNumber();
    let limitingReactantMass = getVal(
      $("#js-reactant-rounded-mass" + limitingReactantTableNumber),
    );
    $("#js-reactant-mass" + limitingReactantTableNumber).val(
      limitingReactantMass,
    );
  });
}

/**
 * Event listener for raw amount changes. Calculates and autofills the amount for a reactant, reagent, solvent, or product.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - Retrieves the limiting reactant's amount.
 * - Gets the equivalent value for the specified component.
 * - If the component is a product, it converts the amount based on unit differences.
 * - Updates the corresponding input field for the component's amount.
 */
function autofillAmount(component, changedParameter, loopValue) {
  $(changedParameter).on("input change", function () {
    const limitingReactantTableNumber = getLimitingReactantTableNumber();
    const firstReactantAmount = getVal(
      $("#js-reactant-amount" + limitingReactantTableNumber),
    );
    let equivalentSelector = $("#js-" + component + "-equivalent" + loopValue);
    const equivalent =
      getVal(equivalentSelector) || equivalentSelector.attr("placeholder");
    let amount;
    // product can have different units to reactant, hence different equation
    if (component === "product") {
      let reactantAmountUnit = getVal($("#js-amount-unit"));
      let productAmountUnit = getVal($("#js-product-amount-unit"));
      amount =
        (firstReactantAmount * equivalent * amountFactor[reactantAmountUnit]) /
        amountFactor[productAmountUnit];
    } else {
      amount = firstReactantAmount * equivalent;
    }
    $("#js-" + component + "-amount" + loopValue).val(amount);
  });
}

/**
 * Event listener for rounded amount changes. Calculates and autofills the rounded moles for a reactant, reagent, solvent, or product.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - Retrieves the current amount for the specified component.
 * - Rounds the amount using `roundedNumber()`.
 * - Updates the corresponding input field for the rounded amount.
 */
function autofillRoundedAmount(component, changedParameter, loopValue) {
  // autofills the rounded amount for reactant, reagent, or product
  $(changedParameter).on("input change", function () {
    let amount = getNum($("#js-" + component + "-amount" + loopValue));
    let roundedAmount = roundedNumber(amount);
    $("#js-" + component + "-rounded-amount" + loopValue).val(roundedAmount);
  });
}

/**
 * Event listener for volume changes. Calculates and autofills the volume for a reactant, reagent, solvent, or product.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - Retrieves the mass, amount, density, and concentration for the specified component.
 * - Uses the `calcVolume()` function to calculate the volume based on the retrieved values.
 * - Updates the corresponding input field for the component's volume.
 */
function autofillVolume(component, changedParameter, loopValue) {
  // autofills the volume for reactant or reagent
  $(changedParameter).on("input change", function () {
    const mass = getVal($("#js-" + component + "-rounded-mass" + loopValue));
    let amount = getVal($("#js-" + component + "-amount" + loopValue));
    let density = getVal($("#js-" + component + "-density" + loopValue));
    let concentration = getVal(
      $("#js-" + component + "-concentration" + loopValue),
    );
    const volume = calcVolume(density, mass, concentration, amount);
    $("#js-" + component + "-volume" + loopValue).val(volume);
  });
}

/**
 * Event listener for volume changes. Calculates and autofills the rounded volume for a reactant, reagent, solvent, or product.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - Retrieves the current volume for the specified component.
 * - Rounds the volume using the `roundedNumber()` function.
 * - Updates the corresponding input field for the rounded volume.
 */
function autofillRoundedVolume(component, changedParameter, loopValue) {
  // autofills the rounded volume for reactant, reagent, or product
  $(changedParameter).on("input change", function () {
    let volume = getNum($("#js-" + component + "-volume" + loopValue));
    let roundedVolume = roundedNumber(volume);
    $("#js-" + component + "-rounded-volume" + loopValue).val(roundedVolume);
  });
}

/**
 * Event listener for mass changes. Calculates and autofills the mass for a reactant, reagent, or product.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - For products, calculates the mass based on the molecular weight, amount, and units.
 * - For reactants and reagents, calculates the mass based on the molecular weight, equivalent, and limiting reactant mass.
 * - Updates the corresponding input field for the component's mass.
 */
function autofillMass(component, changedParameter, loopValue) {
  $(changedParameter).on("input change", function () {
    let molecularWeight = getVal(
      $("#js-" + component + "-molecular-weight" + loopValue),
    );
    let mass;
    if (component === "product") {
      // product has different mass and amount units, hence need to include these in calculation
      let productAmountUnit = getVal($("#js-product-amount-unit"));
      let productMassUnit = getVal($("#js-product-mass-unit"));
      let productAmount = getVal($("#js-product-amount" + loopValue));
      let polymerMnSelector = $("#js-product-mn" + loopValue);
      let polymerMn =
        getVal(polymerMnSelector) || polymerMnSelector.attr("placeholder");
      if (polymerMn) {
        // if polymer
        molecularWeight = polymerMn;
      }
      mass =
        (molecularWeight * productAmount * amountFactor[productAmountUnit]) /
        massFactor[productMassUnit];
    } else {
      const limitingReactantTableNumber = getLimitingReactantTableNumber();
      let firstReactantMolecularWeight = getVal(
        $("#js-reactant-molecular-weight" + limitingReactantTableNumber),
      );
      let firstReactantMass = getVal(
        $("#js-reactant-rounded-mass" + limitingReactantTableNumber),
      );
      let equivalent = getVal(
        $("#js-" + component + "-equivalent" + loopValue),
      );
      let polymerMnSelector = $("#js-reactant-mn" + loopValue);
      let polymerMn =
        getVal(polymerMnSelector) || polymerMnSelector.attr("placeholder");
      if (polymerMn) {
        // if polymer
        molecularWeight = polymerMn;
      }
      mass =
        (molecularWeight * equivalent * firstReactantMass) /
        firstReactantMolecularWeight;
    }
    $("#js-" + component + "-mass" + loopValue).val(mass);
  });
}

/**
 * Event listener for mass changes. Calculates and autofills the rounded mass for a reactant, reagent, or product.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - Retrieves the current mass for the specified component.
 * - Rounds the mass using the `roundedNumber()` function.
 * - Updates the corresponding input field for the rounded mass.
 */
function autofillRoundedMass(component, changedParameter, loopValue) {
  // autofills the rounded mass for reactant, reagent, or product
  $(changedParameter).on("input change", function () {
    let mass = getNum($("#js-" + component + "-mass" + loopValue));
    let roundedMass = roundedNumber(mass);
    $("#js-" + component + "-rounded-mass" + loopValue).val(roundedMass);
  });
}

/**
 * Event listener for volume removal. Clears the rounded volume if both the density and concentration are empty.
 *
 * @param {string} component - The type of component being updated (`reactant`, `solvent`, `product`, or `reagent`).
 * @param {string} changedParameter - The jQuery selector that triggered the change (`#id` or `.class`).
 * @param {string} loopValue - The numeric identifier combined with the component name to form the HTML element ID (e.g., `reactant1`).
 *
 * This function:
 * - Listens for `input` or `change` events on the specified `changedParameter`.
 * - Checks if both density and concentration fields for the component are empty.
 * - If both fields are empty, clears the rounded volume field.
 */
function autoRemoveVolume(component, changedParameter, loopValue) {
  $(changedParameter).on("input change", function () {
    // if reactant density + mol conc are both empty then set volume to also be empty not just '0'
    let density = getVal($("#js-" + component + "-density" + loopValue));
    let concentration = getVal(
      $("#js-" + component + "-concentration" + loopValue),
    );
    if (concentration === "" && density === "") {
      $("#js-" + component + "-rounded-volume" + loopValue).val("");
    }
  });
}

/**
 * Sets up event listeners for reactant rows in the reaction table.
 * This function is called in `initialiseReactionTable()`.
 *
 * It autofills the fields for reactants, including volume, amount, mass, and equivalent,
 * based on changes in related parameters such as density, concentration, and units.
 * The function also handles specific behavior for the limiting reactant, ensuring its fields are properly updated.
 *
 * @returns {void}
 *
 * This function:
 * - Iterates over each reactant and sets up listeners for volume, mass, and amount fields.
 * - Applies calculations based on limiting reactant settings.
 * - Utilizes various autofill helper functions for each reactant field.
 */
function autofillReactantFields2() {
  // autofills each reactant field depending on number of reactants
  const component = "reactant";
  for (let i = 1; i < reactionTable.numberOfReactants + 1; i++) {
    let limitingReactantTableNumber = getLimitingReactantTableNumber();
    let limitingReactantMassID =
      "#js-reactant-rounded-mass" + limitingReactantTableNumber;
    let reactantDensityID = "#js-reactant-density" + i;
    let reactantConcentrationID = "#js-reactant-concentration" + i;
    let reactantEquivalentID = "#js-reactant-equivalent" + i;
    // iterate through to autofill for each reactant
    let volumeChangeParameters = [
      limitingReactantMassID,
      reactantDensityID,
      reactantConcentrationID,
      reactantEquivalentID,
      "#js-volume-unit",
      "#js-mass-unit",
    ];
    for (const param of volumeChangeParameters) {
      autofillVolume(component, param, i);
      autofillRoundedVolume(component, param, i);
    }
    autoRemoveVolume(component, reactantConcentrationID, i);
    autoRemoveVolume(component, reactantDensityID, i);
    // Stops auto-filling amount/mass for the limiting reactant - this uses a separate function defined earlier
    if (i !== Number(limitingReactantTableNumber)) {
      let amountMassChangeParameters = [
        reactantEquivalentID,
        limitingReactantMassID,
        "#js-amount-unit",
        "#js-mass-unit",
        "#js-reactant-mn" + i,
      ];
      for (const param of amountMassChangeParameters) {
        autofillAmount(component, param, i);
        autofillRoundedAmount(component, param, i);
        autofillMass(component, param, i);
        autofillRoundedMass(component, param, i);
      }
    }
  }
}

/**
 * Sets up an event listener for changes in a reagent field, which triggers the posting of reagent data.
 * The function is called for each reagent row and handles the dynamic updating of reagent data.
 *
 * @param {number} x - The index number of the reagent, used to generate the specific reagent ID (e.g., `js-reagent1`).
 *
 * This function:
 * - Listens for `keyup` and `change` events on the reagent input field.
 * - Retrieves the reagent name from the field.
 * - Posts the reagent data using the `postReagentData` function.
 */
function autofillReagentData(x) {
  let reagentID = "#js-reagent" + x;
  $(reagentID).on("keyup change", function () {
    let reagentName = getVal($(this));
    postReagentData(reagentName, x);
  });
}

/**
 * Posts reagent data to the server and handles the response to either update the reagent fields or handle errors.
 * This function is called when the reagent name is updated in the input field.
 * It handles cases such as reagent not found, CAS number not found, or successful reagent retrieval.
 *
 * @param {string} reagentName - The name of the reagent entered by the user.
 * @param {number} x - The index number of the reagent, used to generate the specific reagent element IDs (e.g., `js-reagent1`).
 *
 * @returns {Promise} Resolves with a message indicating the result of the operation:
 * - "CAS not found" if the reagent’s CAS number was not found.
 * - "Reagent not found" if the reagent is not found, and partial matches are returned.
 * - "Reagent found" if the reagent data is found and the fields are populated.
 *
 * This function:
 * - Sends an AJAX `POST` request to the server with the reagent data.
 * - If no reagent data is found, it returns early without making changes.
 * - If the CAS number is not found, it triggers the `novelCompoundInput` function for manual input.
 * - If the reagent is not found, it populates a `datalist` with partial matches.
 * - If the reagent is found, it populates various fields with the reagent's attributes and applies styling.
 * - Calls `autoSaveCheck()` to save the state after the reagent is found.
 */
function postReagentData(reagentName, x) {
  return new Promise(function (resolve) {
    let workbook = getVal($("#js-active-workbook"));
    let workgroup = getVal($("#js-active-workgroup"));
    $.ajax({
      url: "/_reagents",
      type: "post",
      data: {
        reagent: reagentName,
        number: x,
        workbook: workbook,
        workgroup: workgroup,
      },
      dataType: "json",
      success: function (response) {
        // if no reagent data to add dont do anything
        if (response === undefined) {
          return;
        }
        let $reagentID = $("#js-reagent" + x);
        let cas_not_found = response.cas_not_found;
        if (cas_not_found === true) {
          // now need to input new reagent
          novelCompoundInput("reagent", response.reagent, x);
          resolve("CAS not found");
          return;
        }
        let reagent_not_found = response.reagent_not_found;
        if (reagent_not_found === true) {
          // If the exact reagent is not found we populate the datalist with the partial matches.
          $("#js2-reagent" + x).empty();
          $.each(response.identifiers, function (i, item) {
            $("#js2-reagent" + x).append($("<option>").attr("value", item));
          });
          resolve("Reagent not found");
        } else {
          // if reagent found then apply styling and populate fields
          $("#js1-reagent" + x).empty();
          $("#js2-reagent" + x).empty(); // empties the list when a reagent is found
          $reagentID.val(response.name);
          $reagentID
            .attr("readonly", true)
            .removeClass("editable-cell add-highlight-unfilled-cell")
            .addClass("readonly-cell remove-highlight-filled-cell");
          let y = response.number;
          $("#js-report-reagent" + y).show();
          // fill new fields with data from response containing reagent attributes
          let reagentFieldList = [
            "density",
            "concentration",
            "molecular-weight",
            "hazards",
            "primary-key",
            "smiles",
          ];
          fillData(response, "reagent", reagentFieldList, y);
          autoChangeRequiredStylingValidCompound("reagent", y);
          autoSaveCheck();
          resolve("Reagent found");
        }
      },
    });
  });
}

/**
 * Sets up event listeners and autofills fields for a reagent row in the reaction table.
 * This function configures all the relevant input fields for a reagent based on the limiting reactant and various parameters.
 *
 * It handles autofilling the amount, mass, volume, and related values based on changes in mass, amount, unit, equivalents, density, and concentration.
 *
 * @param {number} i - The index of the reagent row to configure (e.g., reagent1, reagent2, etc.).
 *
 * This function:
 * - Sets up autofill functionality for amount, rounded amount, mass, and rounded mass for the reagent.
 * - Sets up autofill functionality for volume and rounded volume based on density and concentration.
 * - Removes volume values if both density and concentration are empty for the reagent.
 *
 * It performs the following steps:
 * 1. Initializes the relevant reagent field IDs and limiting reactant field IDs.
 * 2. Sets up autofill for amount, rounded amount, mass, and rounded mass based on changes in mass, amount, and equivalent values.
 * 3. Sets up autofill for volume and rounded volume based on changes in density, concentration, and other relevant parameters.
 * 4. Calls `autoRemoveVolume` to clear the volume field if both concentration and density are empty for the reagent.
 */
function autofillReagentFields2(i) {
  const component = "reagent";
  let limitingReactantTableNumber = getLimitingReactantTableNumber();
  let limitingReactantMassID =
    "#js-reactant-rounded-mass" + limitingReactantTableNumber;
  let reagentEquivalentID = "#js-reagent-equivalent" + i;
  let reagentDensityID = "#js-reagent-density" + i;
  let reagentConcentrationID = "#js-reagent-concentration" + i;
  const amountMassChangeParameters = [
    limitingReactantMassID,
    reagentEquivalentID,
    "#js-amount-unit",
    "#js-mass-unit",
  ];
  for (const param of amountMassChangeParameters) {
    autofillAmount(component, param, i);
    autofillRoundedAmount(component, param, i);
    autofillMass(component, param, i);
    autofillRoundedMass(component, param, i);
  }
  const volumeChangeParameters = [
    limitingReactantMassID,
    reagentDensityID,
    reagentConcentrationID,
    reagentEquivalentID,
    "#js-volume-unit",
    "#js-mass-unit",
  ];
  for (const param of volumeChangeParameters) {
    autofillVolume(component, param, i);
    autofillRoundedVolume(component, param, i);
  }
  autoRemoveVolume("reagent", "#js-reagent-concentration" + i, i);
  autoRemoveVolume("reagent", "#js-reagent-density" + i, i);
}

function addNewReagent() {
  // get current number of reagents and plus one
  let $reagentNumber = $("#js-number-of-reagents");
  let reagentNumber = getNum($reagentNumber);
  reagentNumber++;
  let $reactantNumber = $("#js-number-of-reactants");
  let numberOfReactants = getNum($reactantNumber);
  let reagentTableNumber = numberOfReactants + reagentNumber;
  // updates reagent number html hidden input
  $reagentNumber.val(reagentNumber);
  $reactantNumber.val(numberOfReactants);
  // makes variable markup for new row to be appended
  let markup = $("#js-reagent-table-new-row")
    .html()
    .replace(/-x-/g, reagentNumber)
    .slice(8, -8);
  $("tbody#js-reagent-table").append(markup);
  $("#js-reagent-table-number" + reagentNumber).val(reagentTableNumber);
  let reagentPhysicalFormID = "js-reagent-physical-form" + reagentNumber;
  $("#js-physical-form-dropdown")
    .clone()
    .prop("id", reagentPhysicalFormID)
    .appendTo("#js-reagent-physical-form-dropdown-cell" + reagentNumber);
  autofillReagentData(reagentNumber);
  autofillReagentFields2(reagentNumber);
  updateStyling();
  updateProductTableNumber();
  // update solvent table numbers
  updateSolventTableNumbers();
  postReagentData("", reagentNumber);
}

/**
 * Removes a reagent from the reaction table and updates the relevant table numbers and field IDs for remaining reagents.
 * This function is triggered when the user clicks the "remove reagent" button.
 * It also handles updating the table numbers for reagents, solvents, and products, and reassigns IDs and event listeners to remaining reagents.
 *
 * @param {number} removedReagentNumber - The index (1-based) of the reagent to be removed from the table.
 *
 * This function:
 * - Removes the HTML row of the reagent being deleted.
 * - Loops through the remaining reagents, adjusting the table number and field IDs.
 * - Replaces the HTML elements with updated IDs and reattaches event listeners.
 * - Updates the reagent, solvent, and product table numbers.
 * - Calls `updateSolventTableNumbers()` and `updateProductTableNumber()` to ensure accurate table numbering across the entire form.
 */
async function removeReagent(removedReagentNumber) {
  // called from remove reagent button
  removedReagentNumber = Number(removedReagentNumber);
  let reactantNumber = getNum("#js-number-of-reactants");
  let reagentNumber = getNum("#js-number-of-reagents");
  // remove the reagent table row from the html
  let reagentTableRowID = "#js-reagent-table-row" + removedReagentNumber;
  $(reagentTableRowID).remove();
  // list of reagent ids that need to be updated
  const reagentIDsToUpdate = [
    "js-reagent-table-number",
    "remove-reagent",
    "js-reagent-table-row",
    "js-reagent",
    "js2-reagent",
    "js-reagent-molecular-weight",
    "js-reagent-density",
    "js-reagent-concentration",
    "js-reagent-equivalent",
    "js-reagent-rounded-amount",
    "js-reagent-rounded-volume",
    "js-reagent-rounded-mass",
    "js-reagent-physical-form",
    "js-reagent-hazards",
    "js-report-reagent",
    "js-reagent-amount",
    "js-reagent-volume",
    "js-reagent-mass",
    "js-reagent-smiles",
    "js-reagent-primary-key",
  ];
  // skip if it is the reagent being removed
  if (reagentNumber !== removedReagentNumber) {
    // for each reagent after the reagent being removed, reduce id by 1
    for (let i = removedReagentNumber + 1; i < reagentNumber + 1; i++) {
      let j = String(i - 1);
      let iStr = String(i);
      $("#js-reagent-table-number" + iStr).val(i + reactantNumber - 1);
      $("#remove-reagent" + iStr).val(j);
      $("#js-reagent" + iStr).attr("list", "js2-reagent" + j);
      for (let reagentIDToUpdate of reagentIDsToUpdate) {
        // remove listener events by cloning elements
        let old_element = document.getElementById(reagentIDToUpdate + iStr);
        let new_element = old_element.cloneNode(true);
        old_element.parentNode.replaceChild(new_element, old_element);
        $("#" + reagentIDToUpdate + iStr).attr("id", reagentIDToUpdate + j);
      }
      // recreate listener events on the cloned elements with updated ids
      autofillReagentData(j);
      autofillReagentFields2(j);
      updateStyling();
    }
  }
  // update the reagent number before the product+solvent table numbers
  reagentNumber--;
  await updateTableNumber(reagentNumber, "#js-number-of-reagents");
  // update solvent table numbers
  updateSolventTableNumbers();
  // update product table numbers
  updateProductTableNumber();
}

/**
 * Updates the table numbers for the solvent rows based on the current number of reactants and reagents.
 * This function is called when the table numbers for reactants or reagents are updated, and it adjusts the
 * solvent table numbers accordingly.
 *
 * The solvent table numbers are calculated by adding the number of reactants and reagents to the current
 * solvent row index.
 *
 * This function will only update the solvent table numbers if the number of solvents is greater than zero.
 */
function updateSolventTableNumbers() {
  let reactantNumber = getNum($("#js-number-of-reactants"));
  let reagentNumber = getNum($("#js-number-of-reagents"));
  let solventNumber = getNum($("#js-number-of-solvents"));
  if (solventNumber !== 0) {
    for (let loopValue = 1; loopValue < solventNumber + 1; loopValue++) {
      let solventTableNumber = reactantNumber + reagentNumber + loopValue;
      $("#js-solvent-table-number" + loopValue).val(solventTableNumber);
    }
  }
}

/**
 * CARRY ON HERE
 * @param compoundNumber
 * @param idToUpdate
 * @returns {Promise<unknown>}
 */
function updateTableNumber(compoundNumber, idToUpdate) {
  return new Promise(function (resolve) {
    if (compoundNumber < 0) {
      compoundNumber = 0;
    }
    resolve($(idToUpdate).val(compoundNumber));
  });
}

function updateProductTableNumber() {
  let numberOfProducts = getNum($("#js-number-of-products"));
  let reagentNumber = getNum($("#js-number-of-reagents"));
  let solventNumber = getNum($("#js-number-of-solvents"));
  for (let loopValue = 1; loopValue < numberOfProducts + 1; loopValue++) {
    let productTableNumber =
      reactionTable.numberOfReactants +
      reagentNumber +
      solventNumber +
      loopValue;
    let mainProductTableNumber =
      reactionTable.numberOfReactants +
      reagentNumber +
      solventNumber +
      loopValue;
    $("#js-product-table-number" + loopValue).val(productTableNumber);

    $("#js-main-product" + loopValue).val(mainProductTableNumber);
  }
}

// function to set up lister events for solvent datalist
function datalist_initiate(solventInputID, solventDatalistID, solventNumber) {
  autofillSolventData(solventNumber);
  autofillSolventFields2();
  updateStyling();
  for (let option of document.getElementById(solventDatalistID).options) {
    option.onclick = function () {
      document.getElementById(solventInputID).value = option.value;
      document.getElementById(solventDatalistID).style.display = "none";
      document.getElementById(solventInputID).style.borderRadius = "5px";
      postSolventData(option.value, solventNumber, true);
    };
  }
  // hide the solvent dropdown if clicking on non-dropdown/input element
  let $solventInputID = "#" + solventInputID;
  $(document).on("click", function (event) {
    let $target = $(event.target);
    let $solventDatalistID = "#" + solventDatalistID;
    if (
      !$target.closest($solventDatalistID).length &&
      !$target.closest($solventInputID).length &&
      $($solventDatalistID).is(":visible")
    ) {
      $($solventDatalistID).hide();
    } else {
    }
  });

  // Filter the options based on input values
  function handler() {
    document.getElementById(solventDatalistID).style.display = "block";
    let userInput = document.getElementById(solventInputID).value.toUpperCase();
    for (let option of document.getElementById(solventDatalistID).options) {
      if (option.value.toUpperCase().indexOf(userInput) > -1) {
        option.style.display = "block";
      } else {
        option.style.display = "none";
      }
    }
  }

  document.getElementById(solventInputID).oninput = handler;
  document.getElementById(solventInputID).onclick = handler;
  // Tracking the options via keyboard up and down arrow keys and select an option using enter
  let currentFocus = -1;
  document.getElementById(solventInputID).onkeydown = function (e) {
    if (e.keyCode === 40) {
      do {
        currentFocus++;
        if (
          currentFocus >=
          document.getElementById(solventDatalistID).options.length
        )
          currentFocus = 0;
      } while (
        document.getElementById(solventDatalistID).options[currentFocus].style
          .display === "none"
      );
      addActive(document.getElementById(solventDatalistID).options);
    } else if (e.keyCode === 38) {
      do {
        currentFocus--;
        if (currentFocus < 0)
          currentFocus =
            document.getElementById(solventDatalistID).options.length - 1;
      } while (
        document.getElementById(solventDatalistID).options[currentFocus].style
          .display === "none"
      );
      addActive(document.getElementById(solventDatalistID).options);
    } else if (e.keyCode === 13) {
      e.preventDefault();
      if (currentFocus > -1) {
        if (document.getElementById(solventDatalistID).options)
          document
            .getElementById(solventDatalistID)
            .options[currentFocus].click();
      }
      currentFocus = -1;
    }
  };

  function addActive(x) {
    if (!x) return false;
    removeActive(x);
    if (currentFocus >= x.length) currentFocus = 0;
    if (currentFocus < 0) currentFocus = x.length - 1;
    x[currentFocus].classList.add("active");
  }

  function removeActive(x) {
    for (let i = 0; i < x.length; i++) {
      x[i].classList.remove("active");
    }
  }
}

function addNewSolvent() {
  // Get current number of solvents and plus one
  let solventNumberID = $("#js-number-of-solvents");
  let solventNumber = getVal(solventNumberID);
  solventNumber++;
  solventNumberID.val(solventNumber);
  // new row stored as the variable markup and replaced/sliced, and then appended to the solvent table
  let markup = $("#js-solvent-table-new-row")
    .html()
    .replace(/-x-/g, String(solventNumber))
    .slice(8, -8);
  $("tbody#js-solvent-table").append(markup);
  // clones the solvent datalist to the new row
  let solventInputID = "js-solvent" + solventNumber;
  let solventDatalistID = "js-solvent-datalist" + solventNumber;
  $("#js-solvent-datalist")
    .clone()
    .prop("id", solventDatalistID)
    .appendTo("#js-solvent-datalist-cell" + solventNumber);
  setColours();
  let solventPhysicalFormID = "js-solvent-physical-form" + solventNumber;
  // clones the physical form dropdown and appends
  $("#js-physical-form-dropdown")
    .clone()
    .prop("id", solventPhysicalFormID)
    .appendTo("#js-solvent-physical-form-dropdown-cell" + solventNumber);
  // initiate datalist for added solvent
  datalist_initiate(solventInputID, solventDatalistID, solventNumber);
  // update table number
  let solventTableNumber =
    getNum($("#js-number-of-reactants")) +
    getNum($("#js-number-of-reagents")) +
    Number(solventNumber);
  $("#js-solvent-table-number" + solventNumber).val(solventTableNumber);
  updateProductTableNumber();
}

async function removeSolvent(removedSolventNumber) {
  let solventNumber = getVal($("#js-number-of-solvents"));
  let reactantNumber = getNum($("#js-number-of-reactants"));
  let reagentNumber = getNum($("#js-number-of-reagents"));
  let solventTableRowID = "#js-solvent-table-row" + removedSolventNumber;
  $(solventTableRowID).remove();
  // list of reagent ids that need to be updated
  const solventIDsToUpdate = [
    "js-solvent-table-number",
    "remove-solvent",
    "js-solvent-table-row",
    "js-solvent",
    "js-solvent-molecular-weight",
    "js-solvent-density",
    "js-solvent-concentration",
    "js-solvent-physical-form-dropdown-cell",
    "js-solvent-physical-form",
    "js-solvent-hazards",
    "js-report-solvent",
    "js-solvent-volume",
    "js-solvent-datalist-cell",
    "js-solvent-datalist",
    "js-solvent-rounded-concentration",
    "go-to-solvent-guide",
    "js-solvent-primary-key",
  ];
  // skip if it is the reagent being removed
  if (solventNumber !== removedSolventNumber) {
    // for each solvent after the solvent being removed, reduce id by 1
    for (
      let i = Number(removedSolventNumber) + 1;
      i < Number(solventNumber) + 1;
      i++
    ) {
      let j = String(i - 1);
      let iStr = String(i);
      // remove listener events for all with solvent number greater than one removed
      let old_element = document.getElementById("js-solvent" + iStr);
      let new_element = old_element.cloneNode(true);
      old_element.parentNode.replaceChild(new_element, old_element);
      $("#js-solvent-table-number" + iStr).val(
        i + reactantNumber + reagentNumber - 1,
      );
      $("#remove-solvent" + iStr).val(j);
      $("#go-to-solvent-guide" + iStr).val(j);
      for (let solventIDToUpdate of solventIDsToUpdate) {
        $("#" + solventIDToUpdate + iStr).attr("id", solventIDToUpdate + j);
      }
    }
  }
  // reinitiate datalists
  if (solventNumber !== removedSolventNumber) {
    // for each solvent after the solvent being removed, reduce id by 1
    for (
      let i = Number(removedSolventNumber) + 1;
      i < Number(solventNumber) + 1;
      i++
    ) {
      let j = String(i - 1);
      datalist_initiate("js-solvent" + j, "js-solvent-datalist" + j, j);
    }
  }
  solventNumber--;
  await updateTableNumber(solventNumber, "#js-number-of-solvents");
  // update product table numbers
  updateProductTableNumber();
}

function initialiseReagentListeners() {
  let numberOfReagents = getNum($("#js-number-of-reagents"));
  for (let i = 1; i < numberOfReagents + 1; i++) {
    autofillReagentFields2(i);
  }
}

function initialiseSolventListeners() {
  let numberOfSolvents = getNum($("#js-number-of-solvents"));
  for (let i = 1; i < numberOfSolvents + 1; i++) {
    autofillSolventData(i);
    autofillSolventFields2(i);
    let solventInputID = "js-solvent" + i;
    let solventDatalistID = "js-solvent-datalist" + i;
    $("#js-solvent-datalist")
      .clone()
      .prop("id", solventDatalistID)
      .appendTo("#js-solvent-datalist-cell" + i);
    setColours();
    datalist_initiate(solventInputID, solventDatalistID, i);
  }
}

// Autofill solvent flag (colour code) and hazards depending on selected solvent
function autofillSolventData(x) {
  setColours();
  let solventID = "#js-solvent" + x;
  $(solventID).on("keyup change", function () {
    let solventName = getVal($(this));
    let oldValue = $(this).attr("oldValue");
    $(this).attr("oldValue", solventName);
    // post data if new value is different to previous value to prevent duplicate ajax calls
    if (oldValue !== solventName) {
      postSolventData(solventName, x);
    }
  });
}
function postSolventData(solventName, x, reload = false) {
  return new Promise(function (resolve) {
    let workbook = getVal($("#js-active-workbook"));
    let workgroup = getVal($("#js-active-workgroup"));
    $.ajax({
      url: "/_solvents",
      type: "post",
      data: {
        solvent: solventName,
        number: x,
        workbook: workbook,
        workgroup: workgroup,
      },
      dataType: "json",
      success: function (response) {
        if (response === undefined) {
          resolve("undefined");
          return;
        }
        if (!reload) {
          checkPCASolvents(solventName, x);
        }
        let y = response.num;
        let solvent = response.solvent;
        let newSolvent = response.new_solvent;
        let alertMessage = response.alert_message;
        if (alertMessage != "") {
          alert(alertMessage);
        }
        if (newSolvent) {
          novelCompoundInput("solvent", solvent, y);
          return;
        }
        let solventFieldList = ["hazards", "primary-key"];
        fillData(response, "solvent", solventFieldList, y);
        let solventID = "#js-solvent" + y;
        $(solventID).attr("value", solvent);
        $(solventID).val(solvent);
        $(solventID)
          .removeClass()
          .addClass(response.flag + " remove-highlight-filled-cell");
        setColours();
        autoChangeRequiredStylingValidCompound("solvent", y);
        autoSaveCheck();
        $("#js-solvent-datalist" + y).hide();
        resolve("solventFound");
      },
    });
  });
}

//Autofill solvent table number depending on number of reagents
function autofillSolventTableNumber(x, changedParameter) {
  $(changedParameter).click(function () {
    let numberOfReagents = getNum($("#js-number-of-reagents"));
    let solventTableNumberID = "#js-solvent-table-number" + x;
    let solventTableNumber =
      reactionTable.numberOfReactants + numberOfReagents + x;
    $(solventTableNumberID).val(solventTableNumber);
  });
}

// Solvent autofill functions
//Autofill solvent concentration in the hidden field
function autofillSolventConcentration(
  x,
  changedParameter,
  limitingReactantTableNumber,
) {
  $(changedParameter).on("input change", function () {
    let firstReactantAmount = getVal(
      $("#js-reactant-amount" + limitingReactantTableNumber),
    );
    let solventVolumeID = "#js-solvent-volume" + x;
    let solventVolume = getVal($(solventVolumeID));
    let solventConcentration = 0;
    const reactantAmountUnit = getVal($("#js-amount-unit"));
    const solventVolumeUnit = getVal($("#js-solvent-volume-unit"));
    if (solventVolume > 0) {
      // multiply by 1000 to go from concentration in mL to L / decimetres
      solventConcentration =
        (amountFactor[reactantAmountUnit] * firstReactantAmount * 1000) /
        (solventVolume * volumeFactor[solventVolumeUnit]);
    }
    let solventConcentrationID = "#js-solvent-concentration" + x;
    $(solventConcentrationID).val(solventConcentration);
  });
}

//Autofill rounded solvent concentration in the explicit field
function autofillRoundedSolventConcentration(x, changedParameter) {
  $(changedParameter).on("input change", function () {
    let solventConcentrationID = "#js-solvent-concentration" + x;
    let solventConcentration = getNum($(solventConcentrationID));
    let roundedSolventConcentration = roundedNumber(solventConcentration);
    let roundedSolventConcentrationID = "#js-solvent-rounded-concentration" + x;
    $(roundedSolventConcentrationID).val(roundedSolventConcentration);
  });
}

function autofillSolventFields2() {
  let numberOfSolvents = getVal($("#js-number-of-solvents"));
  if (numberOfSolvents > 0) {
    for (let i = 1; i < numberOfSolvents + 1; i++) {
      let limitingReactantTableNumber = getLimitingReactantTableNumber();
      let limitingReactantMass =
        "#js-reactant-rounded-mass" + limitingReactantTableNumber;
      let solventVolumeID = "#js-solvent-volume" + i;
      const solventConcentrationParameters = [
        limitingReactantMass,
        solventVolumeID,
        "#js-mass-unit",
        "#js-solvent-volume-unit",
        ".js-reactant-limiting",
      ];
      for (const param of solventConcentrationParameters) {
        autofillSolventConcentration(i, param, limitingReactantTableNumber);
        autofillRoundedSolventConcentration(i, param);
      }
      autofillSolventTableNumber(i, ".js-add-reagent");
      autofillSolventTableNumber(i, ".js-remove-reagent");
      autofillSolventTableNumber(i, ".js-add-solvent");
      autofillSolventTableNumber(i, ".js-remove-solvent");
    }
  }
}

function goToSolventGuide(sol_x) {
  let solventName = getVal($("#js-solvent" + sol_x));
  if (solventName) {
    window.open("/solvent_guide/" + solventName, "_blank").focus();
  } else {
    window.open("/solvent_guide", "_blank").focus();
  }
}

function goToSolventSurfer(sol_x) {
  let reactionClass = $("#reaction-class").val();
  let solventName = $("#js-solvent" + sol_x).val();

  localStorage.setItem("solventName", solventName);
  localStorage.setItem("reactionClass", reactionClass);

  window.open("/solvent_PCA");
}

async function checkPCASolvents(solventName, x) {
  let reactionClass = $("#reaction-class").val();
  let response = await fetch("/from_reaction_table/check_solvents", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
    body: JSON.stringify({
      class_selected: reactionClass,
      name_selected: solventName,
    }),
  });
  let PCASolvents = await response.json();
  let btn = document.getElementById("go-to-solvent-surfer");

  if (!PCASolvents.solvents.includes(solventName.toUpperCase())) {
    btn.disabled = true;
    btn.title = "Selected solvent is not included in the solvent surfer";
  } else {
    btn.disabled = false;
    btn.title = "Open solvent surfer in new tab";

    if (PCASolvents.substitutions) {
      if (
        confirm(
          "Suggested sustainable solvent switch:\n\n" +
            solventName +
            " \t " +
            PCASolvents.alternatives +
            "\n\nProceed to Solvent Surfer?",
        ) == true
      ) {
        goToSolventSurfer(x);
      }
    }
  }
}

// Product autofill
//Autofill product table number
function autofillProductTableNumber(changedParameter, loopValue) {
  $(changedParameter).click(function () {
    let numberOfReagents = getNum($("#js-number-of-reagents"));
    let numberOfSolvents = getNum($("#js-number-of-solvents"));
    let productTableNumber =
      reactionTable.numberOfReactants +
      numberOfReagents +
      numberOfSolvents +
      loopValue;
    let mainProductTableNumber =
      reactionTable.numberOfReactants +
      numberOfReagents +
      numberOfSolvents +
      loopValue;
    $("#js-product-table-number" + loopValue).val(productTableNumber);
    $("#js-main-product" + loopValue).val(mainProductTableNumber);
  });
}

function autofillProductFields2() {
  let productNumber = getNum($("#js-number-of-products"));
  let limitingReactantTableNumber = getLimitingReactantTableNumber();
  let limitingReactantMassID =
    "#js-reactant-rounded-mass" + limitingReactantTableNumber;
  let limitingReactantEquivalentID =
    "#js-reactant-equivalent" + limitingReactantTableNumber;
  const component = "product";
  for (let i = 1; i < productNumber + 1; i++) {
    let productAmountParameters = [
      limitingReactantEquivalentID,
      limitingReactantMassID,
      "#js-product-equivalent" + i,
      "#js-mass-unit",
      "#js-product-amount-unit",
      "#js-product-mass-unit",
      ".js-reactant-limiting",
      "#js-product-mn" + i,
    ];
    for (const param of productAmountParameters) {
      autofillAmount(component, param, i);
      autofillRoundedAmount(component, param, i);
      autofillMass(component, param, i);
      autofillRoundedMass(component, param, i);
    }
    autofillProductTableNumber(".js-add-reagent", i);
    autofillProductTableNumber(".js-remove-reagent", i);
    autofillProductTableNumber(".js-add-solvent", i);
    autofillProductTableNumber(".js-remove-solvent", i);
  }
}

function updateRequiredStylingLimited() {
  let limitingReactantTableNumber = getLimitingReactantTableNumber();
  let RoundedReactantMassID =
    "#js-reactant-rounded-mass" + limitingReactantTableNumber;
  autoChangeRequiredStyling2(RoundedReactantMassID);
  for (let i = 1; i <= reactionTable.numberOfReactants; i++) {
    let reactantEquivalentField = $("#js-reactant-equivalent" + i);
    autoChangeRequiredStyling2(reactantEquivalentField);
  }
}

/**
 * Highlights red/removes red highlights when essential cell is not filled in to draw user attention.
 * @param {string} styleParameterID - # then element ID, e.g., "#element-1" to be used in JQuery Selector
 * @param {Array} [excludedNullValues=[]] - Typically null values which should not be treated as null for this parameter
 */
function autoChangeRequiredStyling2(styleParameterID, excludedNullValues = []) {
  // doesnt require change parameter
  const defaultNulLValues = ["-select-", "0", "", "-", 0];
  // remove excluded nullValues from the default nullValues
  let nullValues = defaultNulLValues.filter(
    (item) => !excludedNullValues.includes(item),
  );
  let parameterValue = getVal($(styleParameterID));
  if (nullValues.includes(parameterValue)) {
    $(styleParameterID)
      .removeClass("remove-highlight-filled-cell")
      .addClass("add-highlight-unfilled-cell");
  } else {
    $(styleParameterID)
      .removeClass("add-highlight-unfilled-cell")
      .addClass("remove-highlight-filled-cell");
  }
}

/**
 * Calls the function to update red highlight when there is an input change to the cell changed parameter
 * @param {string} changedParameter - # then element ID, e.g., "#element-1" to be used in JQuery Selector
 * @param {Array} [excludedNullValues=[]] Typically null values which should not be treated as null for this parameter
 */
function autoChangeRequiredStyling(changedParameter, excludedNullValues = []) {
  // update styling of current reaction table (for reaction table reload)
  autoChangeRequiredStyling2(changedParameter, excludedNullValues);
  // then set up the listener for future changes
  $(changedParameter).on("input change", function () {
    autoChangeRequiredStyling2(changedParameter, excludedNullValues);
  });
}

/**
 * Function called to change styling on reagent and solvet changing
 * @param changedParameter
 * @param changedStyling
 */
function styleValidReagent(changedParameter, changedStyling) {
  if (getVal($(changedParameter)) === "") {
    $(changedStyling)
      .removeClass("remove-highlight-filled-cell")
      .addClass("add-highlight-unfilled-cell");
  } else {
    $(changedStyling)
      .removeClass("add-highlight-unfilled-cell")
      .addClass("readonly-cell remove-highlight-filled-cell");
  }
}

function autoChangeRequiredStylingValidCompound(component, loop_value) {
  // Catches partially filled in reagent/solvent name box - will highlight red unless a valid compound is entered
  let changedParameter = "#js-" + component + "-hazards";
  changedParameter = changedParameter.concat(String(loop_value));
  let changedStyling = "#js-" + component;
  changedStyling = changedStyling.concat(String(loop_value));
  styleValidReagent(changedParameter, changedStyling);
  $(changedStyling).on("input change", function () {
    styleValidReagent(changedParameter, changedStyling);
  });

  if (component === "solvent") {
    $(changedStyling).removeClass("readonly-cell");
  }
}

function updateStyling() {
  for (let i = 1; i < reactionTable.numberOfReactants + 1; i++) {
    autoChangeRequiredStyling("#js-reactant-physical-form" + i);
    autoChangeRequiredStyling("#js-reactant-rounded-mass" + i);
    autoChangeRequiredStyling("#js-reactant-equivalent" + i);
  }
  let numberOfReagents = getNum($("#js-number-of-reagents"));
  for (let i = 1; i < numberOfReagents + 1; i++) {
    autoChangeRequiredStyling("#js-reagent-physical-form" + i);
    autoChangeRequiredStyling("#js-reagent-equivalent" + i);
    autoChangeRequiredStyling("js-reagent-hazards" + i);
    autoChangeRequiredStylingValidCompound("reagent", i);
  }
  let numberOfSolvents = getNum($("#js-number-of-solvents"));
  for (let i = 1; i < numberOfSolvents + 1; i++) {
    autoChangeRequiredStyling("#js-solvent-physical-form" + i);
    autoChangeRequiredStyling("#js-solvent-volume" + i);
    autoChangeRequiredStyling("#js-solvent" + i);
    autoChangeRequiredStylingValidCompound("solvent", i);
  }
  for (let i = 1; i < reactionTable.numberOfProducts + 1; i++) {
    autoChangeRequiredStyling("#js-product-physical-form" + i);
  }
}

// utility functions
function commonGlobalVariables() {
  // returns constant global variables as part of the reactionTable object accessed by reactionTable.numberOfReactants
  return [
    {
      // ReactionTable object
      numberOfReactants: getNum($("#js-number-of-reactants")),
      numberOfProducts: getNum($("#js-number-of-products")),
    },
    // Unit factors objects- amount, mass, volume
    { mol: 1, mmol: 10 ** -3, μmol: 10 ** -6 },
    { g: 1, mg: 10 ** -3, μg: 10 ** -6 },
    { mL: 1, μL: 10 ** -3 },
  ];
}

/**
 * Rounds a number to a useful number of decimal places
 *
 * @param x {number|string} the number being rounded can be type number or type string. returns original string if NaN
 * @returns {number} the rounded number
 */
function roundedNumber(x) {
  // rounds numbers depending on magnitude and returns the original variable if not a number
  x = Number(x);
  if (isNaN(x)) {
    return x;
  }
  let roundedX;
  if (x < 10) {
    roundedX = x.toFixed(2);
  } else if (x >= 10 && x < 100) {
    roundedX = x.toFixed(1);
  } else {
    roundedX = x.toFixed();
  }
  return roundedX;
}

/**
 * We take the mass of compound we need and calculate and return volume of compound that is needed to get this mass
 * Either density or concentration of the compound is used with preference to density if both are used.
 *
 * @param density {number}
 * @param mass {Number}
 * @param concentration {Number}
 * @param amount {number} in moles
 * @returns {number}
 */
function calcVolume(density, mass, concentration, amount) {
  // calculates volume for reactant or reagent
  const reactantMassUnit = getVal($("#js-mass-unit"));
  const reactantVolumeUnit = getVal($("#js-volume-unit"));
  const reactantAmountUnit = getVal($("#js-amount-unit"));
  let volume;
  if (density > 0 && concentration == 0) {
    // divide by volumefactor
    volume =
      (mass * massFactor[reactantMassUnit]) /
      (density * volumeFactor[reactantVolumeUnit]);
  } else if (density == 0 && concentration > 0) {
    // divide by 1000 because mL is 1 not 0.001
    volume =
      (amountFactor[reactantAmountUnit] * amount) /
      ((concentration * volumeFactor[reactantVolumeUnit]) / 1000);
  } else {
    volume = 0;
  }
  return volume;
}

/**
 * For a particular compound this function fills in data. The compound is identified by component+y e.g., js-solvent1
 * We create the HTML element ID for each element by adding a field to the compound id. e.g., js-solvent-hazards1
 * We then use the field as the key (with some exceptions covered by IF statements) to get the value from the response.
 *
 * @param response {JSON} the response json from reagents or solvents/routes.py
 * @param component {string} the component type. e.g., reagent or product, etc.
 * @param fieldList {Array<string>} the properties such as hazards, molecular weight to be filled in
 * @param y {string} the number of the current component. as a string to be concatenated with component+field
 */
function fillData(response, component, fieldList, y) {
  // fill fields with response data for solvents or reagents added to the reaction table
  for (let field of fieldList) {
    // eg component='reagent', field=density, x=1 => "#js-reagent-density1"
    let fieldID = "#js-" + component + "-" + field + y;
    if (field === "molecular-weight") {
      $(fieldID).val(response["molWeight"]);
    } else if (field === "primary-key") {
      $(fieldID).val(response["primary_key"]);
    } else {
      $(fieldID).val(response[field]);
    }
  }
}
