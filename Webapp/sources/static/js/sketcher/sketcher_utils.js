/**
 * Checks which sketcher is selected by the radio buttons and calls the function to switch to the selected sketcher.
 */
async function switchActiveEditor() {
  // in future export structures between the two when user changes
  let selectedSketcher = $('input[name="sketcher-select"]:checked').attr("id");
  if (selectedSketcher === "marvin-select") {
    await switchToMarvinSketcher();
  } else if (selectedSketcher === "ketcher-select") {
    await switchToKetcherSketcher();
  }
}

/**
 * Checks which sketcher is selected by the radio buttons and calls the function to export SMILES from the active editor.
 * @returns {Promise<string>} - The SMILES string from the structures drawn in the active editor
 */
async function exportSmilesFromActiveEditor() {
  let selectedSketcher = $('input[name="sketcher-select"]:checked').attr("id");
  let smiles;
  if (selectedSketcher === "marvin-select") {
    smiles = await exportSmilesFromMarvin();
  } else if (selectedSketcher === "ketcher-select") {
    smiles = await exportSmilesFromKetcher();
  }
  return smiles;
}

/**
 * Checks which sketcher is selected by the radio buttons and calls the function to export image from the active editor.
 * @returns {Promise<void>}
 */
async function exportImageFromActiveEditor() {
  let selectedSketcher = $('input[name="sketcher-select"]:checked').attr("id");
  if (selectedSketcher === "marvin-select") {
    exportMarvinImage();
  } else if (selectedSketcher === "ketcher-select") {
    let smiles = await exportSmilesFromKetcher();
    exportKetcherImage(smiles);
  }
}

/**
 * Causes a pause in the code to give time for sketcher code to laod.
 */
function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

/**
 * Replaces the smiles symbols with a word equivalent prior to conserve meaning during url transfer
 * @param {string} reaction - The reaction SMILES string
 * @returns {string} The reaction string with the SMILES symbols replaced with word equivalents
 */
function replaceSmilesSymbols(reaction) {
  // replace the plus, minus and sharp signs from a smiles string to trasfer it properly
  reaction = reaction.replace(/\+/g, "plus");
  reaction = reaction.replace(/-/g, "minus");
  return reaction.replace(/#/g, "sharp");
}
/**
 * Separates the reaction SMILES string into separate components, removing reagents.
 *
 * @param {string} smiles The reaction SMILES string
 * @returns {Array} An array containing reactants and products, each as an array
 */
function reactionSmilesToReactantsAndProductsSmiles(smiles) {
  let reaction = smiles.split(" |")[0];
  reaction = replaceSmilesSymbols(reaction);
  //split into reactants, solvents and product
  let array = reaction.split(">");
  let reactants = array[0].split(".");
  let products = array[2].split(".");
  // if "f:" is at the end of the smiles string, it contains an ionic compound and needs further processing
  if (smiles.includes("f:") === true) {
    [reactants, products] = processSalts(smiles, reactants, products);
  }
  return [reactants, products];
}

/**
 * Processes cxsmiles (Chemaxon SMILES from marvin JS) that contain an ionic compound.
 * @param {string} smiles - Reaction SMILES string
 * @param {Array} reactants - Reactants SMILES array
 * @param {Array} products - Products SMILES array
 * @returns {Array} Array containing reactants and products, each as an array
 */
function processSalts(smiles, reactants, products) {
  // split after the comma to see how many ionic compounds in reaction smiles
  let ionList = smiles.split("f:")[1].split("|")[0].split(",");
  // for each ionic compound, join all ions to make the full salt and replace ions in reactant/product list
  let productSaltDictList = [];
  let reactantSaltDictList = [];
  for (let adjacentIons of ionList) {
    // make list of index, and smiles, then join to make a string of the full salt from the ions
    let ionIndexList = adjacentIons.split(".").map(Number);
    let ionCompoundsSmiles = []; // [Pd+, OAc-, OAc-]
    for (let idx of ionIndexList) {
      ionCompoundsSmiles.push(compounds[idx]);
    }
    let saltSmiles = ionCompoundsSmiles.join(".");
    // determine if reactant or product
    let rxnComponent;
    if (rl > ionIndexList[0]) {
      rxnComponent = "reactant";
    } else {
      rxnComponent = "product";
    }
    // ions in previous salts for only same rxn component (reactant or product) - to determine insert position
    // when splicing into reactant/product list.
    let ionsInPreviousSalts = 0;
    let insertPosition;
    if (rxnComponent === "reactant") {
      for (let previousSaltDict of reactantSaltDictList) {
        ionsInPreviousSalts += previousSaltDict["numberOfIons"];
      }
      insertPosition =
        ionIndexList[0] - ionsInPreviousSalts + reactantSaltDictList.length;
      reactantSaltDictList.push({
        saltSmiles: saltSmiles,
        numberOfIons: ionCompoundsSmiles.length,
        insertPosition: insertPosition,
      });
    }
    if (rxnComponent === "product") {
      for (let previousSaltDict of productSaltDictList) {
        ionsInPreviousSalts += previousSaltDict["numberOfIons"];
      }
      insertPosition =
        ionIndexList[0] - ionsInPreviousSalts + productSaltDictList.length - rl;
      productSaltDictList.push({
        saltSmiles: saltSmiles,
        numberOfIons: ionCompoundsSmiles.length,
        insertPosition: insertPosition,
      });
    }
  }
  // splice the salts into the lists replacing the ions that they contain
  for (let saltDict of reactantSaltDictList) {
    reactants.splice(
      saltDict["insertPosition"],
      saltDict["numberOfIons"],
      saltDict["saltSmiles"],
    );
  }
  for (let saltDict of productSaltDictList) {
    products.splice(
      saltDict["insertPosition"],
      saltDict["numberOfIons"],
      saltDict["saltSmiles"],
    );
  }
  return [reactants, products];
}
/**
 * Removes reagents between the two '>' symbols from the reaction SMILES
 * @param {string} smiles - Reaction SMILES string
 * @returns {string}  Reaction SMILES string without reagents
 */
function removeReagentsFromSmiles(smiles) {
  // remove reagents from reaction smiles
  let smiles2 = smiles.split(">").slice(-1);
  let smiles1 = smiles.split(">")[0];
  return smiles1 + ">>" + smiles2;
}
/**
 * Shows the loading circle for the chemical sketcher
 */
function showSketcherLoadingCircle() {
  $(".loading-circle").css("display", "block");
}
/**
 * Hides the loading circle for the chemical sketcher
 */
function hideSketcherLoadingCircle() {
  $(".loading-circle").css("display", "none");
}
/**
 * If reaction page is not in demo or tutorial mode, saving is enabled
 * @returns {boolean} returns true if saving is enabled, false if saving is disabled.
 */
function checkIfSaveEnabled() {
  let demo_mode = $("#js-demo").val();
  let tutorial_mode = $("#js-tutorial").val();
  return demo_mode !== "demo" && tutorial_mode !== "yes";
}
