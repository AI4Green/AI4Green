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
