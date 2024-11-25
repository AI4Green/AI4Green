/**
 * Makes a ketcher sketcher, hides it and sets up the autosave
 * @param {string} [width="1020px"] - width of the sketcher in pixels
 */
function setupNewKetcherSketcher(width = "1020px") {
  let ketcherHTML = `<iframe id="ketcher-editor" src=/static/ketcher/index.html width=${width} height="480px"></iframe>`;
  $("#ketcher-sketcher").empty().append(ketcherHTML).hide();
  $("#ketcher-select").prop("disabled", false);
}

/**
 * Displays the ketcher sketcher nd exports SMILES from Marvin to Ketcher
 * @return {Promise<void>}
 */
async function switchToKetcherSketcher() {
  await exportSmilesFromMarvinToKetcher();
  displayKetcher();
}

/**
 * Shows the ketcher sketcher, hides marvin and the loading circle
 */
function displayKetcher() {
  $("#marvin-sketcher").hide();
  $("#ketcher-sketcher").show();
  hideSketcherLoadingCircle();
}

/**
 * Autosave when the sketcher is edited
 */
function setupKetcherAutosave() {
  let ketcher = getKetcher();
  ketcher.editor.subscribe("change", function () {
    setTimeout(autoSaveCheck(null, true), 100);
  });
}

/**
 * Exports structures from Marvin to Ketcher vin via SMILES
 */
async function exportSmilesFromMarvinToKetcher() {
  if (marvin) {
    let smiles = await exportSmilesFromMarvin();
    if (smiles) {
      let ketcher = getKetcher();

      ketcher.setMolecule(smiles);
    }
  }
}

/**
 * Exports the current contetns of Ketcher as SMILES
 * @return {Promise<string>} Reaction SMILES string where reagents are optional and later discarded
 * in format: "reactant1.reactantX>reagent1.reagentX>product1.productX"
 */
async function exportSmilesFromKetcher() {
  let ketcher = getKetcher();
  return ketcher.getSmiles();
}

/**
 * Enters the example reaction into the Ketcher editor
 */
function ketcherExampleSmiles() {
  let smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1";
  let ketcher = getKetcher();
  ketcher.setMolecule(smiles);
}

/**
 * Gets the Ketcher editor as an object with access to its methods
 * @return {Object} Ketcher Object
 */
function getKetcher() {
  let ketcherFrame = document.getElementById("ketcher-editor");
  if ("contentDocument" in ketcherFrame)
    return ketcherFrame.contentWindow.ketcher;
}

/**
 * Generates an image from the SMILES string and saves it to a hidden HTML input as a Blob.
 * The image is then shrunk to match the width of the images exported from MarvinJS.
 * @param {string} smiles - The SMILES string representing the molecule.
 * @returns {Promise<string>} A Promise that resolves to a Base64-encoded string representing the generated image.
 */
async function exportKetcherImage(smiles) {
  try {
    const blob = await generateImageBlob(smiles);
    const shrunkenBlob = await shrinkBlobImage(blob, 600, 400, 100);
    return await convertBlobToBase64(shrunkenBlob);
  } catch (error) {
    console.error("Error generating image:", error);
    return "error generating image";
  }
}

/**
 * Generates an image Blob from the provided SMILES string.
 * @param {string} smiles - The SMILES string representing the molecule.
 * @returns {Promise<Blob>} A Promise that resolves to the generated image Blob.
 */
function generateImageBlob(smiles) {
  let ketcher = getKetcher();
  return new Promise((resolve, reject) => {
    ketcher.generateImage(smiles).then(resolve).catch(reject);
  });
}

/**
 * Shrinks the provided Blob image while maintaining the aspect ratio.
 * @param {Blob} blob - The Blob image to be shrunk.
 * @param {number} maxWidth - The maximum width of the shrunk image.
 * @param {number} maxHeight - The maximum height of the shrunk image.
 * @param {number} padding - The amount of whitespace to add above and below the image.
 * @returns {Promise<Blob>} A Promise that resolves to the shrunk Blob image.
 */
function shrinkBlobImage(blob, maxWidth, maxHeight, padding) {
  return new Promise((resolve) => {
    const img = new Image();
    img.onload = function () {
      const canvas = document.createElement("canvas");
      const ctx = canvas.getContext("2d");
      let width = img.width;
      let height = img.height;

      // Calculate new dimensions to maintain the aspect ratio
      if (width > maxWidth) {
        height *= maxWidth / width;
        width = maxWidth;
      }
      if (height > maxHeight) {
        width *= maxHeight / height;
        height = maxHeight;
      }

      // Add padding to the image height
      height += 2 * padding;

      // Set canvas dimensions to match the new image size
      canvas.width = width;
      canvas.height = height;

      // Draw the image with padding on the canvas
      ctx.fillStyle = "white";
      ctx.fillRect(0, 0, width, height);
      ctx.drawImage(img, 0, padding, width, height - 2 * padding);

      // Convert the canvas to a Blob with the new image size
      canvas.toBlob(resolve, "image/jpeg", 1.0);
    };
    img.src = URL.createObjectURL(blob);
  });
}

/**
 * Converts the provided Blob image to a Base64-encoded string.
 * @param {Blob} blob - The Blob image to be converted.
 * @returns {Promise<string>} A Promise that resolves to the Base64-encoded string.
 */
function convertBlobToBase64(blob) {
  const reader = new FileReader();
  return new Promise((resolve) => {
    reader.onloadend = () => resolve(reader.result);
    reader.readAsDataURL(blob);
  });
}
