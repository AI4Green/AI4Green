/**
 * Makes a ketcher sketcher, hides it and sets up the autosave
 * @param {string} [width="1020px"] - width of the sketcher in pixels
 */
function setupNewKetcherSketcher(width="1020px"){
    let ketcherHTML = `<iframe id="ketcher-editor" src=/static/ketcher/index.html width=${width} height="480px"></iframe>`
    $("#ketcher-sketcher").empty().append(ketcherHTML).hide();
}

/**
 * Displays the ketcher sketcher nd exports SMILES from Marvin to Ketcher
 * @return {Promise<void>}
 */
async function switchToKetcherSketcher() {
    await exportSmilesFromMarvinToKetcher()
    displayKetcher()
}

/**
 * Shows the ketcher sketcher, hides marvin and the loading circle
 */
function displayKetcher(){
    $("#marvin-sketcher").hide()
    $("#ketcher-sketcher").show()
    hideSketcherLoadingCircle()
}

/**
 * Autosave when the sketcher is edited
 */
function setupKetcherAutosave(){
    let ketcher = getKetcher()
    ketcher.editor.subscribe('change', function(){
        setTimeout(autoSaveCheck(null, true), 100)
    });
}

/**
 * Exports structures from Marvin to Ketcher vin via SMILES
 */
async function exportSmilesFromMarvinToKetcher(){
    if (marvin) {
        let smiles = await exportSmilesFromMarvin()
        if (smiles){
            let ketcher = getKetcher()

            ketcher.setMolecule(smiles);
        }
    }
}

/**
 * Exports the current contetns of Ketcher as SMILES
 * @return {Promise<string>} Reaction SMILES string where reagents are optional and later discarded
 * in format: "reactant1.reactantX>reagent1.reagentX>product1.productX"
 */
async function exportSmilesFromKetcher(){
    let ketcher = getKetcher()
    return ketcher.getSmiles()
}

/**
 * Enters the example reaction into the Ketcher editor
 */
function ketcherExampleSmiles(){
    let smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1";
    let ketcher = getKetcher()
    ketcher.setMolecule(smiles);
}

/**
 * Gets the Ketcher editor as an object with access to its methods
 * @return {Object} Ketcher Object
 */
function getKetcher(){
    let ketcherFrame = document.getElementById('ketcher-editor');
    if ('contentDocument' in ketcherFrame)
        return ketcherFrame.contentWindow.ketcher;
}

/**
 * Makes an image from the SMILES string and saves to a hidden HTML input as a Blob that has been shrunk to match the
 * width of the images exported from MarvinJS
 * @param {string }smiles
 */
function exportKetcherImage(smiles){
    let ketcher = getKetcher()
    ketcher.generateImage(smiles).then(function (source) {
        // Shrink the blob image and display the shrunken image
        shrinkBlobImage(source, 600, 400, function (shrunkenBlob) {
            sessionStorage.setItem("reactionSchemeImage", URL.createObjectURL(shrunkenBlob))
        });
    });
}

/**
 * // Function to shrink the blob image whilst keeping the aspect ratio
 * @param {Blob} blob - The blob image to be shrunk
 * @param {number} maxWidth - The maximum width of the shrunk image
 * @param {number} maxHeight - The maximum width of the shrunk image
 * @param {function(blob): void} callback - The anonymous function where the shrunken blob is used
 */
function shrinkBlobImage(blob, maxWidth, maxHeight, callback) {
    const img = new Image();
    img.onload = function() {
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
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
    // Set canvas dimensions to match the new image size
    canvas.width = width;
    canvas.height = height;
    // white background
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, width, height);
    // Draw the image on the canvas with the new size
    ctx.drawImage(img, 0, 0, width, height);
    // Convert the canvas to a blob with the new image size
    canvas.toBlob(callback, 'image/jpeg', 1.0);
    };
    img.src = URL.createObjectURL(blob)
}
