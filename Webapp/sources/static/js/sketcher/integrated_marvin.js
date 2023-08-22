let marvin
/**
 * Assigns the global variable 'marvin' to the marvin sketcher object and prepares the sketcher ready for use
 */
async function setupNewMarvinSketcher(){
    let marvinKey = $("#marvin-js-key").val()
    marvin = await newMarvinSketcher("#marvin-sketcher", marvinKey)
    if (marvin) {
        marvin.setDisplaySettings({"toolbars": "education"});
        $("#marvin-sketcher").hide()
    }
}

/**
 * Loads the required scripts and returns the sketcher or aborts if scripts fail to load
 * @param {string} IDSelector - The ID of the new DOM element in the HTML
 * @param {string} marvinKey - The license key tied to the URL to access Marvin JS services
 * @return {Object} MarvinJS Chemical Editor object with access to the MarvinJS methods
 */
async function newMarvinSketcher(IDSelector, marvinKey) {
    let marvinKey1 = await getMarvinKey()
    try {
        await $.getScript("https://marvinjs.chemicalize.com/v1/client.js")
        await $.getScript(`https://marvinjs.chemicalize.com/v1/${marvinKey1}/client-settings.js`)
    } catch(error) {
        console.log(error.status)
        abortMarvinSketcherCreation()
        return null
    }
    return ChemicalizeMarvinJs.createEditor(IDSelector)
}

/**
 * Gets the marvin key from the application backend
 * @return {Promise<string>} - The marvin JS license key for accessing their API services
 */
async function getMarvinKey() {
    let marvinJsKey
    await $.ajax({
        method: "POST",
        url: "/get_marvinjs_key",
    }).then(function (response) {
        marvinJsKey = response.marvinjs_key
    });
    return marvinJsKey
}

/**
 * Shows the Marvin Sketcher and hides Ketcher and exports SMILES from ketcher to Marvin
 * @return {Promise<void>}
 */
async function switchToMarvinSketcher() {
    $("#ketcher-sketcher").hide()
    if (marvin) {
        await exportSmilesFromKetcherToMarvin()
        $("#marvin-sketcher").show()
    } else {
        alert("Marvin JS is temporarily unavailable")
    }
}

/**
 * Switches back to Ketcher, disables further switching and hides the loading circle again
 */
function abortMarvinSketcherCreation(){
    // runs when
    alert("Marvin JS is temporarily unavailable")
    $("#ketcher-select").prop("checked", true)
    $("#marvin-select").prop("disabled", true)
}

/**
 * Autosave when the sketcher is edited
 */
function setupMarvinAutosave(){
    marvin.on('molchange', function () {
        // sketcher arg to true to save just smiles
        setTimeout(autoSaveCheck(null, true), 100)
    });
}

/**
 * Exports structures from ketcher to marvin via SMILES
 * @return {Promise<void>}
 */
async function exportSmilesFromKetcherToMarvin(){
    let smiles
    if (getKetcher()) {
        smiles = await exportSmilesFromKetcher()
        if (smiles){
            marvin.importStructure("cxsmiles", smiles);
        }
    }
}

/**
 * Exports the current contents of Marvin JS as SMILES
 * @return {Promise<string>} Reaction SMILES string where reagents are optional and later discarded
 * in format: "reactant1.reactantX>reagent1.reagentX>product1.productX"
 */
function exportSmilesFromMarvin() {
    return marvin.exportStructure('cxsmiles', {'extra': 'f'}).then(function (smiles) {
        return smiles
    });
}

/**
 * Enters the example reaction into the Marvin editor
 */
function marvinExampleSmiles(){
    let smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1";
    marvin.importStructure("cxsmiles", smiles);
}

/**
 * Makes an image of the current reaction scheme and saves to a hidden HTML input
 */
function exportMarvinImage() {
    let settings = {'width': 600, 'height': 400};
    marvin.exportStructure("jpeg", settings).then(function (source) {
        sessionStorage.setItem("reactionSchemeImage", source)
        //$("#js-reaction-scheme-image").val(source);
    });
}