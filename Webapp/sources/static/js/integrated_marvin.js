// setup MarvinJS sketcher window
$(document).ready(function () {
    initialiseReactionSketcher()
});

async function initialiseReactionSketcher() {
    // make new marvin instance
    let marvin = await newMarvinSketcher("#marvin-test")
    marvin.setDisplaySettings({"toolbars": "education"});
    // if not creator lock editor
    sketcherLockHandler()
    // reload reaction smiles if present
    if (checkIfSaveEnabled()) {
        reloadHandler(marvin)
        autosaveInitialise(marvin)
    }
    //Sends a reaction from the molecular editor to the table
    sketcherEventListeners(marvin)
}

function checkIfSaveEnabled(){
    let demo_mode = $("#js-demo").val();
    let tutorial_mode = $("#js-tutorial").val();
    // returns true if save is enabled, else false
    return demo_mode !== 'demo' && tutorial_mode !== 'yes';
}



function sketcherLockHandler(){
    if (ifCurrentUserIsNotCreator()) {
        $("#marvin-test").css("pointer-events", "none");
    }
    if ($("#js-complete").val() === "complete") {
        $("#marvin-test").css("pointer-events", "none");
    }
}

function demo(marvin) { //Add the example of chemical reaction in the molecular editor
    // let smiles = "CC.OC(=O)C1=CC=CC=C1>>CCNC(=O)C1=CC=CC=C1.CC";
    let smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1";
    // let smiles = "CC(O)CN>>O=[N+]([O-])c1ccc(Cl)c([N+](=O)[O-])c1";
    marvin.importStructure("cxsmiles", smiles);
}

function reloadHandler(marvin){
    let $reactionSmiles = $("#js-reaction-smiles")
    let reloadedReaction = $reactionSmiles.val();
    marvin.importStructure("cxsmiles", reloadedReaction);
    $("#js-load-status").val("loaded")
    // show compound if reloading a reaction where the reaction table has previously been loaded and description is in json
    let js_reaction_table_data = JSON.parse(document.getElementById('js-reaction-table-data').value)
    if (Object.keys(js_reaction_table_data).includes("reaction_description")) {
        setTimeout(function(){
                fillReactionTable(marvin)}, 500)
    }
}

function autosaveInitialise(marvin){
    // autosave when the sketcher is edited
    marvin.on('molchange', function () {
        // save when the sketcher is edited and set sketcher arg to true to save just smiles
        setTimeout(autoSaveCheck(marvin, null, true), 500)
    });
}

function sketcherDataLossHandler(){
    const reactionDiv = document.getElementById("successAlert1")
    const firstChildNode = reactionDiv.childNodes[0]
    let tutorial_mode = $("#js-tutorial").val();
    if (tutorial_mode !== 'yes') {
        if (reactionDiv.childNodes.length !== 0 && firstChildNode.id !== "js-novel-compound-input-form") {
            const text = "Please note that any reaction data already inputted will be lost. Do you wish to continue?";
            if (confirm(text) === false) {
                return;
            }
        }
    }
}

async function fillReactionTable(marvin) {
    // populates reaction table
    sketcherDataLossHandler()
    exportImageString(marvin)
    let smiles = await exportSmilesFromSketcher(marvin)
    let [reactants, products] = processSmiles(smiles)
    // loading
    $(".loading-bar").css("display", "block");
    let smilesNew = removeReagentsFromSmiles(smiles)
    $("#js-reaction-smiles").val(smilesNew);
    let workgroup = $("#js-active-workgroup").val();
    let workbook = $("#js-active-workbook").val();
    let reaction_id = $("#js-reaction_id").val();
    let demo_mode = $("#js-demo").val();
    // Asynchronous request to _process in routes.py
    $.ajax({
        type: 'GET',
        contentType: 'application/json;charset-utf-08',
        dataType: 'json',
        url: '/_process?reactants=' + reactants + '&products=' + products + '&demo=' + demo_mode +
            '&workgroup=' + workgroup + '&workbook=' + workbook + "&reaction_id=" + reaction_id
    })
        .done(function (data) {
            $(".loading-bar").css("display", "none");
            if (data.novelCompound) {
                $('#successAlert1').html(data.reactionTable).show(); // Sends data to the reaction table
                sketcherNovelCompoundInputValidate()
            } else if (data.error) {
                $('#errorAlert').text(data.error).show();
                $('#successAlert1').hide();
            } else if (data.reactionTable === "Demo") {
                alert("One or more reactants or products are not in the database. Adding novel compounds is not available in Demo mode.")
            } else {
                // we can now generate summary table with data provided
                $('#successAlert1').html(data.reactionTable).show(); // Sends data to the reaction table
                initialiseReactionTable()
                $("#action-summary").show()
                $("#js-reaction-name2").val($("#js-reaction-name").val());
                // initial mass and equivalents to highlight
                let limitingReactantTableNumber = $("input[name='reactant-limiting']:checked").val();
                let colorRoundedReactantMassID = "#js-reactant-rounded-mass" + limitingReactantTableNumber;
                $("#js-reactant-equivalent" + limitingReactantTableNumber).val(1)
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
        });
    // when ajax is done
    document.body.className = "";
}


function sketcherEventListeners(marvin){
    document.getElementById("action-button-submit").addEventListener("click", function() {
        fillReactionTable(marvin)
    });
    document.getElementById("demo-button").addEventListener("click",function(){
        demo(marvin)
    });
}
