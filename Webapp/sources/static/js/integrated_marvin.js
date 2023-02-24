// setup MarvinJS sketcher window
$(document).ready(function () {
    ChemicalizeMarvinJs.createEditor("#marvin-test").then(function (marvin) {
        marvin.setDisplaySettings({
            "toolbars": "education",
        });
        // if not creator lock editor
        if (ifCurrentUserIsNotCreator()){
            $("#marvin-test").css("pointer-events", "none");
        }
        if ($("#js-complete").val() === "complete"){
            $("#marvin-test").css("pointer-events", "none");
        }
        // reload reaction smiles if present
        let demo_mode = $("#js-demo").val();
        let tutorial_mode = $("#js-tutorial").val();
        if (demo_mode !== 'demo' && tutorial_mode !== 'yes'){
            let $reactionSmiles = $("#js-reaction-smiles")
            let reloadedReaction = $reactionSmiles.val();
            marvin.importStructure("cxsmiles", reloadedReaction);
            $("#js-load-status").val("loaded")
             // show compound if reloading a reaction where the reaction table has previously been loaded and description is in json
            let js_reaction_table_data = JSON.parse(document.getElementById('js-reaction-table-data').value)
            if (Object.keys(js_reaction_table_data).includes("reaction_description")) {
                setTimeout(showCompound, 500)
             }
            // autosave when the sketcher is edited
            marvin.on('molchange', function () {
            // save when the sketcher is edited and set sketcher arg to true to save just smiles
            setTimeout(autoSaveCheck(null, true), 500)
            });
        }
        //Sends a reaction from the molecular editor to the table
        async function showCompound() {
            // check if reaction table already populated if so warn user
            const reactionDiv = document.getElementById("successAlert1")
            const firstChildNode = reactionDiv.childNodes[0]
            if (tutorial_mode !== 'yes'){
                if (reactionDiv.childNodes.length !== 0 && firstChildNode.id !== "js-novel-compound-input-form") {
                    const text = "Please note that any reaction data already inputted will be lost. Do you wish to continue?";
                    if (confirm(text) === false) {
                        return;
                    }
                }
            }
            storeImageString(marvin)
            let smiles = await exportSmilesFromSketcher()
            let [reactants, products] = processSmiles(smiles)
            // loading
            $(".loading-bar").css("display", "block");
            let smilesNew = removeReagentsFromSmiles(smiles)
            $("#js-reaction-smiles").val(smilesNew);
            let workgroup = $("#js-active-workgroup").val();
            let workbook = $("#js-active-workbook").val();
            let reaction_id = $("#js-reaction_id").val();
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

        // functions that use marvin
        function demo() { //Add the example of chemical reaction in the molecular editor
            // let smiles = "CC.OC(=O)C1=CC=CC=C1>>CCNC(=O)C1=CC=CC=C1.CC";
            let smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1";
            // let smiles = "CC(O)CN>>O=[N+]([O-])c1ccc(Cl)c([N+](=O)[O-])c1";
            marvin.importStructure("cxsmiles", smiles);
        }

        function exportSmilesFromSketcher() {
            return marvin.exportStructure('cxsmiles', {'extra': 'f'}).then(function (smiles) {
                return smiles
            });
        }

       window.sketcherAutoSave = async function sketcherAutoSave() {
            // autosave for when the sketcher has been updated
            let smiles = await exportSmilesFromSketcher()
            let smilesNew = removeReagentsFromSmiles(smiles)
            $("#js-reaction-smiles").val(smilesNew);
            let workgroup = $("#js-active-workgroup").val();
            let workbook = $("#js-active-workbook").val();
            let reactionID = $("#js-reaction-id").val();
            let userEmail = "{{ current_user.email }}";
            $.ajax({
                url: '/_autosave_sketcher',
                type: 'post',
                data: {
                    workgroup: workgroup,
                    workbook: workbook,
                    reactionID: reactionID,
                    userEmail: userEmail,
                    reactionSmiles: smilesNew
                    },
                dataType: "json",
                success: function (){
                    flashUserSaveMessage()
                },
                error: function (){
                    flashUserErrorSavingMessage()
                }
            });
        }
    document.getElementById("action-button-submit").addEventListener("click", showCompound);
    document.getElementById("demo-button").addEventListener("click", demo);
    });
});


function storeImageString(marvin) {
    // saves image to a hidden html input
    let settings = {'width': 600, 'height': 400};
    marvin.exportStructure("jpeg", settings).then(function (source) {
        $("#js-reloaded-image").val(source);
    });
}


function processSmiles(smiles) {
    let reaction = smiles.split(" |")[0];
    //replace the plus, minus and sharp signs from a smiles string to trasfer it properly
    reaction = reaction.replace(/\+/g, 'plus');
    reaction = reaction.replace(/-/g, 'minus');
    reaction = reaction.replace(/#/g, 'sharp');
    //split into reactants, solvents and product
    let array = reaction.split(">");
    let reactants = array[0].split(".");
    let rl = reactants.length;
    let products = array[2].split(".");
    let reagents = array[1].split(".");
    let compounds, rgl, ending, group, element, c0, c, r, rg;
    let i, k, p;
    if (reagents[0] !== "") {
        rgl = reagents.length;
        compounds = reactants.concat(reagents, products);
    } else {
        rgl = 0;
        compounds = reactants.concat(products);
    }
    // if "f:" is at the end of the smiles string, it contains an ionic compound
    if (smiles.includes("f:") == true) {
        // split after the comma to see how many ionic compounds in reaction smiles
        let ionList = smiles.split("f:")[1].split('|')[0].split(',');
        // for each ionic compound, join all ions to make the full salt and replace ions in reactant/product list
        let productSaltDictList = []
        let reactantSaltDictList = []
        for (let adjacentIons of ionList) {
            // make list of index, and smiles, then join to make a string of the full salt from the ions
            let ionIndexList = adjacentIons.split('.').map(Number);
            let ionCompoundsSmiles = [] // [Pd+, OAc-, OAc-]
            for (let idx of ionIndexList) {
                ionCompoundsSmiles.push(compounds[idx])
            }
            let saltSmiles = ionCompoundsSmiles.join(".")
            // determine if reactant or product
            let rxnComponent
            if (rl > ionIndexList[0]) {
                rxnComponent = "reactant"
            } else {
                rxnComponent = "product"
            }
            // ions in previous salts for only same rxn component (reactant or product) - to determine insert position
            // when splicing into reactant/product list.
            let ionsInPreviousSalts = 0
            let insertPosition
            if (rxnComponent === 'reactant') {
                for (let previousSaltDict of reactantSaltDictList) {
                    ionsInPreviousSalts += previousSaltDict["numberOfIons"]
                }
                insertPosition = ionIndexList[0] - ionsInPreviousSalts + reactantSaltDictList.length
                reactantSaltDictList.push({
                    "saltSmiles": saltSmiles,
                    "numberOfIons": ionCompoundsSmiles.length,
                    "insertPosition": insertPosition
                })
            }
            if (rxnComponent === "product") {
                for (let previousSaltDict of productSaltDictList) {
                    ionsInPreviousSalts += previousSaltDict["numberOfIons"]
                }
                insertPosition = ionIndexList[0] - ionsInPreviousSalts + productSaltDictList.length - rl
                productSaltDictList.push({
                    "saltSmiles": saltSmiles,
                    "numberOfIons": ionCompoundsSmiles.length,
                    "insertPosition": insertPosition
                })
            }
        }
        // splice the salts into the lists replacing the ions that they contain
        for (let saltDict of reactantSaltDictList) {
            reactants.splice(saltDict["insertPosition"], saltDict["numberOfIons"], saltDict["saltSmiles"])
        }
        for (let saltDict of productSaltDictList){
                products.splice(saltDict["insertPosition"], saltDict["numberOfIons"], saltDict["saltSmiles"])

        }
        } else {
            ending = "none";
    }
    return [reactants, products]
}

function removeReagentsFromSmiles(smiles){
    // remove reagents from reaction smiles
    let smiles2 = smiles.split('>').slice(-1);
    let smiles1 = smiles.split('>')[0];
    return smiles1 + '>>' + smiles2;
}