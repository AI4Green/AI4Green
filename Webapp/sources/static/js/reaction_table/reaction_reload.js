async function reactionTableReload(){
    /*
    Reloads the data to the reaction table and then calls the function to reload summary if summary data is present
     */
    // assign every value to its saved value
    if ($("#js-demo").val() === 'demo' || $("#js-tutorial").val() === 'yes'){
        return
    }
    let js_reaction_table_data = JSON.parse(document.getElementById('js-reaction-table-data').value)
    if (js_reaction_table_data === "no data") {
        return
    }
    let $loadTracker = $("#js-load-status")
    $loadTracker.val("loading")
    // limiting reactant table number and main product
    var limiting_name = "#js-reactant-limiting" + js_reaction_table_data["limiting_reactant_table_number"];
    $(limiting_name).prop("checked", true);
    var main_name = "#js-main-product" + js_reaction_table_data["main_product"];
    $(main_name).prop("checked", true);
    let numberOfReactants = Number($("#js-number-of-reactants").val())
    let numberOfProducts = Number($("#js-number-of-products").val())
    // define lists
    const styleElements = ["reactant-rounded-mass", "reactant-equivalent", "reactant-physical-form",
                        "reagent-equivalent", "reagent-physical-form",
                        "solvent-name", "solvent-volume", "solvent-physical-form", "product-physical-form"]
    // Field dicts to relate HTML IDs to js_reactant_table_data keys
    const productFields = {
        "amounts": "-rounded-amount", "amounts_raw": "-amount", "masses": "-rounded-mass",
        "masses_raw": "-mass", "physical_forms": "-physical-form",
    }
    // object.assign concatenates the two arguments - all productFields are also found in reactant and reagent
    const reactantFields =  Object.assign({"equivalents": "-equivalent", "concentrations": "-concentration",
        "densities": "-density", "volumes_raw": "-volume", "volumes": "-rounded-volume"}, productFields)

    const reagentFields = Object.assign({"names": "", "equivalents": "-equivalent", "concentrations": "-concentration",
        "densities": "-density", "volumes_raw": "volume-", "volumes": "-rounded-volume",
        "molecular_weights": "-molecular-weight", "hazards": "-hazards"}, productFields)

    const solventFields = {"volumes": "-volume", "concentrations": "-rounded-concentration",
         "names": "", "hazards": "-hazards", "physical_forms": "-physical-form"}

    // iterate through and reload reactants data
    for (let i=0; i < numberOfReactants; i++){
        let j = String(i + 1)
        // iterate through the dictionary for each reactant field
        for (const [jsonField, fieldID] of Object.entries(reactantFields)){
            let reactantFieldID = "#js-reactant" + fieldID + j
            let jsonID = "reactant_" + jsonField
            fillInputField(reactantFieldID, jsonID, i)
        }
    }
    // iterate through and reload products data
    for (let i=0; i < numberOfProducts; i++) {
        let j = String(i + 1)
        // iterate through the dictionary for each product field
        for (const [jsonField, fieldID] of Object.entries(productFields)) {
            let productFieldID = "#js-product" + fieldID + j
            let jsonID = "product_" + jsonField
            fillInputField(productFieldID, jsonID, i)
        }
    }
    // iterate through and reload reagents data
    let numberOfReagents = js_reaction_table_data["reagent_names"].length
    for (let i=0; i < numberOfReagents; i++) {
        let j = String(i + 1)
        document.querySelector('.js-add-reagent').click();
        // wait for reagent data to be filled in from database before proceeding
        await postReagentData(js_reaction_table_data["reagent_names"][i], i + 1)
        // iterate through the dictionary for each reagent field
        for (const [jsonField, fieldID] of Object.entries(reagentFields)) {
            let reagentFieldID = "#js-reagent" + fieldID + j
            let jsonID = "reagent_" + jsonField
            fillInputField(reagentFieldID, jsonID, i)
        }
        updateProductTableNumber();
    }
    // iterate through and reload solvents data
    let numberOfSolvents = js_reaction_table_data["solvent_names"].length
    for (let i=0; i < numberOfSolvents; i++) {
        let j = String(i + 1)
        document.querySelector('.js-add-solvent').click();
        // wait for solvent data to be filled in from database before proceeding
        await postSolventData(js_reaction_table_data["solvent_names"][i], i + 1)
        // iterate through the dictionary for each solvent field
        for (const [jsonField, fieldID] of Object.entries(solventFields)) {
            let solventFieldID = "#js-solvent" + fieldID + j
            let jsonID = "solvent_" + jsonField
            fillInputField(solventFieldID, jsonID, i)
        }
        let numberOfReagents = Number($("#js-number-of-reagents").val());
        let numberOfReactants = Number($("#js-number-of-reactants").val());
        let solventTableNumberID = "#js-solvent-table-number" + j;
        let solventTableNumber = numberOfReactants + numberOfReagents + Number(j);
        $(solventTableNumberID).val(solventTableNumber);
        updateProductTableNumber();
    }
    // reload units
    const unitFields = Object.assign({"mass_units": "js-mass-unit", "amount_units": "js-amount-unit",
        "product_amount_units": "js-product-amount-unit", "solvent_volume_units": "js-solvent-volume-unit",
        "volume_units": "js-volume-unit", "product_mass_units": "js-product-mass-unit"})
    for (const [jsonField, fieldID] of Object.entries(unitFields)) {
        $("#" + fieldID).val(js_reaction_table_data[jsonField])
    }
    // reaction name
    $("#js-reaction-name2").val($("#js-reaction-name").val());
    // reaction description
    $("#js-reaction-description").val(js_reaction_table_data["reaction_description"]);
    // get summary table
    let js_summary_table_data = JSON.parse($("#js-summary-table-data").val());
    // load summary table if it has previously been loaded, element sustainability is used because this is autofilled upon load.
    if (js_summary_table_data["element_sustainability"] !== 'undefined'){
            setTimeout(showSummary(), 1000);
        } else{
        $loadTracker.val("loaded")
    }
    //  disable editing of the reaction if not owner
    if (ifCurrentUserIsNotCreator()){
        controlNonCreatorFunctionality()
    }
    // auxiliary reload functions that use the reaction table json
    function fillInputField(fieldID, jsonID, i){
        // either use dropdown or value depending if field is physical form or not
        if (fieldID.includes("physical-form")){
            $(fieldID).prop("selectedIndex", js_reaction_table_data[jsonID][i]);
        }
        else {
            $(fieldID).val(js_reaction_table_data[jsonID][i])
        }
        // checks if the field should be styled (highlights red if blank)
        if (styleElements.some(styleElement => fieldID.includes(styleElement))) {
            checkStyleElements(fieldID)
        }
    }
    function checkStyleElements(fieldID){
        // want to remove limiting reactant equiv and non-limiting reactant masses from styling required elements list
        if (fieldID.includes("reactant") || fieldID.includes("reagents")){
            let limitingReactantTableNumber = js_reaction_table_data["limiting_reactant_table_number"];
            if (fieldID.slice(-1) === limitingReactantTableNumber && fieldID.includes("equivalent") && fieldID.includes("reactant")){
                return
            }
            if (fieldID.slice(-1) !== limitingReactantTableNumber && fieldID.includes("rounded-mass")){
                return
            }
        } else {
            // most elements go here to apply styling function
            autoChangeRequiredStyling2(fieldID)
        }
    }
}

/**
 * When a user reloads a reaction they are not the creator of we disable inputs to prevent them editing the reaction
 * We enable the options to print a PDF summary and to download/view file attachments.
 */
function controlNonCreatorFunctionality(){
    $("#page-contents :input").prop("disabled", true);
    $("#print-pdf").prop("disabled", false);
    $("#file-list").find("*").prop("disabled", false);
    $(".delete-file-button").prop("disabled", true);
    $("#reaction-note-button").hide()
}