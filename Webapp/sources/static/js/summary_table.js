// Functions to reload summary & control demo mode.
$(document).ready(function() {initialise()});

async function initialise(){
    checkDemoMode()
    setColours()
    await reloadSummary()
    $("#js-load-status").val('loaded');
}

function checkDemoMode(){
    if ($("#js-demo").val() === "demo"){
        document.getElementById("save-reaction-div").style.display = "none";
        document.getElementById("mark-complete-div").style.display = "none";
    }
}

function reloadSummary() {
    if ($("#js-tutorial").val() === "yes" || $("#js-demo").val() === 'demo') {
        return
    }
    let js_summary_table_data = JSON.parse(document.getElementById("js-summary-table-data").value)
    if (js_summary_table_data === "no data") {
        $("#js-load-status").val('loaded')
        return
    }
    // real product mass
    $("#js-real-product-mass").val(js_summary_table_data["real_product_mass"]);
    if (js_summary_table_data["real_product_mass"] !== "") {
        calculateRealProductMass();
        let massEfficiency = js_summary_table_data["mass_efficiency"];
        $("#js-me").val(massEfficiency);
        let meFlag;
        if (massEfficiency < 70) {
            meFlag = "hazard-hazardous";
        } else if (massEfficiency >= 90) {
            meFlag = "hazard-acceptable";
        } else {
            meFlag = "hazard-warning";
        }
        $("#js-me").attr('class', meFlag);
        $("#js-me-cell").attr('class', meFlag);
        setColours();
    }
    // unreacted reactant mass
    $("#js-unreacted-reactant-mass").val(js_summary_table_data["unreacted_reactant_mass"]);
    if (js_summary_table_data["unreacted_reactant_mass"] !== ""){
        // reload conversion
        flagConversion();
        // refill selectivity
        flagSelectivity();
    }
    // reaction temperature
    $("#js-temperature").val(js_summary_table_data["reaction_temperature"]);
    if (js_summary_table_data["reaction_temperature"] !== "") {
        flagTemperature()
    }
    // reaction elements
    // property is auto calculated - use autocalculation
    if (js_summary_table_data["element_sustainability"] !== 'undefined') {
        $("#js-elements").prop("selectedIndex", js_summary_table_data["element_sustainability"]);
        if (js_summary_table_data["element_sustainability"] !== "") {
            flagElementSustainability();
        }
    }
    // batch or flow
    $("#js-batch-flow").val(js_summary_table_data["batch_flow"]);
    if (js_summary_table_data["batch_flow"] !== "") {
        flagBatchFlow();
    }
    // isolation method
    $("#js-isolation").prop("selectedIndex", js_summary_table_data["isolation_method"]);
    if (js_summary_table_data["isolation_method"] !== "") {
        flagIsolation();
    }
    // catalyst
    $("#js-catalyst").val(js_summary_table_data["catalyst_used"]);
    if (js_summary_table_data["catalyst_used"] !== "") {
        flagCatalyst();
    }
    // recovery
    $("#js-recovery").val(js_summary_table_data["catalyst_recovered"]);
    if (js_summary_table_data["catalyst_recovered"] !== "") {
        flagRecovery();
    }
    // custom protocols
    $("#field1-text").val(js_summary_table_data["custom_protocol1"]);
    $("#field2-text").val(js_summary_table_data["custom_protocol2"]);
    // radio buttons
    for (let i = 0; i < js_summary_table_data["radio_buttons"].length; i++) {
        let radio_name = js_summary_table_data["radio_buttons"][i];
        document.getElementById(radio_name).checked = true;
        calculateRiskScore();
    }
    // other hazard text
    $("#other-risks-textbox").val(js_summary_table_data["other_hazards_text"]);
    // researcher
    $("#js-researcher").val(js_summary_table_data["researcher"]);
    // supervisor
    $("#js-supervisor").val(js_summary_table_data["supervisor"]);
    // refill risk score
    autofillRiskScore();
    // get reaction image
    exportImage();
    // setTimeout(exportImage, 2000);
    $("#complete-reaction-div").show();
    $("#print-pdf").show();
    return new Promise(function (resolve) {
        resolve("reload complete")
    });

}

function calculateRealProductMass(){
    let realProductMass = $("#js-real-product-mass").val();
    let mainTheoreticalProductMass = $("#js-main-product-mass").val();
    let percentageYield = 100*realProductMass/mainTheoreticalProductMass;
    percentageYield = percentageYield.toFixed();
    $("#js-percentage-yield").val(percentageYield);
    $("#js-yield").val(percentageYield);
    let yieldFlag
    if (percentageYield < 70){
        yieldFlag = "hazard-hazardous";
    } else if (percentageYield >= 90){
        yieldFlag = "hazard-acceptable";
    } else {
        yieldFlag = "hazard-warning";
    }
    $("#js-yield").attr('class', yieldFlag);
    $("#js-yield-cell").attr('class', yieldFlag);
    setColours();
}
function flagTemperature(){
    let temperature = $("#js-temperature").val();
    let temperatureFlag
    if (temperature >= 0 && temperature < 70){
        temperatureFlag = "hazard-acceptable";
    } else if (temperature < -20 || temperature > 140){
        temperatureFlag = "hazard-hazardous";
    } else {
        temperatureFlag = "hazard-warning"
    }
    $("#js-temperature").attr('class', temperatureFlag);
    $("#js-temperature-cell").attr('class', temperatureFlag);
    setColours();
}
function flagElementSustainability(){
    let elements = $("#js-elements").val();
    let elementsFlag
    if (elements === "5-50 years"){
        elementsFlag = "hazard-hazardous";
    } else if (elements === "50-500 years"){
        elementsFlag = "hazard-warning";
    } else if (elements === "+500 years"){
        elementsFlag = "hazard-acceptable";
    } else {
        elementsFlag = "hazard-reset-hazard";
    }
    $("#js-elements").attr('class', elementsFlag);
    $("#js-elements-cell").attr('class', elementsFlag);
    setColours();
}
function flagBatchFlow(){
    let batchFlow = $("#js-batch-flow").val();
    let batchFlowFlag
    if (batchFlow === "Batch"){
        batchFlowFlag = "hazard-warning";
    } else if (batchFlow === "Flow"){
        batchFlowFlag = "hazard-acceptable";
    } else {
        batchFlowFlag = "hazard-reset-hazard";
    }
    $("#js-batch-flow").attr('class', batchFlowFlag);
    $("#js-batch-flow-cell").attr('class', batchFlowFlag);
    setColours();
}
function flagIsolation(){
    let isolation = $("#js-isolation").val();
    let isolationFlag
    if (["Crystallization", "Filtration", "Distillation < 140℃"].includes(isolation)){
        isolationFlag = "hazard-acceptable";
    } else if (["Column", "HPLC", "Ion exchange", "Multiple recryst.", "Distillation > 140℃"].includes(isolation)){
        isolationFlag = "hazard-hazardous";
    } else {
        isolationFlag = "hazard-reset-hazard";
    }
    $("#js-isolation").attr('class', isolationFlag);
    $("#js-isolation-cell").attr('class', isolationFlag);
    setColours();
}
function flagCatalyst(){
    let catalyst = $("#js-catalyst").val();
    let catalystFlag
    if (catalyst === "No catalyst" || catalyst === "Catalyst or enzyme"){
        catalystFlag = "hazard-acceptable";
    } else if (catalyst === "Stoichiometric reagents"){
        catalystFlag = "hazard-warning";
    } else if (catalyst === "Excess reagents"){
        catalystFlag = "hazard-hazardous";
    } else {
        catalystFlag = "hazard-reset-hazard";
    }
    $("#js-catalyst").attr('class', catalystFlag);
    $("#js-catalyst-cell").attr('class', catalystFlag);
    setColours();
}
function flagRecovery(){
    let recovery = $("#js-recovery").val();
    let recoveryFlag
    if (recovery === "Recovered catalyst"){
        recoveryFlag = "hazard-acceptable";
    } else if (recovery === "Not recovered catalyst"){
        recoveryFlag = "hazard-warning";
    } else {
        recoveryFlag = "hazard-reset-hazard";
    }
    $("#js-recovery").attr('class', recoveryFlag);
    $("#js-recovery-cell").attr('class', recoveryFlag);
    setColours();
}
//Autofill Risk Score
function calculateRiskScore(){
        let hazard = $(".js-hazard:checked").val();
        let risk = $(".js-risk:checked").val();
        let consequences = $(".js-consequences:checked").val();
        if(typeof hazard === 'undefined'){
            return;
        }
        if(typeof risk === 'undefined'){
            return;
        }
        if(typeof consequences === 'undefined'){
            return;
        }
        let riskScore = hazard*risk*consequences;
        $("#js-risk-score").val(riskScore);
        if (riskScore < 3) {
            $("#categoryD").prop("checked", true);
        } else if (riskScore > 2 && riskScore < 6) {
            $("#categoryC").prop("checked", true);
        } else if (riskScore > 5 && riskScore < 10) {
            $("#categoryB").prop("checked", true);
        } else if (riskScore > 9 && riskScore < 28) {
            $("#categoryA").prop("checked", true);
        }
}

//Autofill conversion
function flagConversion() {
    let unreactedReactantMass = $("#js-unreacted-reactant-mass").val();
    let limitedReactantMass = $("#js-limited-reactant-mass").val();
    let conversion = (1 - unreactedReactantMass / limitedReactantMass) * 100;
    let conversionFlag
    if (conversion < 70) {
        conversionFlag = "hazard-hazardous";
    } else if (conversion >= 90) {
        conversionFlag = "hazard-acceptable";
    } else {
        conversionFlag = "hazard-warning";
    }
    conversion = conversion.toFixed();
    $("#js-conversion").val(conversion);
    $("#js-conversion").attr('class', conversionFlag);
    $("#js-conversion-cell").attr('class', conversionFlag);
    setColours();
}

function flagSelectivity() {
    let conversion = $("#js-conversion").val();
    let percentageYield = $("#js-yield").val();
    if (!conversion || !percentageYield){
        $("#js-selectivity").val("");
        $("#js-selectivity").attr('class',"hazard-reset-hazard");
        $("#js-selectivity-cell").attr('class',"hazard-reset-hazard");
    }
    else{
        let selectivity = 100*percentageYield/conversion;
        let selectivityFlag
        if (selectivity < 70){
            selectivityFlag = "hazard-hazardous";
        } else if (selectivity >= 90){
            selectivityFlag = "hazard-acceptable";
        } else {
            selectivityFlag = "hazard-warning";
        }
        selectivity = selectivity.toFixed();
        $("#js-selectivity").val(selectivity);
        $("#js-selectivity").attr('class',selectivityFlag);
        $("#js-selectivity-cell").attr('class',selectivityFlag);
        setColours();
    }
}

function flagMassEfficiency(){
// sum of reactants and reagents masses / actual product mass
        let massFactor = {'g': 1, 'mg': 10 ** (-3), 'μg': 10 ** (-6)};
        let reactantReagentMassUnit = $("#js-mass-unit").val();
        let reactantReagentMassFactor = massFactor[reactantReagentMassUnit]

        let numberOfReactants = $("#js-number-of-reactants").val();
        let massSum = 0;
        for (let i=1; i <= numberOfReactants; i++){
            // Mass of each reactant
            let reactantMassID = "#js-reactant-mass" + i;
            let reactantMass = parseFloat($(reactantMassID).val());
            massSum += reactantMass;
        }
        let numberOfReagents = parseInt($("#js-number-of-reagents").val());
        for (let i=1; i <= numberOfReagents; i++){
            // Mass of each reagent
            let reagentMassID = "#js-reagent-mass" + i;
            let reagentMass = parseFloat($(reagentMassID).val());
            massSum += reagentMass;
        }

        let productMass = parseFloat($("#js-real-product-mass").val());
        let productUnit = $("#js-product-mass-unit").val();
        let productMassFactor = massFactor[productUnit]
        let massEfficiency = 100 * (productMass * productMassFactor) / (massSum * reactantReagentMassFactor);
        let meFlag;
        if (massEfficiency < 70) {
            meFlag = "hazard-hazardous";
        } else if (massEfficiency >= 90) {
            meFlag = "hazard-acceptable";
        } else {
            meFlag = "hazard-warning";
        }
        massEfficiency = massEfficiency.toFixed();
        $("#js-me").val(massEfficiency);
        $("#js-me").attr('class', meFlag);
        $("#js-me-cell").attr('class', meFlag);
        setColours();
}



