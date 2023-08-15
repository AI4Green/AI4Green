// Functions for autofilling & flagging the sustainability metrics. Many functions referenced found in summary_table.js
$(document).ready(function() {initialise()});

function initialise(){
    autofillYield()
    autoFillTemperature()
    autoFillElementSustainability()
    autofillBatchFlow()
    autoFillIsolation()
    autoFillCatalyst()
    autoFillRecovery()
    autofillRiskScoreObserver()
    autoFillConversion()
    autoFillSelectivityObserver()
    autofillMassEfficiency()
}

function autoFillSelectivityObserver(){
    autofillSelectivity("#js-unreacted-reactant-mass");
    autofillSelectivity("#js-real-product-mass");
}

function autofillRiskScoreObserver(){
    autofillRiskScore(".js-hazard");
    autofillRiskScore(".js-risk");
    autofillRiskScore(".js-consequences");
}
function autofillYield(){
    // autofills the yield & % yield
    $("#js-real-product-mass").on('input', function(){
        let realProductMass = $("#js-real-product-mass").val();
        if (!realProductMass){
            $("#js-percentage-yield").val("");
            $("#js-yield").val("");
            $("#js-yield").attr('class', "hazard-reset-hazard");
            $("#js-yield-cell").attr('class', "hazard-reset-hazard");
            setColours();
        }
        else {
            calculateRealProductMass()
        }
    });
}

function autoFillTemperature() {
    // autofill temperature flag
    $("#js-temperature").on('input', function () {
        flagTemperature()
    });
}
function autoFillElementSustainability(){
    // autofill element sustainability flag
    $("#js-elements").on('input',function(){
        flagElementSustainability()
    });
}

function autofillBatchFlow() {
    // autofill batch/flow flag
    $("#js-batch-flow").on('input', function () {
        flagBatchFlow()
    });
}

function autoFillIsolation(){
    $("#js-isolation").on('input', function () {
        flagIsolation()
    });
}
function autoFillCatalyst() {
    $("#js-catalyst").on('input', function () {
       flagCatalyst()
    });
}
function autoFillRecovery() {
    $("#js-recovery").on('input', function () {
        flagRecovery()
    });
}
//Autofill Risk Score
function autofillRiskScore(changedParameter){
    $(changedParameter).on('input',function(){
        calculateRiskScore()
    });
}

function autoFillConversion() {
//Autofill conversion
    $("#js-unreacted-reactant-mass").on('input', function () {
        let unreactantReactantMass = $("#js-unreacted-reactant-mass").val();
        if (!unreactantReactantMass) {
            $("#js-conversion").val("");
            $("#js-conversion").attr('class', "hazard-reset-hazard");
            $("#js-conversion-cell").attr('class', "hazard-reset-hazard");
            setColours();
        } else {
            flagConversion()
        }
    });
}
//Autofill selectivity
function autofillSelectivity(changedParameter){
    $(changedParameter).on('input',function(){
        flagSelectivity()
    });
}

function autofillMassEfficiency() {
    $("#js-real-product-mass").on('input', function () {
        // sum of reactants and reagents masses / actual product mass
        flagMassEfficiency()
    });
}
