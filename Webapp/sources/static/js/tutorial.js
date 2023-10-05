$(function() {
    const tutorial = $("#js-tutorial").val()
    if (tutorial === "yes") {
        // disable page and scrolling
        setTimeout(() => {
            window.scrollTo(0, 0);
        }, "200")
        $('#page-contents').css('pointer-events', 'none');
        // show tutorial text
        document.getElementById("tutorial-1").style.display = "block";
    }
});

function tutorialNext(number){
    document.getElementById("tutorial-" + number).style.display = "none";
    number = number + 1
    document.getElementById("tutorial-" + number).style.display = "block";
}

function tutorialBack(number){
    document.getElementById("tutorial-" + number).style.display = "none";
    number = number - 1
    document.getElementById("tutorial-" + number).style.display = "block";
}

function tutorial_2(){
    window.scrollTo(0, 0);
    $('#name-id').removeClass('highlight');
}

function tutorial_3(){
    window.scrollTo(0, 0);
    $('#name-id').addClass('highlight');
    $('#sketchers-div').removeClass('highlight');
}

function tutorial_4(){
    window.scrollTo(0, document.body.scrollHeight);
    $('#name-id').removeClass('highlight');
    $('#editor-instructions').removeClass('highlight');
    $('#sketchers-div').addClass('highlight');
}

function tutorial_5(){
    window.scrollTo(0, 0);
    $('#sketchers-div').removeClass('highlight');
    $('#editor-instructions').addClass('highlight');
}

function tutorial_6(){
    window.scrollTo(0, 0);
    $('#action-button-submit').removeClass('highlight');
    $('#editor-instructions').addClass('highlight');
    // remove demo
    document.getElementById("clear-button").click();
}

function tutorial_7(){
    $('#reaction-table-div').html("").hide();
    document.getElementById("summary-div").style.display = "none";
    document.getElementById("demo-button").click();
    $('#editor-instructions').removeClass('highlight');
    $('#action-button-submit').addClass('highlight');
    $('#reaction-name-description').removeClass('highlight');
    window.scrollTo(0, document.body.scrollHeight);
}

function tutorial_8(){
    document.getElementById("summary-div").style.display = "block";
    document.getElementById("action-button-submit").click();
    setTimeout(() => {
        window.scrollTo(0, document.body.scrollHeight);
        $('#action-button-submit').removeClass('highlight');
        $('#reaction-name-description').addClass('highlight');
        }, "500")
}

function tutorial_9(){
    $('#reaction-name-description').removeClass('highlight');
    $('#reaction-table-div').addClass('highlight');
}

function tutorial_10(){
    $('#js-reactant-rounded-mass1').val("").trigger('change');
    $('#js-reactant-equivalent2').val("").trigger('change');
    $("#js-reactant-physical-form1").prop('selectedIndex', 0).trigger('change');
    $("#js-reactant-physical-form2").prop('selectedIndex', 0).trigger('change');
    $("#js-product-physical-form1").prop('selectedIndex', 0).trigger('change');
    $("#js-reactant-rounded-amount1").val("-")
    $("#js-reactant-rounded-volume1").val("-")
    $("#js-reactant-rounded-amount2").val("-")
    $("#js-reactant-rounded-volume2").val("-")
    $("#js-reactant-rounded-mass2").val("-")
    $("#js-product-rounded-amount1").val("-")
    $("#js-product-rounded-mass1").val("-")
    $('#reagent-row').removeClass('highlight');
    $('#reaction-table-div').addClass('highlight');
}

function tutorial_11(){
    $('#js-reactant-rounded-mass1').val(3).trigger('change');
    $('#js-reactant-equivalent2').val(1).trigger('change');
    $("#js-reactant-physical-form1").prop('selectedIndex', 4).trigger('change');
    $("#js-reactant-physical-form2").prop('selectedIndex', 8).trigger('change');
    $("#js-product-physical-form1").prop('selectedIndex', 4).trigger('change');
    $('#reaction-table-div').removeClass('highlight');
    $('#reagent-row').addClass('highlight');
    if ($("#js-number-of-reagents").val() !== "") {
        if ($("#js-number-of-reagents").val() !== "0") {
            document.getElementById("remove-reagent1").click();
        }
    }
}

function tutorial_12(){
    if ($("#js-number-of-solvents").val() !== ""){
        if ($("#js-number-of-solvents").val() !== "0") {
            document.getElementById("remove-solvent1").click();
        }
    }
    if ($("#js-number-of-reagents").val() !== "") {
        if ($("#js-number-of-reagents").val() !== "0") {
            document.getElementById("remove-reagent1").click();
        }
    }
    document.getElementsByClassName("js-add-reagent")[0].click();
    $('#reagent-row').removeClass('highlight');
    $('#js-reagent-table-row1').addClass('highlight');
}

function tutorial_13(){
    $('#js-reagent-table-row1').removeClass('highlight');
    $('#js-reagent1').val("4128-37-4").trigger('change');
    $('#js-reagent-equivalent1').val(1).trigger('change')
    $("#js-reagent-physical-form1").prop('selectedIndex', 3).trigger('change');
    if ($("#js-number-of-solvents").val() !== ""){
        if ($("#js-number-of-solvents").val() !== "0") {
            document.getElementById("remove-solvent1").click();
        }
    }
    document.getElementsByClassName("js-add-solvent")[0].click();
    document.getElementById("js-solvent1").click();
    $('#js-solvent-table-row1').addClass('highlight');
    setTimeout(() => {
        document.getElementById("js-reaction-name2").scrollIntoView();
        }, "100")
}

function tutorial_14(){
    $("#js-solvent1").val("Diethyl ether").trigger("keyup");
    $("#js-solvent1").css("border-color", "black");
    $("#js-solvent1").css("border-width", "thin");
    document.getElementById("js-solvent-datalist1").style.display = "none";
    $('#js-solvent-volume1').val("10").trigger('change');
    $("#js-solvent-physical-form1").prop('selectedIndex', 8).trigger('change');
    $('#action-summary').removeClass('highlight');
    $('#js-solvent-table-row1').addClass('highlight');
}

function tutorial_15(){
    document.getElementById("print-container").style.display = "none";
    document.getElementById("complete-reaction-div").style.display = "none";
    document.getElementById("reaction-file-attachments").style.display = "none";
    document.getElementById("print-pdf").style.display = "none";
    window.scrollTo(0, document.body.scrollHeight);
    $('#js-solvent-table-row1').removeClass('highlight');
    $('#action-summary').addClass('highlight');
}

function tutorial_16(){
    document.getElementById("print-container").style.display = "block";
    document.getElementById("action-summary").click();
    $('#action-summary').removeClass('highlight');
    setTimeout(() => {
        document.getElementById("action-summary").scrollIntoView();
        flagElementSustainability();
        document.getElementById("complete-reaction-div").style.display = "block";
        document.getElementById("print-pdf").style.display = "block";
        }, "200")
}

function tutorial_17(){
    document.getElementById("summary-product-cell").scrollIntoView();
    $('#js-temperature').val("80").trigger('input');
    $("#js-batch-flow").prop('selectedIndex', 1).trigger('input');
    $("#js-isolation").prop('selectedIndex', 1).trigger('input');
    $("#js-catalyst").prop('selectedIndex', 1).trigger('input');
    $('#js-unreacted-reactant-mass').val("").trigger('input');
    $('#js-real-product-mass').val("").trigger('input');
}

function tutorial_18(){
    document.getElementById("js-product1").scrollIntoView();
    $('#js-unreacted-reactant-mass').val("0.5").trigger('input');
    $('#js-real-product-mass').val("3").trigger('input');
}

function tutorial_19(){
    window.scrollTo(0, document.body.scrollHeight);
}

function tutorial_20(){
    window.scrollTo(0, 0);
}