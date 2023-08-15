$(function() {
    if ($("#hazard_disclaimer").is(':checked') && $("#privacy").is(':checked')) {
        $("#submit").removeAttr('disabled');
    } else{
        $("#submit").prop('disabled', true);
    }
    if ($("#hazard_disclaimer").is(':checked')){
        $("#hazard_disclaimer").removeAttr('disabled')
    }
    if ($("#privacy").is(':checked')){
        $("#privacy").removeAttr('disabled')
    }

});

// enable the checkboxes once the user has clicked on the hyperlink to the relevant information
function enableHazardCheckbox(){
    $("#hazard_disclaimer").removeAttr('disabled');
}
function enablePrivacyCheckbox(){
    $("#privacy").removeAttr('disabled');
}
// enable the submit button once both checkboxes are enabled.
function enableSubmit(){
    if ($("#hazard_disclaimer").is(':checked') && $("#privacy").is(':checked')){
        $("#submit").removeAttr('disabled');
    } else {
         $("#submit").prop('disabled', true);
    }
}

// reaction test_reaction_reload; test_save_reaction; test_complete_reaction_locked; test_add_new_reagent_cas; test_solvent_selection