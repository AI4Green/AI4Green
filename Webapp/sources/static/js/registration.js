// enable the submit button once both checkboxes are enabled.
function enableSubmit() {
  if ($("#hazard").is(":checked") && $("#privacy").is(":checked")) {
    $("#submit").removeAttr("disabled");
  } else {
    $("#submit").prop("disabled", true);
  }
}

// reaction test_reaction_reload; test_save_reaction; test_complete_reaction_locked; test_add_new_reagent_cas; test_solvent_selection
