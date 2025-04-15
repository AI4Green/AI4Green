// Functions for adding new compounds to the database - reagents and those drawn by sketcher
function novelCompoundDataFromSketcher() {
  // sends data to the backend routes
  let name = $("#js-new-compound-name").val();
  let molWeight = $("#js-new-compound-mw").val();
  let hPhrase = $("#js-new-compound-hazards").val();
  let cas = $("#js-new-compound-cas").val();
  let density = $("#js-new-compound-density").val();
  let concentration = $("#js-new-compound-concentration").val();
  let smiles = $("#js-new-compound-smiles").val();
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  // let polymerMode = $('input[id="polymer-mode-select"]').prop("checked");
  // let polymer = $("#js-polymer").val();
  let rxn = $("#js-reaction-rxn").val();
  let polymerIndices = identifyPolymers(rxn);
  console.log(polymerIndices);
  let requestData = {
    name: name,
    molWeight: molWeight,
    hPhrase: hPhrase,
    density: density,
    concentration: concentration,
    smiles: smiles,
    workbook: workbook,
    workgroup: workgroup,
    component: "component",
    source: "sketcher",
    // polymerMode: polymerMode,
    // polymer: polymer, // from reaction_table/routes via html
    polymerIndices: polymerIndices,
  };

  // sends polymers to polymer db
  if (polymerIndices.length) {
    novelCompoundURL = "/_polymer_novel_compound";
  } else {
    novelCompoundURL = "/_novel_compound";
    requestData.cas = cas;
  }

  $.ajax({
    url: novelCompoundURL,
    type: "post",
    dataType: "json",
    data: requestData,
  }).done(function (data) {
    alert(data.feedback);
    if (data.feedback === "Compound added to the database") {
      $("#js-novel-compound-input-form").hide();
      $("#js-load-status").val("loaded");
      $("#action-button-submit").click();
    }
    let reactionID = $("#js-reaction-id").val();
    let reactionSmiles = $("#js-reaction-smiles").val();
    updateReactionTable(reactionSmiles, workgroup, workbook, reactionID);
  });
}

function sketcherNovelCompoundInputValidate() {
  let toValidate = $("#js-new-compound-name");
  let submitButtom = $("#js-new-compound-submit");
  if (toValidate.val() == "") {
    submitButtom.attr("disabled", "disabled");
  }
  toValidate.keyup(function () {
    let empty = false;
    toValidate.each(function () {
      if ($(this).val() == "") {
        empty = true;
      }
    });
    if (empty) {
      submitButtom.attr("disabled", "disabled");
    } else {
      submitButtom.removeAttr("disabled");
    }
  });
}

function makeSolventInput() {
  let newSolventNumber = Number($("#js-number-of-solvents").val()) + 1;
  novelCompoundInput("solvent", null, newSolventNumber);
}

function makeReagentInput() {
  let newReagentNumber = Number($("#js-number-of-reagents").val()) + 1;
  novelCompoundInput("reagent", null, newReagentNumber);
}

function novelCompoundInput(component, cas = null, x = null) {
  // shows the novel reagent or solvent input fields
  toggleInputButton(component);
  $(".not-novel-compound-input").css("opacity", 0.4);
  $(`#js-input-${component}1`).show();
  $(`#js-input-${component}`).val(x);
  $(`#js-input-${component}2`).show();
  $(`#js-add-new-${component}-by-table`).hide();
  $(`.js-add-${component}`).hide();
  $(`.js-remove-${component}`).hide();
  $(`#js-input-${component}-cas`).val(cas);
  autoChangeRequiredStyling(`#js-input-${component}-name`);
}

function cancelNovelCompoundInput(component) {
  $(".not-novel-compound-input").css("opacity", 1);
  $(`#js-input-${component}1`).hide();
  $(`#js-input-${component}2`).hide();
  $(`#js-add-new-${component}-by-table`).show();
  $(`.js-add-${component}`).show();
  $(`.js-remove-${component}`).show();
}

//enabling/disabling input button when all necessary fields are filled/not filled
function toggleInputButton(component) {
  $(`#js-input-${component}`).attr("disabled", "disabled");
  let toValidate = $(`#js-input-${component}-name`);
  toValidate.keyup(function () {
    let empty = false;
    toValidate.each(function () {
      if ($(this).val() == "") {
        empty = true;
      }
    });
    if (empty) {
      $(`#js-input-${component}`).attr("disabled", "disabled");
    } else {
      $(`#js-input-${component}`).removeAttr("disabled");
    }
  });
}

function submitNovelCompoundViaTable(component) {
  let newCompoundInputNumber = $(`#js-input-${component}`).val();
  let name = getVal(`#js-input-${component}-name`);
  let smiles = $(`#js-input-${component}-smiles`).val();
  let hazards = $(`#js-input-${component}-hazards`).val();
  let cas = $(`#js-input-${component}-cas`).val();
  let density = $(`#js-input-${component}-density`).val();
  let molWeight = $(`#js-input-${component}-mw`).val();
  let concentration = $(`#js-input-${component}-concentration`).val();
  let workbook = $("#js-active-workbook").val();
  let workgroup = $("#js-active-workgroup").val();
  let confirmVal = true;
  if (molWeight === "" && smiles === "") {
    confirmVal = confirm(
      "If you are sure the compound has no molecular weight, press OK to proceed",
    );
  }
  if (confirmVal === true) {
    $.ajax({
      url: "/_novel_compound",
      type: "post",
      dataType: "json",
      data: {
        name: name,
        smiles: smiles,
        hPhrase: hazards,
        cas: cas,
        density: density,
        molWeight: molWeight,
        concentration: concentration,
        workgroup: workgroup,
        workbook: workbook,
        component: component,
        source: "table",
      },
      success: function (response) {
        alert(response.feedback);
        if (response.feedback === "Compound added to the database") {
          // restore form and show buttons again
          $(`.not-novel-compound-input`).css("opacity", 1);
          $(`#js-input-${component}-name`).val("");
          $(`#js-input-${component}-smiles`).val("");
          $(`#js-input-${component}-hazards`).val("");
          $(`#js-input-${component}-cas`).val("");
          $(`#js-input-${component}-density`).val("");
          $(`#js-input-${component}-mw`).val("");
          $(`#js-input-${component}-concentration`).val("");
          $(`#js-input-${component}1`).hide();
          $(`#js-input-${component}2`).hide();
          $(`#js-add-new-${component}-by-table`).show();
          $(`.js-add-${component}`).show();
          $(`.js-remove-${component}`).show();
          // check if input for novel compound exists
          let $searchID = $(`#js-${component}${newCompoundInputNumber}`);
          // add reagent/solvent if the required search id does not exist already
          if (!$searchID.length) {
            $(`.js-add-${component}`).click();
          }
          // enter name and post data to either reagent or solvent routes
          $searchID.val(name);
          if (component === "reagent") {
            postReagentData(name, newCompoundInputNumber);
          } else if (component === "solvent") {
            let numberOfSolvents = getNum("#js-number-of-solvents");
            for (let i = 1; i < numberOfSolvents + 1; i++) {
              $(`#js-solvent-datalist${i}`).append(
                `<option value="${name}" class="datalistOption hazard-reset-hazard">${name}`,
              );
            }
            $("#js-solvent-datalist").append(
              `<option value="${name}" class="datalistOption hazard-reset-hazard">${name}`,
            );
            postSolventData(name, newCompoundInputNumber);
          }
          let inputName = $(`#js-input-${component}-name`);
          inputName
            .removeClass("remove-highlight-filled-cell")
            .addClass("add-highlight-unfilled-cell");
        }
      },
    });
  }
}
