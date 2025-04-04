$(async function () {
  showSketcherLoadingCircle();
  await setupNewKetcherSketcher("100%");

  // sleep used to allow sketchers to load scripts and make js Objects
  await sleep(1000);
  // sketcher changes when user clicks radio button
  await switchActiveEditor();
  $('input[name="sketcher-select"]').click(function () {
    switchActiveEditor();
  });
  hideSketcherLoadingCircle();
  updateSelectedWorkGroup();
  await setupNewMarvinSketcher();
});

async function structureSearch(searchType) {
  let workgroup = $("#active-workgroup").val();
  let workbook = $("#active-workbook").val();
  let smiles = await exportSmilesFromActiveEditor();
  let mol = await exportMOLFromActiveEditor();
  $.ajax({
    url: "/structure_search_handler",
    type: "post",
    datatype: "json",
    data: {
      workgroup: workgroup,
      workbook: workbook,
      searchType: searchType,
      smiles: smiles,
      mol: mol,
    },
    success: function (response) {
      $("#search-results-message").text(response.message);
      if (response.status === "success") {
        showSearchReactions(response);
      } else if (response.status === "fail") {
        $("#search-results-contents").empty();
      } else if (response.status === "error message") {
        $("#search-results-contents").empty();
        alert(response.message);
      }
    },
  });
}

async function showSearchReactions(response) {
  $("#search-results-contents").html(response.search_results).show();
  for (const [idx, imgSource] of response.images.entries()) {
    let idx1 = idx + 1;
    $(`#image${idx1}`).attr("src", imgSource);
  }
  document.getElementById("export-div").style.display = "none";
}

function updateSelectedWorkGroup() {
  // updates the workbook+creator dropdowns after change to selected workgroup
  let workgroup = $("#active-workgroup").val();
  $.ajax({
    url: "/updated_workgroup_dropdown",
    type: "post",
    datatype: "json",
    data: { workgroup: workgroup },
    success: function (response) {
      // update dropdown with workbook options
      let dropdownSelect = $("#active-workbook");
      dropdownSelect.empty();
      for (let workbook of response.workbooks) {
        let option = document.createElement("option");
        option.text = workbook;
        option.value = workbook;
        dropdownSelect.append(option);
      }
    },
  });
}
