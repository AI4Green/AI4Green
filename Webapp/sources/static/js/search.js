$(async function() {
    showSketcherLoadingCircle()
    await setupNewKetcherSketcher("100%");
    await setupNewMarvinSketcher();
    // sleep used to allow sketchers to load scripts and make js Objects
    await sleep(1000)
    // sketcher changes when user clicks radio button
    await switchActiveEditor()
    // setTimeout(switchActiveEditor, 500);
    $('input[name="sketcher-select"]').click(function () {
        switchActiveEditor();
    });
    hideSketcherLoadingCircle()
    updateSelectedWorkGroup()
});

async function structureSearch(){
    let workgroup = $("#active-workgroup").val()
    let workbook = $("#active-workbook").val()
    let searchType = "exact_structure"
    let smiles = await exportSmilesFromActiveEditor()
    $.ajax({
        url: '/structure_search_handler',
        type: 'post',
        datatype: 'json',
        data: {
            workgroup: workgroup,
            workbook: workbook,
            searchType: searchType,
            smiles: smiles
        },
        success: function(response){
            $("#search-results-message").text(response.message)

            if (response.status === 'success'){
                showSearchReactions(response)
            } else if (response.status === 'failure'){
                alert("testing if this still appears")
                // show failure message
            }
        }
    })
}

async function showSearchReactions(response){
    $('#search-results-contents').html(response.search_results).show();
    for (const [idx, scheme] of response.schemes.entries()){
        let idx1 = idx + 1
        $(`#image${idx1}`).append($('<div>').html(scheme))
    }
    document.getElementById("export-div").style.display = "none";
}

function updateSelectedWorkGroup() {
    // updates the workbook+creator dropdowns after change to selected workgroup
    let workgroup = $('#active-workgroup').val()
    $.ajax({
        url: '/updated_workgroup_dropdown',
        type: 'post',
        datatype: 'json',
        data: {workgroup: workgroup},
        success: function(response) {
            // update dropdown with workbook options
            let dropdownSelect = $("#active-workbook")
            dropdownSelect.empty()
            for (let workbook of response.workbooks){
                let option = document.createElement('option');
                option.text = workbook
                option.value = workbook
                dropdownSelect.append(option)
            }
        }
    });
}


