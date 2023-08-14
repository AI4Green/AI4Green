$(document).ready(function() {
    initialiseSearchPage()
});
let searchMarvin

function initialiseSearchPage(){
    initialiseSketcher()
    updateSelectedWorkGroup()
}

async function initialiseSketcher(){
    searchMarvin = await newMarvinSketcher("#marvin-search")
    searchMarvin.setDisplaySettings({"toolbars": "education"});
}

async function structureSearch(){
    let workgroup = $("#workgroup-select").val()
    let workbook = $("#workbook-select").val()
    let searchType = "exact_structure"
    let smiles = await exportSmilesFromSketcher(searchMarvin)
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
    document.getElementById("export-div").style.display = "block";
}

function updateSelectedWorkGroup() {
    // updates the workbook+creator dropdowns after change to selected workgroup
    let workgroup = $('#workgroup-select').val()
    $.ajax({
        url: '/updated_workgroup_dropdown',
        type: 'post',
        datatype: 'json',
        data: {workgroup: workgroup},
        success: function(response) {
            // update dropdown with workbook options
            let dropdownSelect = $("#workbook-select")
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


