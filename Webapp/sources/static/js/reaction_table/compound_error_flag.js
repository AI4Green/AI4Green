// disable report button until an option is selected
$("#compound-data-report-dropdown").on('change', function() {
    if ($(this).find('option:selected').text() === "-select-"){
        $("#report-compound-data-submit").attr('disabled',true)
    }
    else{
        $("#report-compound-data-submit").attr('disabled',false)
    }
}).trigger('change');

function reportWindow(clickedID){
    // reset modal window and then find and put compound name into report window
    $("#compound-data-report-additional-info").val('')
    $("#modal-report-inner-body").show()
    $("#report-compound-data-submit").show()
    $("#compound-data-report-dropdown").val('-select-');
    // e.g., reactant1, or product3, etc.
    let componentName = clickedID.split("js-report-")[1]
    let compoundNameID = "#js-" + componentName
    let compoundName = $(compoundNameID).val()
    let componentNumber = componentName.slice(-1)
    let componentType = componentName.slice(0, -1)
    let compoundIDID = "#js-" + componentType + '-primary-key' + componentNumber
    let compoundID = $(compoundIDID).val()
    $("#compound-data-report-id").val(compoundID)
    $("#compound-data-report-text").html("Do you wish to report an error in the compound data for " + compoundName + "?")
    $("#compound-data-report-modal-title").html(compoundName)
}

function reportCompoundSubmit(){
    // post data to the backend
    let compoundName = $("#compound-data-report-modal-title").html()
    let errorType = $("#compound-data-report-dropdown").val()
    let additionalInfo = $("#compound-data-report-additional-info").val()
    let compoundID = $("#compound-data-report-id").val()
    $.ajax({
        url: '/compound_data_error_report',
        type: 'post',
        datatype: 'json',
        data: {compoundName: compoundName, compoundID: compoundID, errorType: errorType, additionalInfo: additionalInfo},
        success: function(){
            $("#modal-report-inner-body").hide()
            $("#compound-data-report-text").html("Thank you for your report")
            $("#report-compound-data-submit").hide().attr('disabled', true)
        }
    })
}

