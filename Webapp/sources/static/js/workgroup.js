function getReactions(sort_crit){
    if ($('#WB-select').val() === "No Workbooks to Show"){
        document.getElementById("export-div").style.display = "none";
    }
    else {
        document.getElementById("reaction-content").style.display = "block";
        document.getElementById("reaction-column").style.display = "block";
        document.getElementById("no-reactions").style.display = "none";
        let workbook = $("#WB-select").val()
        let workgroup = $('#workgroup_selected').val()
        $.ajax({
            url: '/get_reactions',
            type: 'post',
            dataType: 'json',
            data: {
                workbook: workbook,
                workgroup: workgroup,
                sort_crit: sort_crit
            },
            success: function (data) {
                $('#reaction-details').html(data.reactionDetails).show(); // Sends data to the reaction table
            }
        })
    }
}

function updateSelectedWorkbook () {
    // update inputs and then get reactions for this workbook
    getReactions("AZ");
}

function newReactionModalWindow() {
    // assign reaction id that corresponds to workbook and clear other existing fields in modal window
    let reactionIDs = $("#workbook_corresponding_next_reaction_ids").val()
    let reactionIDsDic = JSON.parse(reactionIDs)
    let workbook = $("#WB-select").val()
    let activeReactionID = reactionIDsDic[workbook]
    $("#new-reaction-id").val(activeReactionID)
    $("#new-reaction-name").val('')
    $("#error-warning-new-reaction").html('')
}

function newReactionCreate() {
    // creates new reaction if name and ID pass validation in the backend routes.
    let workgroup = $("#workgroup_selected").val()
    let workbook = $("#WB-select").val()
    let reactionName = $("#new-reaction-name").val()
    let reactionID = $("#new-reaction-id").val()
    $.ajax({
        url: "/new_reaction",
        type: 'post',
        datatype: 'json',
        data: {
            reactionName: reactionName,
            reactionID: reactionID,
            workgroup: workgroup,
            workbook: workbook
        },
        success: function (response) {
            if (response.feedback === 'New reaction made') {
                window.location.href = "/sketcher/" + workgroup + "/" + workbook + "/" + reactionID + "/no"
            }
            else {
                $("#error-warning-new-reaction").html(response.feedback)
            }
        }
    })
}
