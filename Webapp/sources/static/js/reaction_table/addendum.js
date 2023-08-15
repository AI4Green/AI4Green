function saveNewReactionNoteToDatabase(){
    let reactionNoteText = $("#new-reaction-note-text").val()
    let workgroup = $("#js-active-workgroup").val()
    let workbook = $("#js-active-workbook").val()
    let reactionID = $("#js-reaction-id").val()
    $.ajax({
        url: '/_save_reaction_note',
        type: 'post',
        data: {
            reactionNoteText: reactionNoteText,
            workgroup: workgroup,
            workbook: workbook,
            reactionID: reactionID
        },
        dataType: 'json',
        success: function(response) {
            showReactionNote(response.reaction_note)
        }
    });
}


function showReactionNote(reactionNote){
    // get number of existing cards and add one to get id for new note
    let newNoteNumber = $("#addenda-cardholder").children().length + 1
    // use template, prepend it to reaction note section and add in details
    let newNote = $("#blank-addendum").html().replace(/-x-/g, newNoteNumber)
    $("#addenda-cardholder").prepend(newNote)
    $(`#addendum-text${newNoteNumber}`).text(reactionNote.text)
    $(`#addendum-author${newNoteNumber}`).text(reactionNote.author.user.fullname)
    $(`#addendum-time${newNoteNumber}`).text(reactionNote.time_of_creation)
    $(`#addendum-card${newNoteNumber}`).css('display', 'block');
    $("#addenda-content").children().show()
}
