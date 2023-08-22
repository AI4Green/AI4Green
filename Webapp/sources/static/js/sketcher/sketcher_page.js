$(function() {
    if ($("#js-demo").val() === 'not demo') {
        lockReactionButtonHandler()
        reactionNoteButtonHandler()
        hideAddendaIfNoNotes()
    }
});

function lockReactionButtonHandler() {
    // disable complete reaction button once reaction is completed
    if ($("#js-complete").val() === "not complete") {
        $("#mark-complete").attr("disabled", false)
    }
}

function reactionNoteButtonHandler(){
    if ($("#js-complete").val() === "complete") {
        $("#reaction-note-button").show();
    }
    else {
        $("#reaction-note-button").hide();
    }
}

function hideAddendaIfNoNotes(){
    // hide the addenda background if there are no notes
    if ($("#addenda-cardholder").children().length == 0){
        $("#addenda-content").children().hide()
    }
}
