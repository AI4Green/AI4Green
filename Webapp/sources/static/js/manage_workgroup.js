function makeChange(id, mode, current_status){

    if (mode === "remove"){
        const completeConfirm = window.confirm("Are you sure you want to remove this user entirely from the Workgroup? They will also be removed from any Workbooks they are part of within this Workgroup.");
        if (completeConfirm === false) {
            return;
        }
    }
    let email = $("#" + id).val()
    let workgroup = $("#current_workgroup").val();
    fetch('/manage_workgroup/make_change/' + workgroup + '/' + email + '/' + mode + '/' + current_status).then(function(response) {
        response.json().then(function(data) {
            alert(data.feedback);
            window.location.href = "/manage_workgroup/" + workgroup;
        });
    });
}

function make_change_request(id, mode, decision){
    let email = $("#" + id).val()
    let workgroup = $("#current_workgroup").val();
    fetch('/manage_workgroup/change_status_request/' + workgroup + '/' + email + '/' + mode + '/' + decision).then(function(response) {
        response.json().then(function(data) {
            alert(data.feedback);
            window.location.href = "/manage_workgroup/" + workgroup + "/yes";
        });
    });
}

function addUserByEmailModal(){
    modal = $("#add-user-modal").modal("show")
}

