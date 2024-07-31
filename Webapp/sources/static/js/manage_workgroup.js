/**
 * Handles changes to a user's status in a workgroup, including removal of the user.
 *
 * @param {string} id - The ID of the HTML element containing the user's email.
 * @param {string} mode - The mode of change to apply, such as "remove".
 * @param {string} current_status - The current status of the user within the workgroup.
 *
 */
function make_change(id, mode, current_status){

    if (mode === "remove"){
        const completeConfirm = window.confirm("Are you sure you want to remove this user entirely from the Workgroup? They will also be removed from any Workbooks they are part of within this Workgroup.");
        if (completeConfirm === false) {
            return;
        }
    }
    let email = getVal($("#" + id))
    let workgroup = getVal($("#current_workgroup"));
    fetch('/manage_workgroup/make_change/' + workgroup + '/' + email + '/' + mode + '/' + current_status).then(function(response) {
        response.json().then(function(data) {
            alert(data.feedback);
            window.location.href = "/manage_workgroup/" + workgroup;
        });
    });
}

/**
 * Handles a request to change a user's status in a workgroup.
 *
 * @param {string} id - The ID of the HTML element containing the user's email.
 * @param {string} mode - The mode of change to apply, such as "approve" or "reject".
 * @param {string} decision - The decision to be made regarding the user's status.
 *
 */
function make_change_request(id, mode, decision) {
    let email = getVal($("#" + id))
    let workgroup = getval($("#current_workgroup"));
    fetch('/manage_workgroup/change_status_request/' + workgroup + '/' + email + '/' + mode + '/' + decision).then(function(response) {
        response.json().then(function(data) {
            alert(data.feedback);
            window.location.href = "/manage_workgroup/" + workgroup + "/yes";
        });
    });
}

/**
 * Displays the modal for adding a user by email.
 */
function addUserByEmailModal(){
    modal = $("#add-user-modal").modal("show")
}

/**
 * Generates a QR code for the current workgroup and displays it in a modal.
 *
 */
function generateQRCode(){
    let currentWorkgroup = getVal($("#current_workgroup"));

    fetch("/generate_qr_code/" + currentWorkgroup, {
        headers: {
        "Content-Type": "application/json",
      },
      method: "POST",
      body: JSON.stringify({
        workgroup: currentWorkgroup
        })
    })
    .then(function(response) { return response.json() })
    .then(function(item) {
        document.getElementById("qr-code").setAttribute("src", "data:image/png;base64," + item)
        $("#qr-code-modal").modal("show")
    })
}

/**
 * Initiates the print dialog to print the QR code.
 */
function printQRCode(){
    window.print()
}