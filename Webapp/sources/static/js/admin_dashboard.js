function denyWorkgroup(loopDenyClicked) {
    // function sends an alert to confirm the admin wishes to deny workgroup
    let loopNumber = loopDenyClicked.slice(-1)
    let workgroupId = "#workgroup-name" + loopNumber
    let workgroup = $(workgroupId).html().trim()
    let institution = $("#institution-name" + loopNumber).html().trim()
    let text = "Press OK if you are sure you mean to delete the workgroup:" + workgroup;
        if (confirm(text) == true) {
        window.location.href = 'admin_dashboard/' + institution + '/' + workgroup + '/deny'
      } else {
      }
}

function deleteUser(userInfo) {
    let info = userInfo.split("-")
    // function sends an alert to confirm the admin wishes to delete a user
    let text = "Are you sure you want to delete this user:\n" + "Full Name: " + info[2] + "\nEmail: " + info[1] + "\nUsername: " + info[3] + "\nRole: " + info[4];
    if (confirm(text) == true) {
        if (confirm("Are you sure? This cannot be undone.") == true) {
            window.location.href = '/admin_delete_user/' + info[0] + "#userFocus"
        }
    }
}

function userFocus() {
    if(window.location.hash) {
        $('.nav-tabs a[href="#users_tab"]').tab('show')
    }
}