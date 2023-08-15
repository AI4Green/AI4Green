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