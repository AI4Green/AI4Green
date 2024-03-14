$(function () {
  updateSelectedWorkGroup("export_data");
  checkExportPermissions();
});

/**
 * Checks export permissions for the current workbook.
 * Sends a POST request to "/export_permission" endpoint with the selected workgroup information.
 * Updates the export permission result message accordingly.
 */
function checkExportPermissions() {
  fetch("/export_permission", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({ workgroup: $("#active-workgroup").val() }),
  })
    .then((response) => response.json())
    .then((permissionResult) => {
      if (
        permissionResult["permission"] ===
        "user has export permission for this workbook"
      ) {
        $("#export-permission-result").html(
          'You have permission to export data from this workbook <i class="bi bi-check" style="color:green"></i>',
        );
      } else if (
        permissionResult["permission"] ===
        "user does not have export permission for this workbook"
      ) {
        $("#export-permission-result").html(
          'Only a principal investigator of a workgroup has permission to export data from a workbook. <i class="bi bi-cross" style="color:red"></i>',
        );
      } else {
        $("#export-permission-result").html();
      }
    });
}

function export_data() {
  let workgroup = $("#active-workgroup").val();
  let workbook = $("#active-workbook").val();

  $.ajax({
    url: "/export_data_eln_file",
    type: "post",
    datatype: "json",
    data: {
      workgroup: workgroup,
      workbook: workbook,
      // exportType: exportType,
    },
    success: function (response) {
      if (response.status === "approved") {
        showSuccessMessage(response);
      } else if (response.status === "not approved") {
        showFailureMessage(response);
      }
    },
  });
}

function showSuccessMessage() {
  $("#export-response").html(response.success_message).show();
}

function showFailureMessage() {}
