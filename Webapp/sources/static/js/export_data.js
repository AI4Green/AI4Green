$(function () {
  updateSelectedWorkGroup("export_data");
});

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
