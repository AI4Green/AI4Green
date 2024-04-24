$(function () {
  let currentURL = window.location.href;
  // These functions run on load of export_data/home
  if (currentURL.includes("export_data/home")) {
    updateSelectedWorkGroup("export_data");
    checkExportPermissions();
    initiateReactionLists();
    $('[data-toggle="tooltip"]').tooltip();
  }
});

function initiateReactionLists() {
  // These functions enable user to move items between the lists in the reaction selection modal window.
  $(document).on("click", "#excluded-list .list-group-item", function () {
    const $item = $(this);
    $item.detach().appendTo("#included-list");
    sortList("#included-list");
  });

  $(document).on("click", "#included-list .list-group-item", function () {
    const $item = $(this);
    $item.detach().appendTo("#excluded-list");
    sortList("#excluded-list");
  });
}

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
    body: JSON.stringify({
      workgroup: $("#active-workgroup").val(),
      workbook: $("#active-workbook").val(),
    }),
  })
    .then((response) => response.json())
    .then((permissionResult) => {
      if (permissionResult === "permission accepted") {
        $("#export-permission-result").html(
          'You have permission to export data from this workbook <i class="bi bi-check" style="color:green"></i>',
        );
      } else if (permissionResult === "permission denied") {
        $("#export-permission-result").html(
          'Only a principal investigator of a workgroup has permission to export data from a workbook. <i class="bi bi-cross" style="color:red"></i>',
        );
      } else {
        $("#export-permission-result").html();
      }
    });
}

function getReactionIDList() {
  return fetch("/get_reaction_id_list", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      workgroup: $("#active-workgroup").val(),
      workbook: $("#active-workbook").val(),
    }),
  })
    .then((response) => {
      if (!response.ok) {
        throw new Error("Failed to fetch reaction ID list");
      }
      return response.json();
    })
    .catch((error) => {
      console.error("Error:", error);
      throw error;
    });
}

function openReactionModal() {
  // if anything in left column don't fetch.
  const $excludedList = $("#excluded-list");
  if ($excludedList.children().length === 0) {
    getReactionIDList().then((data) => {
      // Populate the included column with all reaction IDs initially
      const $includedList = $("#included-list");
      $includedList.empty(); // Clear existing content
      data.forEach((item) => {
        const $listItem = $("<a>")
          .addClass("list-group-item list-group-item-action")
          .text(item);
        $includedList.append($listItem);
      });
    });
  }
  $("#reactionsModal").modal("show");
}

/**
 * Sends data to the backend to make a new data export request
 * @param exportFormat {string} - the format
 */
function makeDataExportRequest(exportFormat) {
  let workgroup = $("#active-workgroup").val();
  let workbook = $("#active-workbook").val();
  // get reaction_id from the included list. we use map to get the text of each element which is the reaction id.
  let reactionIDList = $("#included-list")
    .children(".list-group-item")
    .map(function () {
      return $(this).text(); //
    })
    .get();

  fetch("/export_data/new_request", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      workgroup: workgroup,
      workbook: workbook,
      reactionIDList: reactionIDList,
      exportFormat: exportFormat,
    }),
  })
    .then((response) => {
      if (!response.ok) {
        throw new Error("Response was not ok");
      }
      return response.json();
    })
    .then((data) => {
      handleNewExportRequestResponse(data, exportFormat);
    })
    .catch((error) => {
      console.error("Error:", error);
    });
}

/**
 * Handles the response
 * @param feedback {string} - the result from the backend when creating the export request
 * @param exportFormat {string} - the data format of the requested export
 */
function handleNewExportRequestResponse(feedback, exportFormat) {
  let $requestFeedback = $("#request-feedback");
  switch (feedback) {
    case "permission accepted":
      $requestFeedback.text(
        `${exportFormat} Data export request created. Your export will be accessible once all of the workgroup principal investigators have approved the request.`,
      );
      $(`#${exportFormat.toLowerCase()}-button`)
        .prop("disabled", true)
        .tooltip("hide");
      break;
    case "permission denied":
      $requestFeedback.text(
        "You do not have permission to export this data. You must be either the principal investigator or workbook member",
      );
      break;
    default:
      $requestFeedback.text("Request Failed Unexpectedly");
      break;
  }
}
/**
 *
 * @param listId {string} - jquery id selector. Expected to be #list-id
 */
function sortList(listId) {
  const $list = $(listId);
  const items = $list.children(".list-group-item").get();
  items.sort(function (reactionA, reactionB) {
    // get the numbers at the end of the reaction id and sort by these
    const numberA = $(reactionA).text().split("-").pop();
    const numberB = $(reactionB).text().split("-").pop();
    return parseInt(numberA) - parseInt(numberB);
  });
  $list.empty().append(items);
}

/**
 *
 * @param exportID {str} the data export request ID
 */
function dataExportRequestDeny(exportID) {
  fetch("/export_data/export_denied", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      exportID: exportID,
    }),
  }).then((response) => {
    if (response.redirected) {
      const redirectUrl = new URL(response.url);
      redirectUrl.searchParams.set("message", "Data export rejected");
      window.location.href = redirectUrl.toString();
    }
  });
}

/**
 *
 * @param exportID the data export request ID
 */
function dataExportRequestApprove(exportID) {
  fetch("/export_data/export_approved", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      exportID: exportID,
    }),
  }).then((response) => {
    if (response.redirected) {
      const redirectUrl = new URL(response.url);
      redirectUrl.searchParams.set("message", "Data export approved");
      window.location.href = redirectUrl.toString();
    }
  });
}

// function export_data() {
//   let workgroup = $("#active-workgroup").val();
//   let workbook = $("#active-workbook").val();
//
//   $.ajax({
//     url: "/export_data_eln_file",
//     type: "post",
//     datatype: "json",
//     data: {
//       workgroup: workgroup,
//       workbook: workbook,
//       // exportType: exportType,
//     },
//     success: function (response) {
//       if (response.status === "approved") {
//         showSuccessMessage(response);
//       } else if (response.status === "not approved") {
//         showFailureMessage(response);
//       }
//     },
//   });
// }

function showSuccessMessage() {
  $("#export-response").html(response.success_message).show();
}

function showFailureMessage() {}
