$(function () {
  updateSelectedWorkGroup("export_data");
  checkExportPermissions();

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
    body: JSON.stringify({
      workgroup: $("#active-workgroup").val(),
      workbook: $("#active-workbook").val(),
    }),
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
      handleNewExportRequestResponse(data.feedback);
    })
    .catch((error) => {
      console.error("Error:", error);
    });
}

/**
 * Handles the response
 * @param feedback
 */
function handleNewExportRequestResponse(feedback) {
  switch (feedback) {
    case "Data export request created. Your export will be accessible once all of the workgroup principal investigators have approved the request.":
      console.log("Show success");
      break;
    case "You do not have permission to export this data. You must be either the principal investigator or workbook member":
      console.log("Show failed");
      break;
    default:
      console.log("Show problem");
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