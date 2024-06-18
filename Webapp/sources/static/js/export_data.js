$(async function () {
  let currentURL = window.location.href;
  // These functions run on load of export_data/home
  if (currentURL.includes("export_data/home")) {
    // set up the page
    await updateSelectedWorkGroup("export_data");
    checkExportPermissions();
    initiateReactionLists();
    populateReactionIDList();
    $('[data-toggle="tooltip"]').tooltip();
    // need to update workbooks when the workgroup is changed
    $("#active-workgroup").change(async function () {
      await updateSelectedWorkGroup("export_data");
      $("#active-workbook").trigger("change");
    });
    // need to update reaction list when the workbook is updated
    $("#active-workbook").change(function () {
      if (checkValidWorkbookAndWorkgroup()) {
        checkExportPermissions();
        populateReactionIDList();
      }
    });
  }
});

/**
 * Checks neither workgroup or workbook are null and hence invalid
 * @return {boolean} true if the workbook/group combo is valid, false if not
 */
function checkValidWorkbookAndWorkgroup() {
  if (
    $("#active-workbook").val() === null ||
    $("#active-workgroup").val() === null
  ) {
    $("#included-list").empty();
    $("#excluded-list").empty();
    $("#export-permission-result").html(
      'Please enter a valid workgroup and workbook <i class="bi bi-x" style="color:red"></i>',
    );
    return false;
  }
  return true;
}

/**
 *  Enables user to move items between the lists in the reaction selection modal window.
 */
function initiateReactionLists() {
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
  let workgroup = $("#active-workgroup").val();
  let workbook = $("#active-workbook").val();
  fetch("/export_permission", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      workgroup: workgroup,
      workbook: workbook,
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
          'Only a principal investigator of a workgroup has permission to export data from a workbook. <i class="bi bi-x" style="color:red"></i>',
        );
      } else {
        $("#export-permission-result").html();
      }
    });
}

/**
 * Gets the list of valid reactions (have a recorded yield) to populate the list in the modal window
 * @return {Promise<any>}
 */
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

/**
 * Empties the current lists in the modal window and populates with reactions from the newly selected workbook
 */
function populateReactionIDList() {
  getReactionIDList().then((data) => {
    // Populate the included column with all reaction IDs initially
    $("#excluded-list").empty(); // Clear existing content
    const $includedList = $("#included-list");
    $includedList.empty(); // Clear existing content
    data.forEach((item) => {
      const $listItem = $("<a>")
        .addClass("list-group-item list-group-item-action")
        .text(item);
      $includedList.append($listItem);
    });
    sortList("#included-list");
  });
}

/**
 * Opens the modal window to see the reaction list
 */
function openReactionModal() {
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
        `${exportFormat} Data export request created. You will receive an email notification regarding your export
        once all of the workgroup principal investigators have approved the request.`,
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
      $requestFeedback.text(
        "Request failed. This could happen if there are no reactions selected for export. " +
          "For further support get in touch through our help page.",
      );
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
