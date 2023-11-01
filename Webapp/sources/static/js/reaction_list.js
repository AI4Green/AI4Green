/**
 * Gets reactions for specific workbook/workgroup via an AJAX request and adds to the reaction list HTML
 * @param sortCriteria {string} - the sort criteria will be either a-z or time
 */
function getReactions(sortCriteria) {
  let workbook = getVal("#active-workbook");
  if (workbook === "No Workbooks to Show" || workbook === "None") {
    document.getElementById("export-div").style.display = "none";
    document.getElementById("new-reaction").style.display = "none";
  } else {
    document.getElementById("reaction-content").style.display = "block";
    document.getElementById("reaction-column").style.display = "block";
    document.getElementById("no-reactions").style.display = "none";
    let workgroup = getVal("#active-workgroup");
    let workbook = getVal("#active-workbook");
    $.ajax({
      url: "/get_reactions",
      type: "post",
      dataType: "json",
      data: {
        workbook: workbook,
        workgroup: workgroup,
        sortCriteria: sortCriteria,
      },
      success: function (data) {
        $("#reaction-details").html(data.reactionDetails).show(); // sends data to the reaction list
      },
    });
  }
}

/**
 * Calls the getReactions() func if the workbook dropdown is changed, which updates the reaction list.
 */
function updateSelectedWorkbook() {
  getReactions("AZ");
}

/**
 * loads and resets the modal window to make a new reaction
 * and assigns the reaction id corresponding to the active workbook
 */
function newReactionModalWindow() {
  let reactionIDs = getVal("#workbook_corresponding_next_reaction_ids");
  let reactionIDsDic = JSON.parse(reactionIDs);
  let workbook = getVal("#active-workbook");
  let activeReactionID = reactionIDsDic[workbook];
  $("#new-reaction-id").val(activeReactionID);
  $("#new-reaction-name").val("");
  $("#error-warning-new-reaction").html("");
}

/**
 * Validates the input of the modal window.
 * Upon successful validation creates a new reaction and redirects browser to the page for the new reaction
 */
function newReactionCreate() {
  // creates new reaction if name and ID pass validation in the backend routes.
  let workgroup = getVal("#active-workgroup");
  let workbook = getVal("#active-workbook");
  let reactionName = getVal("#new-reaction-name");
  let reactionID = getVal("#new-reaction-id");

  $.ajax({
    url: "/new_reaction",
    type: "post",
    datatype: "json",
    data: {
      reactionName: reactionName,
      reactionID: reactionID,
      workgroup: workgroup,
      workbook: workbook,
    },
    success: function (response) {
      if (response.feedback === "New reaction made") {
        window.location.href = `/sketcher/${workgroup}/${workbook}/${reactionID}/no`;
      } else {
        $("#error-warning-new-reaction").html(response.feedback);
      }
    },
  });
}

/**
 * Calls the delete_reaction routes after the reaction after confirmation from the user
 * @param reaction {HTMLElement}
 */
function deleteReaction(reaction) {
  const completeConfirm = window.confirm(
    "Are you sure you want to permanently delete this reaction?",
  );
  if (completeConfirm === false) {
    return;
  }
  let workbook = getData(reaction, "workbook");
  let workgroup = getData(reaction, "workgroup");
  // id of reaction element is the reaction_id of the reaction
  location.href = `/delete_reaction/${reaction.id}/${workgroup}/${workbook}`;
}

/**
 * Redirects the browser to the page of the reloaded reaction via the main/sketcher routes
 * @param reaction {HTMLElement}
 */
function redirectToReloadReaction(reaction) {
  let workbook = getData(reaction, "workbook");
  let workgroup = getData(reaction, "workgroup");
  // search results open in new tab, otherwise reloaded reactions open in current tab
  if (window.location.pathname.split("/")[1] === "search") {
    window.open(
      `/sketcher/${workgroup}/${workbook}/${reaction.id}/no`,
      "_blank",
    );
  } else {
    window.location.href = `/sketcher/${workgroup}/${workbook}/${reaction.id}/no`;
  }
}

/**
 *
 * @param sortCriteria {string} - the sort criteria, either a-z or time
 * @param workbook {string} - the active workbook name
 * @param workgroup {string} - the active workgroup name
 * @param size {string} - the size of the reaction scheme image
 * @return {Promise<Array>} - array of reaction scheme images formatted as an svg
 */
function getSchemata(sortCriteria, workbook, workgroup, size) {
  return new Promise(function (resolve, reject) {
    // post to get_schemata and get the schemes for reaction images
    $.ajax({
      method: "POST",
      url: "/get_schemata",
      dataType: "json",
      data: {
        sortCriteria: sortCriteria,
        workgroup: workgroup,
        workbook: workbook,
        size: size,
      },
      success: function (data) {
        resolve(data.schemes);
      },
      error: function () {
        reject("error accessing /get_schemata");
      },
    });
  });
}

/**
 * Sorts reactions alphabetically, reordering the cards in the list
 */
function sortReactionsAlphabetically() {
  let reactionCards = $(".reaction-card");
  // sorts by reaction name
  reactionCards.sort((a, b) =>
    a
      .querySelector(".reaction-name")
      .firstChild.innerHTML.toLowerCase()
      .localeCompare(
        b.querySelector(".reaction-name").firstChild.innerHTML.toLowerCase(),
      ),
  );
  $("#reaction-list").empty();
  $("#reaction-list").append(reactionCards);
}

/**
 * Sorts reactions in date order, reordering the cards in the list
 */
function sortReactionsByTime() {
  let reactionCards = $(".reaction-card");
  reactionCards.sort(function (a, b) {
    return (
      new Date(b.querySelector(".reaction-time").innerHTML) -
      new Date(a.querySelector(".reaction-time").innerHTML)
    );
  });
  $("#reaction-list").empty();
  $("#reaction-list").append(reactionCards);
}

/**
 * Calls the export_data_pdf routes to open a page where all reactions schemes are in the PDF for the workbook
 */
function getpdf() {
  let workgroup = getVal("#active-workgroup");
  let workbook = getVal("#active-workbook");
  let sortCriteria = getVal("#js-sort-crit");
  window
    .open(`/export_data_pdf/${workgroup}/${workbook}/${sortCriteria}`, "_blank")
    .focus();
}

/**
 * Calls the export_data_csv routes to prompt download of a csv with all reaction data for the workbook
 */
function getcsv() {
  let workgroup = getVal("#active-workgroup");
  let workbook = getVal("#active-workbook");
  window.location = `/export_data_csv/${workgroup}/${workbook}`;
}

/**
 * Loads the saved reactions using the selected sort criteria including images of reaction schemes
 * @return {Promise<void>}
 */
async function showSavedReactionsSchemes() {
  let sortCriteria = getVal("#js-sort-crit");
  let workgroup = getVal("#active-workgroup");
  let workbook = getVal("#active-workbook");
  let schemes = await getSchemata(sortCriteria, workbook, workgroup, "small");
  for (const [idx, scheme] of schemes.entries()) {
    let idx1 = idx + 1;
    $(`#image${idx1}`).append($("<div>").html(scheme));
  }
}
