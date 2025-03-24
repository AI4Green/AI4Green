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
  $("#new-reaction-id-input").val(activeReactionID);
  $("#new-reaction-name").val("");
  $("#error-warning-new-reaction").html("");
}

function cloneReactionModalWindow(reaction) {
  let workbook = getData(reaction, "workbook");
  let workgroup = getData(reaction, "workgroup");
  let oldReactionID = reaction.id;
  let name = "Repeat of Reaction " + oldReactionID;

  fetch("/get_new_reaction_id", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
    body: JSON.stringify({
      workgroup: workgroup,
      workbook: workbook,
    }),
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (newReactionID) {
      if (newReactionID !== "Bad Request") {
        $("#new-reaction-id-input").val(newReactionID);
        $("#old-reaction-id").val(oldReactionID);
        $("#new-reaction-name").val(name);
        $("#new-reaction-data-submit").attr("onclick", "cloneReaction()");
        $("#new-reaction-modal").modal("show");
      } else {
        window.alert("Something went wrong. Please try again.");
      }
    });
}

function cloneReaction() {
  let reactionWorkgroup = getVal("#workgroup");
  let reactionWorkbook = getVal("#workbook");
  let reactionName = getVal("#new-reaction-name");
  let oldReactionID = getVal("#old-reaction-id");
  let newReactionID = getVal("#new-reaction-id-input");
  $.ajax({
    url: "/clone_reaction",
    type: "post",
    datatype: "json",
    data: {
      reactionName: reactionName,
      workbook: reactionWorkbook,
      workgroup: reactionWorkgroup,
      reactionID: oldReactionID,
      newReactionID: newReactionID,
    },
    success: function (response) {
      if (response.feedback === "New reaction made") {
        window.location.href = `/sketcher/${reactionWorkgroup}/${reactionWorkbook}/${newReactionID}/no`;
      } else {
        $("#error-warning-new-reaction").html(response.feedback);
      }
    },
  });
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
  let reactionID = getVal("#new-reaction-id-input");

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
 * Get images of reaction scheme from database
 * @param sortCriteria {string} - the sort criteria, either a-z or time
 * @param workbook {string} - the active workbook name
 * @param workgroup {string} - the active workgroup name
 * @return {Promise<Array>} - array of reaction images formatted as base64 (text)
 */
function getReactionImages(sortCriteria, workbook, workgroup) {
  return new Promise(function (resolve, reject) {
    // post to get_reaction_images and get the reaction images
    $.ajax({
      method: "POST",
      url: "/get_reaction_images",
      dataType: "json",
      data: {
        sortCriteria: sortCriteria,
        workgroup: workgroup,
        workbook: workbook,
      },
      success: function (data) {
        resolve(data.reaction_images);
      },
      error: function () {
        reject("error accessing /get_reaction_images");
      },
    });
  });
}

/**
 * Get reaction SMILES string from database
 * @param sortCriteria {string} - the sort criteria, either a-z or time
 * @param workbook {string} - the active workbook name
 * @param workgroup {string} - the active workgroup name
 * @return {Promise<Array>} - array of reaction smiles strings
 */
function getSmiles(sortCriteria, workbook, workgroup) {
  return new Promise(function (resolve, reject) {
    // post to get_smiles and get the rxn smiles
    $.ajax({
      method: "POST",
      url: "/get_smiles",
      dataType: "json",
      data: {
        sortCriteria: sortCriteria,
        workgroup: workgroup,
        workbook: workbook,
      },
      success: function (data) {
        resolve(data.smiles);
      },
      error: function () {
        reject("error accessing /get_smiles");
      },
    });
  });
}

/**
 * Get reaction RXN string from database
 * @param sortCriteria {string} - the sort criteria, either a-z or time
 * @param workbook {string} - the active workbook name
 * @param workgroup {string} - the active workgroup name
 * @return {Promise<Array>} - array of reaction RXN strings
 */
function getRXNs(sortCriteria, workbook, workgroup) {
  return new Promise(function (resolve, reject) {
    // post to get_smiles and get the rxn smiles
    $.ajax({
      method: "POST",
      url: "/get_rxns",
      dataType: "json",
      data: {
        sortCriteria: sortCriteria,
        workgroup: workgroup,
        workbook: workbook,
      },
      success: function (data) {
        resolve(data.rxns);
      },
      error: function () {
        reject("error accessing /get_rxns");
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
 * Parse dates in "Created: " format
 */
function parseDate(dateString) {
  // Remove "Created: " prefix and trim whitespace
  const cleanDateString = dateString.replace("Created: ", "").trim();
  return new Date(cleanDateString); // Parse as a date
}

/**
 * Sorts reactions in date order, reordering the cards in the list
 */
function sortReactionsByTime() {
  let reactionCards = $(".reaction-card");
  reactionCards.sort(function (a, b) {
    return (
      parseDate(b.querySelector(".reaction-time").innerHTML) -
      parseDate(a.querySelector(".reaction-time").innerHTML)
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
  window.location = `/export_data/csv/${workgroup}/${workbook}`;
}

/**
 * Save images - updates reaction dict
 */
function saveNewImages(sortCriteria, workbook, workgroup, images) {
  $.ajax({
    url: "/_save_new_images",
    type: "post",
    data: {
      sortCriteria: sortCriteria,
      workgroup: workgroup,
      workbook: workbook,
      images: JSON.stringify(images),
    },
    success: function (response) {},
    error: function (error) {
      // Handle error
      console.error(error);
    },
  });
}

/**
 * Generate missing images with hidden ketcher element
 * @param sortCriteria {string} - the sort criteria, either a-z or time
 * @param workbook {string} - the active workbook name
 * @param workgroup {string} - the active workgroup name
 * @param images {Array} - incomplete list of image strings, where missing are ""
 * @return {Promise<Array>} - Returns list of images
 */
async function regenerateImages(sortCriteria, workbook, workgroup, images) {
  let ketcherFrame = document.getElementById("ketcher-editor");
  await new Promise((resolve) => {
    ketcherFrame.onload = () => resolve();
  });

  let rxns = await getRXNs(sortCriteria, workbook, workgroup);
  console.log(rxns);
  let smiles = await getSmiles(sortCriteria, workbook, workgroup);
  console.log(smiles);

  //indices of empty images
  const missingImagesIdx = images
    .map((entry, index) => (entry === "" ? index : -1))
    .filter((index) => index !== -1);
  let missingImages = [];

  for (const [idx, smile] of smiles.entries()) {
    if (!missingImagesIdx.includes(idx)) {
      // only generate missing imgs
      continue;
    }
    if (smile === "") {
      missingImages.push("Empty Reaction");
      continue;
    }
    try {
      try {
        let rxn = rxns[idx];
        let source = await ketcherFrame.contentWindow.ketcher.generateImage(
          rxn,
        );
        const shrunkenBlob = await shrinkBlobImage(source, 600, 400, 50);
        let imgSource = await convertBlobToBase64(shrunkenBlob);
        missingImages.push(imgSource);
      } catch (error) {
        // catch reactions with no rxn file, use smiles instead.
        let source = await ketcherFrame.contentWindow.ketcher.generateImage(
          smile,
        );
        const shrunkenBlob = await shrinkBlobImage(source, 600, 400, 50);
        let imgSource = await convertBlobToBase64(shrunkenBlob);
        missingImages.push(imgSource);
      }
    } catch (error) {
      missingImages.push("Error generating image");
    }
  }
  // Replace "" with the corresponding image
  for (let i = 0; i < images.length; i++) {
    if (images[i] === "") {
      images[i] = missingImages[0];
      //remove missingImages[0] from list
      missingImages.shift();
    }
  }
  saveNewImages(sortCriteria, workbook, workgroup, images);
  return images;
}

/**
 * Loads the saved reactions using the selected sort criteria, then gets images from database.
 * If no image in database, generate with ketcher.
 * @return {Promise<void>}
 */
async function showSavedReactionsImages() {
  let sortCriteria = getVal("#js-sort-crit");
  let workgroup = getVal("#active-workgroup");
  let workbook = getVal("#active-workbook");

  let images = await getReactionImages(sortCriteria, workbook, workgroup);
  if (images.some((entry) => entry === "")) {
    // if any of the images arent in db, generate now
    images = await regenerateImages(sortCriteria, workbook, workgroup, images);
  }

  for (const [idx, imgSource] of images.entries()) {
    let idx1 = idx + 1;
    let $image = $(`#image${idx1}`);

    if (imgSource === "Empty Reaction") {
      // change img div to text div
      $image.replaceWith(
        '<div id="image${idx1}" style="text-align: center; padding: 20px;">[Empty Reaction]</div>',
      );
    } else if (imgSource === "Error generating image") {
      // change img div to text div
      $image.replaceWith(
        '<div id="image${idx1}" style="text-align: center; padding: 20px;">[Error generating image]</div>',
      );
    }

    $image.attr("src", imgSource);
  }
}
