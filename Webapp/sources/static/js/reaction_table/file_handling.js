function downloadFile(clickedFile) {
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let reactionID = $("#js-reaction-id").val();
  let uuid = clickedFile.value;
  $.ajax({
    url: "/_download_reaction_attachment",
    type: "post",
    data: {
      workgroup: workgroup,
      workbook: workbook,
      reactionID: reactionID,
      uuid: uuid,
    },
    success: function (response) {
      let index = clickedFile.id.split("download-file-button").pop();
      let filename = $(`#view-file-button${index}`).html();
      let url = `data:${response.mimetype};base64,` + response.stream;
      let a = $("<a />", {
        href: url,
        download: filename,
      })
        .appendTo("body")
        .get(0)
        .click();
    },
  });
}

function viewFile(clicked_file) {
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let reactionID = $("#js-reaction-id").val();
  let uuid = clicked_file.value;
  $.ajax({
    url: "/view_reaction_attachment",
    type: "get",
    data: {
      workgroup: workgroup,
      workbook: workbook,
      reactionID: reactionID,
      uuid: uuid,
    },
    success: function (response) {
      // open new tab, set the name, then create a new html element for the file attachment, then go to new tab
      let newWindow = window.open("", "_blank");
      newWindow.document.title = response.name;
      let newElement = makeNewElement(response);
      newWindow.document.body.appendChild(newElement);
      newWindow.focus();
    },
  });
}

function makeNewElement(response) {
  // Makes a pdf or image from the file upload
  let newElement;
  if (response.mimetype.includes("pdf")) {
    newElement = makeNewPDFElement(response.stream);
  } else if (response.mimetype.includes("image")) {
    newElement = makeNewImageElement(response.stream);
  }
  return newElement;
}

function makeNewPDFElement(stream) {
  let newElement = document.createElement("embed");
  newElement.style.width = "100%";
  newElement.style.height = "100%";
  newElement.src = stream;
  return newElement;
}

function makeNewImageElement(stream) {
  let newElement = document.createElement("img");
  newElement.src = stream;
  return newElement;
}

function checkNumberOfUploads(filesToUpload) {
  if (filesToUpload.length === 0) {
    alert(
      "No files selected.\nSelect files to upload by clicking the 'Choose files' button.",
    );
    throw "No files selected";
  } else if (filesToUpload.length >= 10) {
    alert("The maximum number of attachments for a reaction is 10.");
    throw "Max file attachments reached";
  }
}

function validateFiles(filesToUpload) {
  let fileFormData = new FormData();
  for (let [idx, file] of Object.entries(filesToUpload)) {
    if (validateFile(file)) {
      fileFormData.append(`file${idx}`, file, file.name);
    } else {
      // exit function
      throw "invalid file";
    }
  }
  return fileFormData;
}

/**
 * Appends reaction data to a formData dictionary. Reaction data required to identify reaction in backend.
 * @param {formData} reactionFormData - a blank or partly populated form.
 */
function appendReactionDataToForm(reactionFormData) {
  const { workgroup, workbook, reactionID } = getReactionFormVariables();
  reactionFormData.append("workgroup", workgroup);
  reactionFormData.append("workbook", workbook);
  reactionFormData.append("reactionID", reactionID);
}

/**
 * Retrieves form-related variables required to identify reaction in backend.
 * @returns {Object} Form-related variables.
 */
function getReactionFormVariables() {
  return {
    workgroup: $("#js-active-workgroup").val(),
    workbook: $("#js-active-workbook").val(),
    reactionID: getVal("#js-reaction-id"),
  };
}

function uploadFiles() {
  // get files from input, basic front-end validation, post to backend for additional validation and save
  maximumNumberOfFileAttachmentsCheck();
  let filesToUpload = $("#upload-files").prop("files");
  checkNumberOfUploads(filesToUpload);
  // iterate through and validate files
  let fileFormData = validateFiles(filesToUpload);
  appendReactionDataToForm(fileFormData);
  $.ajax({
    url: "/_upload_experimental_data",
    type: "post",
    processData: false,
    contentType: false,
    data: fileFormData,
    success: function (response) {
      alert("files successfully uploaded");
      $("#upload-files").val("");
      showFileAttachmentButtons(response.uploaded_files);
    },
    error: function (response) {
      alert("file upload failed");
    },
  });
}

function maximumNumberOfFileAttachmentsCheck() {
  let numberOfFileAttachments = $("#file-list").children().length;
  if (numberOfFileAttachments >= 10) {
    alert("The maximum number of attachments for a reaction is 10.");
    throw "Max file attachments reached";
  }
}

/**
 * Makes a view, delete, and download multibutton for each file attachment in a reaction from the given list.
 * @param uploadedFiles {Array[Object]} - A list of objects with file name and uuid. One item in list per file attachment
 */
function showFileAttachmentButtons(uploadedFiles) {
  let numberOfFileAttachments = $("#file-list").children().length + 1;
  for (let [idx, file] of uploadedFiles.entries()) {
    let newIdx = numberOfFileAttachments + idx;
    let newFileButtonGroup = $("#blank-file-attachment")
      .html()
      .replace(/-x-/g, newIdx);
    $("#file-list").append(newFileButtonGroup);
    $(`#view-file-button${newIdx}`).val(file.uuid).html(file.name);
    $(`#delete-file-button${newIdx}`).val(file.uuid);
    $(`#download-file-button${newIdx}`).val(file.uuid);
    $(`#file-button-group${newIdx}`).show();
    if (file.autogenerated === true) {
      $(`#delete-file-button${newIdx}`).prop("disabled", true);
    }
  }
}

function validateFile(file) {
  let fileSizeValidation = validateFileSize(file);
  if (fileSizeValidation === "failed") {
    alert("Files must be  below 1 mb");
    return false;
  }
  return true;
}

function validateFileSize(file) {
  // if larger than 1 mb
  if (file.size > 1000000) {
    return "failed";
  }
  return "success";
}

function deleteFile(clickedFile) {
  let index = clickedFile.id.split("delete-file-button").pop();
  if (!confirmFileDeletion(index)) {
    return;
  }
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let reactionID = $("#js-reaction-id").val();
  let uuid = clickedFile.value;
  $.ajax({
    url: "/_delete_reaction_attachment",
    type: "delete",
    data: {
      workgroup: workgroup,
      workbook: workbook,
      reactionID: reactionID,
      uuid: uuid,
    },
    success: function (response) {
      // update files shown
      removeDeletedFileAttachment(index);
      updateRemainingIDs(index);
    },
  });
}

function removeDeletedFileAttachment(index) {
  $(`#file-button-group${index}`).remove();
}

function updateRemainingIDs(index) {
  // when a file is deleted, the ids of all subsequent elements must be updated
  let idsToUpdate = [
    "delete-file-button",
    "view-file-button",
    "download-file-button",
    "file-button-group",
  ];
  let numberOfFileAttachments = $("#file-list").children().length + 1;
  index = Number(index);
  // index + 1 is iteration startpoint and numberOfFileAttachments is endpoint
  for (let i = index + 1; i <= numberOfFileAttachments; i++) {
    let newIndex = i - 1;
    if (i > 15) {
      throw "error number of children exceeded 10";
    }
    for (let elementID of idsToUpdate) {
      // update ID for each element
      let newID = elementID + newIndex;
      $(`#${elementID + i}`).attr("id", newID);
    }
  }
}

function confirmFileDeletion(index) {
  let fileName = $(`#view-file-button${index}`).html();
  let text = `You are about to delete file "${fileName}". Press OK to confirm.`;
  // returns true if user confirms, else returns false
  return confirm(text) === true;
}
