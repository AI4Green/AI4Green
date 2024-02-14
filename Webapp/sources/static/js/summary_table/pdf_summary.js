/**
 * Asynchronously generates a PDF summary of the reaction, sends it to the backend, and saves it.
 * @param {string} mode - either none or "locked" which indicates the summary is generated for the locked reaction
 */
async function makePDF(mode) {
  // Create the main element to be included in the PDF
  const elementToPrint = createPDFElement(mode);

  // Generate a blob from the HTML element using html2pdf
  const blob = await html2pdf()
    .set(getPDFOptions())
    .from(elementToPrint)
    .output("blob");

  // Create FormData and append the file and other variables to it
  const formData = createFormData(blob);

  // Send the FormData to the backend
  await sendPDFToBackend(formData);
}

/**
 * Creates the main HTML element to be included in the PDF summary.
 * @param {string} mode - either none or "locked" to indicate reaction is being locked
 * @returns {HTMLDivElement} The created HTML element.
 */
function createPDFElement(mode) {
  // clone variables to add to PDF summary
  const element = document.getElementById("print-container").cloneNode(true);
  const elementTitle = document.getElementById("name-id").cloneNode(true);
  const elementExperimental = createExperimentalElement();
  const dateTime = getDateTime();
  const elementDateTime = createDateTimeElement(dateTime);

  // make new element to append all desired elements for pdf summary
  const elementToPrint = document.createElement("div");
  elementToPrint.id = "outer-print-div";
  elementToPrint.appendChild(elementTitle);
  elementToPrint.appendChild(element);
  elementToPrint.appendChild(elementExperimental);
  elementToPrint.appendChild(elementDateTime);
  // if locking the reaction also include this in the summary
  if (mode === "locked") {
    const elementLockedTime = document.createElement("div");
    elementLockedTime.innerHTML = `<br><br><b>Reaction locked on: ${dateTime}</b>`;
    elementToPrint.appendChild(elementLockedTime);
  }
  return elementToPrint;
}

/**
 * Creates an HTML element containing the content from the experimental textarea.
 * @returns {HTMLDivElement} The created HTML element.
 */
function createExperimentalElement() {
  const textareaElement = document.getElementById("js-reaction-description");
  const elementExperimental = document.createElement("div");
  elementExperimental.innerHTML = $("<div>").text(textareaElement.value).html();
  elementExperimental.classList.add("break-word");
  return elementExperimental;
}

/**
 * Creates an HTML element containing the formatted date and time.
 * @param {string} dateTime - The formatted date and time.
 * @returns {HTMLDivElement} The created HTML element.
 */
function createDateTimeElement(dateTime) {
  // Create a new div element and put content from experimental textarea in to it because textarea does not work in html2pdf
  const elementDateTime = document.createElement("div");
  elementDateTime.innerHTML = `<br><br><b>PDF Summary generated on: ${dateTime}</b>`;
  return elementDateTime;
}

/**
 * Retrieves the PDF generation options.
 * @returns {Object} The PDF generation options.
 */
function getPDFOptions() {
  return {
    margin: 10,
    filename: "myfile.pdf",
    image: { type: "jpeg", quality: 0.98 },
    html2canvas: { scale: 1.2, scrollY: 0 },
    jsPDF: { orientation: "landscape", format: [420, 594] },
  };
}

/**
 * Creates a FormData and appends the file and other variables to it.
 * @param {Blob} blob - The PDF blob to append.
 * @returns {FormData} The created FormData.
 */
function createFormData(blob) {
  const formData = new FormData();
  appendReactionDataToForm(formData);
  formData.append("pdfFile", blob, `${formData.get("reactionID")}-summary.pdf`);
  return formData;
}

/**
 * Sends FormData to the backend.
 * @param {FormData} formData - The FormData to send.
 * @returns {Promise} A promise representing the completion of the backend request.
 */
async function sendPDFToBackend(formData) {
  return fetch("/pdf", {
    method: "POST",
    body: formData,
  }).then(() => {
    updateFileAttachmentList();
  });
}

/**
 * Fetch request to backend to get the updated reaction attachments in order for the current reaction.
 */
function updateFileAttachmentList() {
  const { workgroup, workbook, reactionID } = getReactionFormVariables();
  return fetch("/get_file_attachment_list", {
    method: "POST",
    body: JSON.stringify({
      workgroup: workgroup,
      workbook: workbook,
      reactionID: reactionID,
    }),
    headers: {
      "Content-Type": "application/json",
    },
  })
    .then((response) => response.json()) // Parse the JSON response
    .then((data) => {
      // empty the list of attachments then repopulate with the response data
      $("#file-list").empty();
      showFileAttachmentButtons(data.file_attachments);
    });
}

/**
 * Gets the current datetime for UK timezone and formats.
 * @return {string} - formatted datetime string.
 */
function getDateTime() {
  // Get the current date and time in UTC
  const currentUtcDate = new Date();
  // Specify the target timezone
  const targetTimeZone = "Europe/London";
  // Create an Intl.DateTimeFormat object with the specified timezone
  const dateTimeFormat = new Intl.DateTimeFormat("en-GB", {
    dateStyle: "full",
    timeStyle: "long",
    timeZone: targetTimeZone,
  });
  // Format the current date and time in the UK timezone
  return dateTimeFormat.format(currentUtcDate);
}

/**
 * Displays a loading overlay with a spinning flask whilst saving a PDF summary of the reaction.
 * @param {string} [message='Please do not exit or refresh the page'] - The message to display on the loading overlay.
 */
function showLoadingOverlay(message) {
  let overlay = document.getElementById("loading-overlay");
  let messageElement = document.getElementById("loading-message");
  let iconElement = document.getElementById("loading-icon");

  if (!overlay) {
    // Create overlay and message elements if not exists
    overlay = document.createElement("div");
    overlay.id = "loading-overlay";
    document.body.appendChild(overlay);

    messageElement = document.createElement("div");
    messageElement.id = "loading-message";
    overlay.appendChild(messageElement);

    iconElement = document.createElement("i");
    iconElement.id = "summary-loading-icon";
    iconElement.className = "fa fa-flask";
    overlay.appendChild(iconElement);
  }

  // Set the message
  messageElement.textContent =
    message + "\nPlease do not exit or refresh the page";

  // Show the loading overlay
  overlay.style.display = "block";
}

/**
 * Hides the loading overlay.
 */
function hideLoadingOverlay() {
  let overlay = document.getElementById("loading-overlay");
  if (overlay) {
    // Hide the loading overlay
    overlay.style.display = "none";
  }
}

/**
 * Shows the loading overlay during the creation of the PDF
 * @param {String} mode - pdf generation mode. Either 'summary' or 'locked'
 */
function displayOverlayWhilstMakingPDF(mode) {
  setTimeout(async () => {
    await makePDF(mode)
      .then(() => {
        // PDF creation is complete, hide loading circle
        hideLoadingOverlay();
      })
      .catch((error) => {
        // Handle any errors that occurred during PDF creation
        console.error("Error making PDF:", error);
        hideLoadingOverlay(); // Ensure loading circle is hidden even on error
      });
  }, 5000);
}
