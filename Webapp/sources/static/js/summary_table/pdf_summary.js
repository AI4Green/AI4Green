// /**
//  * Sends a blob of a PDF summary of the reaction to the backend to be saved
//  */
// async function makePDF() {
//   await sleep(7000);
//
//   let element = document.getElementById("print-container").cloneNode(true);
//   let elementTitle = document.getElementById("name-id").cloneNode(true);
//
//
//   let textareaElement = document.getElementById("js-reaction-description");
//   let elementExperimental = document.createElement("div");
//   elementExperimental.innerHTML = textareaElement.value;
//   elementExperimental.classList.add("break-word");
//   let dateTime = getDateTime();
//   let elementDateTime = document.createElement("div");
//   elementDateTime.innerHTML = `<br><br><b>PDF Summary generated on: ${dateTime}</b>`;
//
//
//   let elementToPrint = document.createElement("div");
//   elementToPrint.id = "outer-print-div";
//   elementToPrint.appendChild(elementTitle);
//   elementToPrint.appendChild(element);
//   elementToPrint.appendChild(elementExperimental);
//   elementToPrint.appendChild(elementDateTime);
//   elementToPrint.classList.add("pdf-print");
//
//   // options for html2pdf
//   const opt = {
//     margin: 10,
//     filename: "myfile.pdf",
//     image: { type: "jpeg", quality: 0.98 },
//     html2canvas: { scale: 1.2, scrollY: 0 },
//     jsPDF: {
//       orientation: "landscape",
//       format: [420, 594],
//     },
//   };
//
//
//   let blob = await html2pdf().set(opt).from(elementToPrint).output("blob");
//
//   // get other variables to send to backend
//   let workgroup = $("#js-active-workgroup").val();
//   let workbook = $("#js-active-workbook").val();
//   let reactionID = getVal("#js-reaction-id");
//
//   // make formData and append the file and others variables to it.
//   const formData = new FormData();
//   formData.append("pdfFile", blob, `${reactionID}-summary.pdf`);
//   formData.append("workgroup", workgroup);
//   formData.append("workbook", workbook);
//   formData.append("reactionID", reactionID);
//   await fetch("/pdf", {
//     method: "POST",
//     body: formData,
//   });
// }

/**
 * Asynchronously generates a PDF summary of the reaction, sends it to the backend, and saves it.
 */
async function makePDF() {
  await sleep(7000);

  // Create the main element to be included in the PDF
  const elementToPrint = createPDFElement();

  // Generate a blob from the HTML element using html2pdf
  const blob = await html2pdf()
    .set(getPDFOptions())
    .from(elementToPrint)
    .output("blob");

  // Get form-related variables
  const { workgroup, workbook, reactionID } = getFormVariables();

  // Create FormData and append the file and other variables to it
  const formData = createFormData(blob, reactionID);

  // Send the FormData to the backend
  await sendFormDataToBackend(formData);
}

/**
 * Creates the main HTML element to be included in the PDF summary.
 * @returns {HTMLDivElement} The created HTML element.
 */
function createPDFElement() {
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
  elementToPrint.classList.add("pdf-print");

  return elementToPrint;
}

/**
 * Creates an HTML element containing the content from the experimental textarea.
 * @returns {HTMLDivElement} The created HTML element.
 */
function createExperimentalElement() {
  const textareaElement = document.getElementById("js-reaction-description");
  const elementExperimental = document.createElement("div");
  elementExperimental.innerHTML = textareaElement.value;
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
 * Retrieves form-related variables.
 * @returns {Object} Form-related variables.
 */
function getFormVariables() {
  return {
    workgroup: $("#js-active-workgroup").val(),
    workbook: $("#js-active-workbook").val(),
    reactionID: getVal("#js-reaction-id"),
  };
}

/**
 * Creates a FormData and appends the file and other variables to it.
 * @param {Blob} blob - The PDF blob to append.
 * @param {string} reactionID - The reaction ID.
 * @returns {FormData} The created FormData.
 */
function createFormData(blob, reactionID) {
  blob;
  const formData = new FormData();
  formData.append("pdfFile", blob, `${reactionID}-summary.pdf`);
  const { workgroup, workbook } = getFormVariables();
  formData.append("workgroup", workgroup);
  formData.append("workbook", workbook);
  formData.append("reactionID", reactionID);
  return formData;
}

/**
 * Sends FormData to the backend.
 * @param {FormData} formData - The FormData to send.
 * @returns {Promise} A promise representing the completion of the backend request.
 */
async function sendFormDataToBackend(formData) {
  return fetch("/pdf", {
    method: "POST",
    body: formData,
  });
}

/**
 * Gets the current datetime for UK timezone and formats.
 * @return {string} - formatted datetime string.
 */
function getDateTime() {
  // Get the current date and time in UTC
  const currentUtcDate = new Date();
  // Specify the target timezone (e.g., 'Europe/London' for UK)
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
