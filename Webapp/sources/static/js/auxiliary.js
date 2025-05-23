window.addEventListener("load", function () {
  checkCookie();
});

function checkCookie() {
  // function to check cookies preferences
  let cookie = getCookie("cookies_preference");
  if (cookie === "accepted") {
    // hide banner if user has previously accepted
    $("#cookie-consent-banner").hide();
  } else {
    $("#cookie-consent-banner").show();
  }
}

function updateCookiePreferences() {
  // function to update/save cookie preferences
  document.cookie = "cookies_preference=accepted";
  checkCookie();
}

function updatePrivacyPolicyAgreement() {
  // patch request to update the date when the user agreed to the privacy policy, then redirect home.
  fetch("/accept_privacy_policy", {
    method: "PATCH",
  }).then((response) => {
    if (response.ok) {
      window.location.href = "/home";
    }
  });
}

function getCookie(cname) {
  // retrieves cookie values from cookie object
  let name = cname + "=";
  let decodedCookie = decodeURIComponent(document.cookie);
  let ca = decodedCookie.split(";");
  for (let i = 0; i < ca.length; i++) {
    let c = ca[i];
    while (c.charAt(0) == " ") {
      c = c.substring(1);
    }
    if (c.indexOf(name) == 0) {
      return c.substring(name.length, c.length);
    }
  }
  return "";
}

/**
 *  Uses JQuery to get the value from an element and escape special characters to prevent XSS vulnerabilities
 * @param jquerySelector {jQuery | HTMLElement | String}- The JQuery Selector used to get the value
 * @returns {String} - The value of HTML element with special characters escaped
 */
function getVal(jquerySelector) {
  const unescapedInput = $(jquerySelector).val();
  return $("<div>").text(unescapedInput).html();
}

/**
 * Uses JQuery to get the value of a data-* attribute from a HTML element and escapes special characters
 * @param jquerySelector {jQuery | HTMLElement | String} - The jQuery Selector used to get the value
 * @param dataAttribute {String} - The specific data attribute e.g., to access data-workbook arg should be 'workbook'
 * @returns {String} - The data-* value of the HTML element with special characters escaped
 */
function getData(jquerySelector, dataAttribute) {
  const unescapedInput = $(jquerySelector).data(dataAttribute);
  return $("<div>").text(unescapedInput).html();
}

/**
 *  Uses JQuery to get the numerical value from an element and escape special characters to prevent XSS vulnerabilities
 * @param jquerySelector {jQuery | HTMLElement | String}- The JQuery Selector used to get the value
 * @returns {Number} - The numerical value of the HTML element with special characters escaped
 */
function getNum(jquerySelector) {
  return Number(getVal(jquerySelector));
}

/**
 * Returns the limiting reactant table number. The user changes this by clicking a radio button
 * @returns {string} - The string of the number, used to form HTML element IDs.
 */
function getLimitingReactantTableNumber() {
  return getVal($("input[name='reactant-limiting']:checked"));
}

function updateSelectedWorkGroup(origin_page = "") {
  // updates the workbook+creator dropdowns after change to selected workgroup
  let workgroup = $("#active-workgroup").val();
  // Return a promise that resolves when the AJAX request completes
  return new Promise((resolve, reject) => {
    $.ajax({
      url: "/updated_workgroup_dropdown",
      type: "post",
      datatype: "json",
      data: { workgroup: workgroup, origin_page: origin_page },
      success: function (response) {
        // update dropdown with workbook options
        let dropdownSelect = $("#active-workbook");
        dropdownSelect.empty();
        for (let workbook of response.workbooks) {
          let option = document.createElement("option");
          option.text = workbook;
          option.value = workbook;
          dropdownSelect.append(option);
        }
        // Resolve the promise once the AJAX request is successful
        resolve();
      },
      error: function (xhr, status, error) {
        // Reject the promise if there's an error in the AJAX request
        reject(error);
      },
    });
  });
}

/**
 * Capitalizes the first character of a string and returns the result.
 * @param {string} str - The input string to capitalize
 * @returns {string} The string with its first character capitalized
 */
function capitalize(str) {
  return str.charAt(0).toUpperCase() + str.slice(1);
}

/**
 * Converts the first character of a string to lowercase and returns the result.
 * @param {string} str - The input string to modify
 * @returns {string} The string with its first character in lowercase
 */
function lower(str) {
  return str.charAt(0).toLowerCase() + str.slice(1);
}
