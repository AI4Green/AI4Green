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
