window.addEventListener('load', function () {
  checkCookie();
})

function checkCookie(){
    // function to check cookies preferences
    let cookie = getCookie('cookies_preference');
    if (cookie === 'accepted') {
        // hide banner if user has previously accepted
        $("#cookie-consent-banner").hide();
    } else {
        $("#cookie-consent-banner").show();
        }
}

function updateCookiePreferences(){
    // function to update/save cookie preferences
    document.cookie = 'cookies_preference=accepted';
    checkCookie();
}

function getCookie(cname) {
    // retrieves cookie values from cookie object
  let name = cname + "=";
  let decodedCookie = decodeURIComponent(document.cookie);
  let ca = decodedCookie.split(';');
  for(let i = 0; i <ca.length; i++) {
    let c = ca[i];
    while (c.charAt(0) == ' ') {
      c = c.substring(1);
    }
    if (c.indexOf(name) == 0) {
      return c.substring(name.length, c.length);
    }
  }
  return "";
}
