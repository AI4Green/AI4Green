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

function getSchemata(sort_crit, workbook, workgroup, size){
    return new Promise(function(resolve, reject){
        // post to get_schemata and get the schemes for reaction images
        $.ajax({
            method: "POST",
            url: "/get_schemata",
            dataType: 'json',
            data: {
                sort_crit: sort_crit,
                workgroup: workgroup,
                workbook: workbook,
                size: size
            },
            success: function (data) {
                resolve(data.schemes)
            },
            error: function () {
                reject("error accessing /get_schemata")
            }
        });
    });
}

function sortReactionsAlphabetically(){
    let reactionCards = $(".reaction-card")
    reactionCards.sort((a, b) => a.querySelector(".reaction-name").firstChild.innerHTML.toLowerCase().localeCompare(b.querySelector(".reaction-name").firstChild.innerHTML.toLowerCase()))
    $("#reaction-list").empty()
    $("#reaction-list").append(reactionCards)
}

function sortReactionsByTime(){
    let reactionCards = $(".reaction-card")
    reactionCards.sort(function(a, b) {return new Date(b.querySelector(".reaction-time").innerHTML) - new Date(a.querySelector(".reaction-time").innerHTML)})
    $("#reaction-list").empty()
    $("#reaction-list").append(reactionCards)
}

function getpdf(){
    let workbook = $("#WB-select").val();
    let workgroup = $('#workgroup_selected').val();
    let sort_crit = $("#js-sort-crit").val();
    window.open("/export_data_pdf/" + workgroup + "/" + workbook + "/" + sort_crit, '_blank').focus();
}

function getcsv(){
    let workbook = $("#WB-select").val();
    let workgroup = $('#workgroup_selected').val();
    window.location ="/export_data_csv/" + workgroup + "/" + workbook;
}