function focusIcon(type, btnId) {
    document.querySelectorAll('.btn.icon.' + type).forEach(btn => {
                btn.classList.remove('focused');
            });

    document.getElementById(btnId).classList.add("focused")
}

function blurIcon(button) {
    button.classList.remove('focused');
}

function toggleIcon(input, type, btnId) {
    button = document.getElementById(btnId)
    if (button.classList.contains('focused')) {
        blurIcon(button)
        unloadIconPanel(type)
    }

    else {
        if (type != "reaction") {
            focusIcon(type, btnId)
            loadIcons(input, type)
        }
        else {
        var workgroup = document.getElementsByClassName("btn icon workgroup focused")[0].id
        var workbook = document.getElementsByClassName("btn icon workbook focused")[0].id
           window.location.href = "/sketcher/" + workgroup + "/" + workbook + "/" + input +"/no";
        }
    }
}

function unloadIconPanel(type){
    if (type == "workgroup") {
        $("#workbook-icons").fadeOut("slow");
        $("#reaction-icons").fadeOut("slow");
    }

    else if (type == "workbook") {
        $("#reaction-icons").fadeOut("slow");
    }

}

function loadIcons(selected, type){
    // loads workbooks in workgroup selected in homepage
    // only needs selected workgroup to load active reactions
    let selectedWorkgroup = document.getElementsByClassName("btn icon workgroup focused")[0].id
    fetch("/load_icons", {
    headers: {
        "Content-Type": "application/json",
      },
      method: "POST",
      body: JSON.stringify({
        input: selected,
        load_type: type,
        activeWorkgroup: selectedWorkgroup
      })
    })
    .then(function(response) {return response.json()} )
    .then(function(item) {
        if (type == "workgroup"){
            $("#workbook-icons").html(item).hide().fadeIn("slow");
            $("#reaction-icons").html(item).hide().fadeOut("slow");
        }
        else if (type == "workbook"){
            $("#reaction-icons").html(item).hide().fadeIn("slow");
        }
    })
}

function newReactionSetup() {
    var activeWorkgroup = document.getElementsByClassName("btn icon workgroup focused")[0].id
    var activeWorkbook = document.getElementsByClassName("btn icon workbook focused")[0].id

     fetch("/get_new_reaction_id", {
    headers: {
        "Content-Type": "application/json",
      },
      method: "POST",
      body: JSON.stringify({
        workgroup: activeWorkgroup,
        workbook: activeWorkbook,
      })
    })
    .then(function (response) { return response.json() })
    .then(function (item) {
      $("#new-reaction-id-input").val(item);
      $("#new-reaction-name").val("");
      $("#error-warning-new-reaction").html("");
    })
}

function newType(type) {
    if (type == "workgroup") {
        window.location.href = "/create_workgroup"
    }
    else if (type == "workbook") {
        var workgroup = document.getElementsByClassName("btn icon workgroup focused")[0].id
        window.location.href = "/create_workbook/" + workgroup
    }

    else {
        newReactionSetup()
    }
}