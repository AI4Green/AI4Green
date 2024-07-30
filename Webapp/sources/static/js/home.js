function focusIcon(type, btnId) {
    document.querySelectorAll('.btn.icon.' + type).forEach(btn => {
                btn.classList.remove('focused');
            });

    button = document.getElementById(btnId)
    button.classList.add("focused")
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
           var activeIcons = findActiveWorGroupAndWorkbook()
           window.location.href = "/sketcher/" + activeIcons[0] + "/" + activeIcons[1] + "/" + input +"/no";
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

function findActiveWorGroupAndWorkbook() {
    // Select all elements with both 'icon' and 'focused' classes
            var focusedIcons = document.querySelectorAll('.icon.focused');
            var workgroup, workbook;

            // Iterate through the focused icons and check their classes
            focusedIcons.forEach(function(icon) {
                if (icon.classList.contains('workgroup')) {
                    workgroup = icon.id;
                }
                if (icon.classList.contains('workbook')) {
                    workbook = icon.id;
                }
            });
            return [workgroup, workbook]
}

function loadIcons(selected, type){
    // loads workbooks in workgroup selected in homepage
    // only needs selected workgroup to load active reactions
    var activeIcons = findActiveWorGroupAndWorkbook()
    fetch("/load_icons", {
    headers: {
        "Content-Type": "application/json",
      },
      method: "POST",
      body: JSON.stringify({
        input: selected,
        load_type: type,
        activeWorkgroup: activeIcons[0]
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
    activeIcons = findActiveWorGroupAndWorkbook()

    fetch("/get_new_reaction_id", {
    headers: {
        "Content-Type": "application/json",
      },
      method: "POST",
      body: JSON.stringify({
        workgroup: activeIcons[0],
        workbook: activeIcons[1],
      })
    })
    .then(function (response) { return response.json() })
    .then(function (item) {
      $("#new-reaction-id-input").val(item);
      $("#new-reaction-name").val("");
      $("#error-warning-new-reaction").html("");
      $("#new-reaction-modal").modal("show");
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