/**
 * Focuses on a specific icon button by adding a "focused" class to it and
 * removing the "focused" class from all other buttons of the same type.
 *
 * @param {string} type - The type of the icon buttons to be focused. Should be either "workgroup", "workbook" or "reaction"
 *
 * @param {string} btnId - The ID of the button to be focused.
 *
 */
function focusIcon(type, btnId) {
  document.querySelectorAll(".btn.icon." + type).forEach((btn) => {
    btn.classList.remove("focused");
  });

  button = document.getElementById(btnId);
  button.classList.add("focused");
}

/**
 * Removes the "focused" class from a specified button.
 *
 * @param {HTMLElement} button - The button element from which the "focused" class will be removed.
 *
 */
function blurIcon(button) {
  button.classList.remove("focused");
}

/**
 * Toggles the "focused" state of a button and performs additional actions based on the type of icon and its current state.
 * If the button is already focused, it blurs the button and unloads the icon panel.
 * If the button is not focused, it either focuses the button and loads icons or redirects to the sketcher based on the type.
 *
 * @param {string} input - label for the icon. Either Workgroup name, Workbook name or Reaction code
 * @param {string} type - The type of the icon buttons. This determines the behavior of the function:
 *                        - If not "reaction", it focuses the button and loads icons.
 *                        - If "reaction", it triggers a URL redirect.
 * @param {string} btnId - The ID of the button to be toggled.
 *
 */
function toggleIcon(input, type, btnId) {
  let button = document.getElementById(btnId);
  if (button.classList.contains("focused")) {
    blurIcon(button);
    unloadIconPanel(type);
  } else {
    if (type != "reaction") {
      focusIcon(type, btnId);
      if (type != "reactwise") {
        loadIcons(input, type);
      }
    } else {
      var activeIcons = findActiveWorkgroupAndWorkbook();
      window.location.href =
        "/sketcher/" +
        activeIcons[0] +
        "/" +
        activeIcons[1] +
        "/" +
        input +
        "/no";
    }
  }
}

/**
 * Unloads the icon panel based on the specified type by fading out the corresponding icons.
 *
 * @param {string} type - The type of icon panel to unload. This determines which elements are faded out:
 *                        - "workgroup": Fades out both the workbook icons and reaction icons.
 *                        - "workbook": Fades out only the reaction icons.
 *
 */
function unloadIconPanel(type) {
  if (type == "workgroup") {
    $("#workbook-icons").fadeOut("slow");
    $("#reaction-icons").fadeOut("slow");
  } else if (type == "workbook") {
    $("#reaction-icons").fadeOut("slow");
  }
}

/**
 * Finds and returns the IDs of the currently focused workgroup and workbook icons.
 * The function looks for elements with both 'icon' and 'focused' classes, and then checks
 * their specific types (workgroup and workbook) to determine their IDs.
 *
 * @returns {[string | undefined, string | undefined]} An array containing the IDs of the focused
 *                                                    workgroup and workbook icons. The first element
 *                                                    of the array is the workgroup ID and the second
 *                                                    is the workbook ID. If an icon of a particular type
 *                                                    is not found, its corresponding array element will
 *                                                    be `undefined`.
 */
function findActiveWorkgroupAndWorkbook() {
  // Select all elements with both 'icon' and 'focused' classes
  var focusedIcons = document.querySelectorAll(".icon.focused");
  var workgroup, workbook;

  // Iterate through the focused icons and check their classes
  focusedIcons.forEach(function (icon) {
    if (icon.classList.contains("workgroup")) {
      workgroup = icon.id;
    }
    if (icon.classList.contains("workbook")) {
      workbook = icon.id;
    }
  });
  return [workgroup, workbook];
}

/**
 * Loads and displays icons based on the specified type and selected workgroup or workbook.
 * Sends a POST request to load icons and updates the UI with the response.
 *
 * @param {string} selected - The ID or identifier of the selected workgroup or workbook.
 * @param {string} type - The type of icons to load. Determines which elements are updated:
 *                        - "workgroup": Loads and displays workbook icons, hides reaction icons.
 *                        - "workbook": Loads and displays reaction icons.
 *
 */
function loadIcons(selected, type) {
  // loads workbooks in workgroup selected in homepage
  // only needs selected workgroup to load active reactions
  var activeIcons = findActiveWorkgroupAndWorkbook();
  fetch("/load_icons", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
    body: JSON.stringify({
      input: selected,
      load_type: type,
      activeWorkgroup: activeIcons[0],
    }),
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      if (type == "workgroup") {
        $("#workbook-icons").html(item).hide().fadeIn("slow");
        $("#reaction-icons").html(item).hide().fadeOut("slow");
      } else if (type == "workbook") {
        $("#reaction-icons").html(item).hide().fadeIn("slow");
      }
    });
}

/**
 * Sets up a new reaction by fetching a new reaction ID and loading the new reaction modal which deals with the creation
 *  of the new reaction
 *
 */
function newReactionSetup() {
  activeIcons = findActiveWorkgroupAndWorkbook();
  $("#active-workgroup").val(activeIcons[0]);
  $("#active-workbook").val(activeIcons[1]);

  fetch("/get_new_reaction_id", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
    body: JSON.stringify({
      workgroup: activeIcons[0],
      workbook: activeIcons[1],
    }),
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      $("#new-reaction-id-input").val(item);
      $("#new-reaction-name").val("");
      $("#error-warning-new-reaction").html("");
      $("#new-reaction-modal").modal("show");
    });
}

/**
 * Redirects to new workgroup/new workbook URL or triggers a new reaction setup based on the type parameter.
 *
 * @param {string} type - The type of action to perform. Determines the following actions:
 *                        - "workgroup": Redirects to the workgroup creation page.
 *                        - "workbook": Redirects to the workbook creation page for the currently focused workgroup.
 *                        - Any other value: Triggers the new reaction setup.
 *
 */
function newType(type) {
  if (type == "workgroup") {
    window.location.href = "/create_workgroup";
  } else if (type == "workbook") {
    var workgroup = document.getElementsByClassName(
      "btn icon workgroup focused",
    )[0].id;
    window.location.href = "/create_workbook/" + workgroup;
  } else {
    newReactionSetup();
  }
}
