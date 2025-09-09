function assignSolvents(unknownSolvents) {
  let numberOfSolvents = Number($("#js-number-of-unknown-solvents").val());
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let assignments = {};
  // get assignments from user input
  for (let i = 1; i <= numberOfSolvents; i++) {
    let unknownSolvent = $("#unknown-solvent-name-" + i)
      .text()
      .trim();
    let assignedSolvent = $("#js-solvent" + i).val();

    if (assignedSolvent) {
      // preserve experiments from backend dict
      assignments[unknownSolvent] = {
        assigned: assignedSolvent,
        reactions: unknownSolvents[unknownSolvent] || [],
      };
    }
  }
  fetch("/assign_reactwise_solvents", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      assignments: assignments,
      workbook_name: workbook,
      workgroup_name: workgroup,
    }),
  })
    .then((response) => response.json()) // or .text() depending on what the backend sends
    .then((data) => {
      updateSolventCards(data);
      updateTracker();
    })
    .catch((err) => console.error("Error:", err));
}

function updateSolventCards(data) {
  // Loop over successful solvents
  data.success.forEach((solventName) => {
    // Find the label with this solvent name
    const label = Array.from(
      document.querySelectorAll("[id^='unknown-solvent-name-']"),
    ).find((el) => el.textContent.trim() === solventName);

    if (label) {
      // Get the card container and remove it
      const card = label.closest(".col-md-6");
      if (card) card.remove();
    }
  });

  // Loop over failed solvents
  data.failed.forEach((solventName) => {
    const label = Array.from(
      document.querySelectorAll("[id^='unknown-solvent-name-']"),
    ).find((el) => el.textContent.trim() === solventName);

    if (label) {
      const card = label.closest(".card");
      if (card) {
        // Highlight failed card
        card.classList.add("border-danger");
        card.style.backgroundColor = "#ffe5e5"; // light red
        showCardToast(card, "Solvent not recognised");
      }
    }
  });

  // Check if any unknown solvent cards remain
  const unknownSolventsCard = document.getElementById("unknown-solvents-card");

  if (unknownSolventsCard) {
    // If no more unknown solvent cards remain, highlight the entire section
    const remainingCards = unknownSolventsCard.querySelectorAll(
      "[id^='unknown-solvent-card-']",
    );
    if (remainingCards.length === 0) {
      // Add green border and light green background
      unknownSolventsCard.classList.add("border-success");
      unknownSolventsCard.style.backgroundColor = "#e6f4ea"; // light green

      // Remove the button
      const saveButton = unknownSolventsCard.querySelector(
        ".btn-success[onclick^='assignSolvents']",
      );
      if (saveButton) saveButton.remove();

      // Show success message with tick
      const cardBody = unknownSolventsCard.querySelector(".card-body");
      if (cardBody) {
        const msg = document.createElement("p");
        msg.innerHTML = "&#10003; All solvents are assigned."; // ✔ tick
        msg.className = "text-success fw-bold mt-3";
        cardBody.appendChild(msg);
      }
    }
  }
}

function showCardToast(card, message) {
  // Make card relative if not already
  card.classList.add("position-relative");

  // Create toast element
  const toast = document.createElement("div");
  toast.className = "card-toast";
  toast.textContent = message;

  card.appendChild(toast);

  // Fade out after 2 seconds
  setTimeout(() => {
    toast.classList.add("fade-out");
  }, 2000);

  // Remove after 3 seconds
  setTimeout(() => {
    toast.remove();
  }, 3000);
}

function assignVariables() {
  let numberOfVariables = Number($("#js-number-of-unknown-variables").val());
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let reactionSetId = $("#js-active-reaction-set").val();
  let expDetails = JSON.parse($("#js-experimental-details").text()); // converts JSON string to JS object
  let assignments = {};
  // get assignments from user input, currently only does amounts
  for (let i = 1; i <= numberOfVariables; i++) {
    let label = $("#unknown-var-" + i);
    let unknownVariable = label.data("variable");
    let variableUnit = label.data("unit");
    // get better name for variable nature
    let variableNature = $("#variable-nature-" + i).val();
    let reactionComponent = $("#variable-assignment-" + i).val();
    let variableType = $("#variable-type-" + i).val();

    assignments[unknownVariable] = {
      nature: variableNature,
      component: reactionComponent,
      type: variableType,
      unit: variableUnit,
    };
  }

  fetch("/assign_reactwise_variables", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      assignments: assignments,
      workbook_name: workbook,
      workgroup_name: workgroup,
      reaction_set_id: reactionSetId,
      exp_details: expDetails,
    }),
  })
    .then((response) => response.json()) // or .text() depending on what the backend sends
    .then((data) => {
      updateVariableCards(data);
      updateTracker();
    })
    .catch((err) => console.error("Error:", err));
}

function updateVariableCards(data) {
  // Remove successful variables
  data.success.forEach((variableName) => {
    const label = Array.from(
      document.querySelectorAll("[id^='unknown-var-']"),
    ).find((el) => el.dataset.variable === variableName);

    if (label) {
      const card = label.closest(".col-md-6");
      if (card) card.remove();
    }
  });

  // Mark failed variables
  data.failed.forEach((variableName) => {
    const label = Array.from(
      document.querySelectorAll("[id^='unknown-var-']"),
    ).find((el) => el.dataset.variable === variableName);

    if (label) {
      const card = label.closest(".card");
      if (card) {
        card.classList.add("border-danger");
        card.style.backgroundColor = "#ffe5e5";
        showCardToast(card, "Variable not recognised");
      }
    }
  });

  // Section-level success check
  const unknownVariablesCard = document.getElementById(
    "unknown-variables-card",
  );
  if (unknownVariablesCard) {
    const remainingCards = unknownVariablesCard.querySelectorAll(
      "[id^='unknown-var-']",
    );
    if (remainingCards.length === 0) {
      unknownVariablesCard.classList.add("border-success");
      unknownVariablesCard.style.backgroundColor = "#e6f4ea";

      // remove save button
      const saveButton = unknownVariablesCard.querySelector(
        ".btn-success[onclick^='assignVariables']",
      );
      if (saveButton) saveButton.remove();

      // success message
      const cardBody = unknownVariablesCard.querySelector(".card-body");
      if (cardBody) {
        const msg = document.createElement("p");
        msg.innerHTML = "&#10003; All variables are assigned.";
        msg.className = "text-success fw-bold mt-3";
        cardBody.appendChild(msg);
      }
    }
  }
}

function postNovelCompound(index) {
  // Select fields by index
  let name = $(`#js-new-compound-name-${index}`).val();
  let molWeight = $(`#js-new-compound-mw-${index}`).val();
  let hPhrase = $(`#js-new-compound-hazards-${index}`).val();
  let cas = $(`#js-new-compound-cas-${index}`).val();
  let density = $(`#js-new-compound-density-${index}`).val();
  let concentration = $(`#js-new-compound-concentration-${index}`).val();
  let smiles = $(`#js-new-compound-smiles-${index}`).val();
  let polymer = $(`#js-polymer-${index}`).val();
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();

  let requestData = {
    name: name,
    molWeight: molWeight,
    hPhrase: hPhrase,
    density: density,
    concentration: concentration,
    smiles: smiles,
    polymer: polymer,
    workbook: workbook,
    workgroup: workgroup,
    component: "reactant",
    source: "card",
  };

  if (!polymer) {
    requestData.cas = cas;
  }

  // Send AJAX request
  $.ajax({
    url: "/_novel_compound",
    type: "POST",
    dataType: "json",
    data: requestData,
  }).done(function (data) {
    console.log(data);
    alert(data.feedback);
    if (data.feedback === "Compound added to the database") {
      removeNovelCompoundCard(index);
      updateTracker();
    }
  });
}

function removeNovelCompoundCard(index) {
  // Find the specific form container by index
  const formContainer = document.getElementById(
    `js-novel-compound-input-form-${index}`,
  );
  if (!formContainer) return;

  // Remove the outer card of this form
  const card = formContainer.closest(".card.p-3");
  if (card) card.remove();

  // Get the section/container holding all novel compound cards
  const novelCompoundsSection = document.getElementById(
    "novel-compounds-section",
  );
  if (!novelCompoundsSection) return;

  // Check if any cards are left
  const remainingCards = novelCompoundsSection.querySelectorAll(".card.p-3");
  if (remainingCards.length === 0) {
    // Update section styling to indicate all compounds processed
    novelCompoundsSection.classList.add("border-success");
    novelCompoundsSection.style.backgroundColor = "#e6f4ea";

    // Remove save button if present
    const saveButton = novelCompoundsSection.querySelector(
      ".btn-success[onclick^='saveNovelCompounds']",
    );
    if (saveButton) saveButton.remove();

    // Add success message
    const cardBody = novelCompoundsSection.querySelector(".card-body");
    if (cardBody) {
      const msg = document.createElement("p");
      msg.innerHTML = "&#10003; All novel compounds have been processed.";
      msg.className = "text-success fw-bold mt-3";
      cardBody.appendChild(msg);
    }
  }
}

function updateTracker() {
  // Count unresolved cards in each section
  const solventCards = document.querySelectorAll(
    '#unknown-solvents-card .card[tabindex="0"]',
  );
  const variableCards = document.querySelectorAll(
    '#unknown-variables-card .card[tabindex="0"]',
  );
  const compoundCards = document.querySelectorAll(
    '#novel-compounds-section .card[tabindex="0"]',
  );

  const unresolved =
    solventCards.length + variableCards.length + compoundCards.length;

  // Get tracker elements
  const icon = document.getElementById("tracker-icon");
  const text = document.getElementById("tracker-text");
  const button = document.getElementById("tracker-button");

  if (!icon || !text || !button) return; // safety check

  if (unresolved > 0) {
    icon.innerHTML = "&#10060;"; // Red cross
    icon.style.color = "red";
    text.textContent = "There are still unrecognised fields";
    button.style.display = "none";
  } else {
    icon.innerHTML = "&#9989;"; // Green tick
    icon.style.color = "green";
    text.textContent = "All fields recognised!";
    button.style.display = "inline-block";
  }
}
