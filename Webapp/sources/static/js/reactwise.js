function assignSolvents(unknownSolvents) {
  let numberOfSolvents = Number($("#js-number-of-unknown-solvents").val());
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let assignments = {};
  console.log(workgroup, workbook);
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
  fetch("/reactwise_assign_solvents", {
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
        ".btn-primary[onclick^='assignSolvents']",
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
