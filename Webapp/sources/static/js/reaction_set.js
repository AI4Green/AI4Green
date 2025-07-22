// Initialize the page with the correct reactor type settings
window.onload = async function () {
  let reactorDimensions = JSON.parse($("#js-reactor-dimensions").val());
  await initialiseReactors(reactorDimensions); // Ensure the page is loaded with the correct display settings

  observer();
};

async function initialiseReactors(reactorDimensions) {
  let reactionSet = JSON.parse($("#js-reaction-set").val());

  // Sort reactions by reaction_id
  reactionSet.reactions.sort((a, b) => {
    if (a.reaction_id < b.reaction_id) return -1;
    if (a.reaction_id > b.reaction_id) return 1;
    return 0;
  });

  updateReactorType(reactorDimensions);
  assignReactions(reactionSet.reactions);

  // focus first well
  await focusWell("0");
}

let currentLoadToken = 0;

async function focusWell(index) {
  $(".focused").removeClass("focused");
  let wellId = `circle-${index}`;
  let well = $(`[data-id='${wellId}']`);
  if (!well.length) return;

  well.addClass("focused");

  const thisLoad = ++currentLoadToken;

  // Set load-status early so autosave will stop
  $("#js-load-status").val("loading");

  const reactionId = well.attr("reaction-id");

  const success = await loadReaction(reactionId, thisLoad);

  if (success && thisLoad === currentLoadToken) {
    $("#js-load-status").val("loaded");
  }
}

async function loadReaction(reactionID, loadToken) {
  const reactionSet = JSON.parse($("#js-reaction-set").val());
  const reactionIndicator = $("#js-reaction-id");

  // Guard against stale loads
  if (loadToken !== currentLoadToken) return false;

  const reaction = reactionSet.reactions.find(
    (r) => r.reaction_id === reactionID,
  );
  if (!reaction) return false;

  // Fill out all fields
  reactionIndicator.val(reactionID).text(reactionID);
  $("#js-reaction-smiles").val(reaction.reaction_smiles);
  $("#js-reaction-name").val(reaction.name);
  $("#js-summary-table-data").val(reaction.summary_table_data);
  $("#js-reaction-rxn").val(reaction.reaction_rxn);

  await reloadReactionTable(JSON.parse(reaction.reaction_table_data));
  await reloadSketcher();

  // Confirm still valid load
  return loadToken === currentLoadToken;
}

function assignReactions(reactions) {
  reactions.forEach((reaction, index) => {
    let wellId = `circle-${index}`;
    let well = $(`[data-id='${wellId}']`);
    if (well) {
      well.attr("reaction-id", reaction.reaction_id);
      well.attr("tabindex", "0"); // Make div focusable
      well.attr("title", reaction.reaction_id);
    }
  });
}

function createGrid(numCols, numRows) {
  const grid = document.getElementById("grid");
  const sidePanel = document.getElementById("side-panel");
  grid.innerHTML = ""; // Clear existing grid
  grid.style.gridTemplateRows = `repeat(${numRows}, 1fr)`;
  grid.style.gridTemplateColumns = `repeat(${numCols}, 1fr)`;
  let circleIndex = -1;

  for (let row = 0; row < numRows; row++) {
    for (let col = 0; col < numCols; col++) {
      const circle = document.createElement("div");
      circle.className = "btn icon";
      circle.style.width = "30px";
      circle.style.height = "30px";
      circle.dataset.id = `circle-${(circleIndex += 1)}`; // Assign a unique ID

      circle.onclick = () => focusWell(`${row}-${col}`);
      //circle.textContent = i + 1;

      // // Tooltip for reaction details
      // const tooltip = document.createElement('div');
      // tooltip.className = 'tooltip';
      // tooltip.textContent = getReactionDetails(row, col);
      // circle.appendChild(tooltip);
      //
      // // Add click event to focus circle
      // circle.addEventListener('click', () => {
      //     // Populate the side panel with reaction details
      //     sidePanel.innerHTML = `
      //         <div class="reaction-header">Reaction Constructor</div>
      //         <p>Details for well at row ${row + 1}, column ${col + 1}:</p>
      //         <p>This is a placeholder for the reaction constructor.</p>
      //     `;
      //     sidePanel.classList.add('active');
      // });

      grid.appendChild(circle);
    }
  }
}

function createCarousel(numberOfReactions) {
  const carousel = document.getElementById("carousel");
  const sidePanel = document.getElementById("side-panel");
  carousel.innerHTML = ""; // Clear existing carousel

  const { centerX, centerY, radius } = getCarouselDimensions();

  for (let i = 0; i < numberOfReactions; i++) {
    // Calculate angle for this point
    const angle = (2 * Math.PI * i) / numberOfReactions;

    // Calculate x and y positions for the point
    const x = centerX + radius * Math.cos(angle);
    const y = centerY + radius * Math.sin(angle);

    // Create the circle element
    const circle = document.createElement("div");
    circle.className = "btn icon";
    circle.style.position = "absolute";
    circle.style.fontsize = "30px";
    circle.style.left = `${x - centerX + radius - 42}px`;
    circle.style.top = `${y - centerY + radius - 42}px`;
    circle.dataset.id = `circle-${i}`; // Assign a unique ID
    circle.onclick = () => focusWell(i);
    circle.textContent = i + 1;

    // // Tooltip for reaction details
    // const tooltip = document.createElement("div");
    // tooltip.className = "tooltip";
    // tooltip.textContent = `Reaction at Point ${i + 1}`;
    // circle.appendChild(tooltip);
    //
    // // Add click event to focus on the circle
    // circle.addEventListener("click", () => {
    //   // Populate the side panel with reaction details
    //   sidePanel.innerHTML = `
    //     <div class="reaction-header">Reaction Constructor</div>
    //     <p>Details for Point ${i + 1}:</p>
    //     <p>This is a placeholder for the reaction constructor.</p>
    // `;
    //   sidePanel.classList.add("active");
    // });

    carousel.appendChild(circle);
  }
}

function getCarouselDimensions() {
  const carousel = document.getElementById("carousel");
  const rect = carousel.getBoundingClientRect(); // Get position and size of the container

  const radius = Math.min(rect.width, rect.height) / 2; // Radius is half the smaller dimension
  const centerX = rect.left + rect.width / 2; // X-coordinate of the center
  const centerY = rect.top + rect.height / 2; // Y-coordinate of the center

  return { centerX: centerX, centerY: centerY, radius: radius };
}

function updateReactorType(reactorDimensions) {
  // Get references to selectors and reactor containers
  const wellplateContainer = document.getElementById("wellplate-selector");
  const carouselContainer = document.getElementById("carousel-selector");
  const multiwellReactor = document.getElementById("wellplate-container");
  const carouselReactor = document.getElementById("carousel-container");
  let reactorType = reactorDimensions["reactorType"];
  if (reactorType === "well-plate") {
    // Show multiwell setup and reactor, hide carousel
    wellplateContainer.style.display = "block";
    multiwellReactor.style.display = "block";

    carouselContainer.style.display = "none";
    carouselReactor.style.display = "none";

    createGrid(reactorDimensions["columns"], reactorDimensions["rows"]);
  } else if (reactorType === "carousel") {
    // Show carousel setup and reactor, hide multiwell
    carouselContainer.style.display = "block";
    carouselReactor.style.display = "block";

    wellplateContainer.style.display = "none";
    multiwellReactor.style.display = "none";
    createCarousel(reactorDimensions["numberOfReactions"]); // Ensure carousel is created/updated
  }
}

function makeReactorDraggable(containerId, reactorId) {
  const container = document.getElementById(containerId); // Draggable area (modal container)
  const reactor = container.querySelector(`#${reactorId}`); // Reactor element within the container
  let selectables = [];
  const state = new Map(); // Tracks which circles are inside the draggable area

  function updateSelectables() {
    selectables = [];
    const selectableElems = [...reactor.querySelectorAll(".btn.icon")];

    for (const selectable of selectableElems) {
      const { x, y, width, height } = selectable.getBoundingClientRect();
      selectables.push({
        x: x + window.scrollX,
        y: y + window.scrollY,
        width,
        height,
        elem: selectable,
      });
    }
  }

  function toggleSelectionStateLive(selectAreaElem) {
    const select = selectAreaElem.getBoundingClientRect();
    const { x, y, height, width } = select;

    for (const selectable of selectables) {
      const isIntersecting = checkRectIntersection(
        {
          x: x + window.scrollX,
          y: y + window.scrollY,
          height,
          width,
        },
        selectable,
      );

      const id = selectable.elem.dataset.id;

      // If intersecting, toggle state
      if (isIntersecting && !state.get(id)) {
        selectable.elem.classList.toggle("selected");
        state.set(id, true);
      }

      // Reset state when not intersecting
      if (!isIntersecting && state.get(id)) {
        state.set(id, false);
      }
    }
  }

  function checkRectIntersection(r1, r2) {
    return !(
      r1.x + r1.width < r2.x ||
      r2.x + r2.width < r1.x ||
      r1.y + r1.height < r2.y ||
      r2.y + r2.height < r1.y
    );
  }

  function initializeDragSelect() {
    container.addEventListener("pointerdown", createSelectAreaDiv);

    function createSelectAreaDiv(event) {
      if (!reactor.contains(event.target)) return;

      updateSelectables(); // Update selectable items on drag start
      event.preventDefault();

      const x = event.pageX;
      const y = event.pageY;

      const div = document.createElement("div");
      div.style.position = "absolute";
      div.style.width = "0";
      div.style.height = "0";
      div.style.left = `${x}px`;
      div.style.top = `${y}px`;
      div.classList.add("drag-select");
      div.style.zIndex = 1060; // Ensure the drag area appears above the modal
      document.body.appendChild(div);

      state.clear(); // Reset selection state for new drag

      function resize(event) {
        const diffX = event.pageX - x;
        const diffY = event.pageY - y;

        div.style.left = diffX < 0 ? `${x + diffX}px` : `${x}px`;
        div.style.top = diffY < 0 ? `${y + diffY}px` : `${y}px`;
        div.style.width = Math.abs(diffX) + "px";
        div.style.height = Math.abs(diffY) + "px";

        toggleSelectionStateLive(div);
      }

      function stopResize() {
        removeEventListener("pointermove", resize);
        removeEventListener("pointerup", stopResize);
        div.remove();
      }

      addEventListener("pointermove", resize);
      addEventListener("pointerup", stopResize);
    }
  }

  initializeDragSelect();
}

function applyToWellModal() {
  let reactorDimensions = JSON.parse($("#js-reactor-dimensions").val());
  let reactorType = reactorDimensions["reactorType"];

  let modalReactorContainer = document.getElementById(
    "add-to-well-reactor-container",
  );
  modalReactorContainer.innerHTML = ""; // Clear existing modal content

  let clonedReactor;
  let modalReactor;
  if (reactorType === "well-plate") {
    // Clone the well plate reactor (adjust the ID or class of the well plate container)
    clonedReactor = document
      .getElementById("wellplate-container")
      .cloneNode(true);
    modalReactor = "modal-grid";
  } else if (reactorType === "carousel") {
    // Clone the carousel reactor
    clonedReactor = document
      .getElementById("carousel-container")
      .cloneNode(true);
    modalReactor = "modal-carousel";
  }

  // Rename IDs in the cloned reactor to ensure uniqueness
  function renameIds(element, prefix) {
    if (element.id) {
      element.id = prefix + element.id;
    }
    for (let child of element.children) {
      renameIds(child, prefix);
    }
  }
  renameIds(clonedReactor, "modal-"); // Add "modal-" prefix to all IDs

  modalReactorContainer.appendChild(clonedReactor);

  // Make the cloned reactor draggable
  makeReactorDraggable("add-to-well-reactor-container", modalReactor);

  // Show the modal
  $("#apply-to-well-modal").modal("show");
}

async function updateReactionSet() {
  // sync front and back end after saving multiple rows
  let setId = JSON.parse($("#js-reaction-set").val())["reaction_id"];
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  fetch("/update_reaction_set", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
    body: JSON.stringify({
      set_id: setId,
      workgroup: workgroup,
      workbook: workbook,
    }),
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      console.log("UPDATE SET");
      console.log(item);
    });
}

function saveReactions(reactions) {
  reactions.forEach((reaction) => {
    console.log(reaction.reaction_id);
    // set each reaction_id as active then post data for each reaction
    $("#js-reaction-id").val(reaction.reaction_id);
    $("#js-reaction-name").val(reaction.name);
    $("#js-summary-table-data").val(reaction.summary_table_data);

    postReactionData();
    // reload reaction set data from backend ----> IS THIS THE BEST IDEA?
  });
}

async function applyToAll() {
  let currentID = $("#js-reaction-id").val();
  let reactionSet = JSON.parse($("#js-reaction-set").val());
  saveReactions(reactionSet.reactions);
  await updateReactionSet();
  // reset reaction ID to what it was before save
  $("#js-reaction-id").val(currentID);
}

async function applyToSelectedWells() {
  let reactionSet = JSON.parse($("#js-reaction-set").val());
  let selectedReactions = reactionSet.reactions.filter((reaction) =>
    $("div.btn.icon.selected")
      .map(function () {
        return $(this).attr("reaction-id");
      })
      .get()
      .includes(reaction.reaction_id),
  );
  console.log(selectedReactions);

  saveReactions(selectedReactions);
  await updateReactionSet();
  $("#apply-to-well-modal").modal("hide");
}
