
// Mock database function to get reaction details
const getReactionDetails = (row, col) => {
    return `Reaction at row ${row + 1}, column ${col + 1}`;
};

function createGrid() {
    const cols = document.getElementById("column-size").value;
    const rows = document.getElementById("row-size").value;
    const grid = document.getElementById('grid');
    const sidePanel = document.getElementById('side-panel');
    grid.innerHTML = ''; // Clear existing grid
    grid.style.gridTemplateRows = `repeat(${rows}, 1fr)`;
    grid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`



    for (let row = 0; row < rows; row++) {
        for (let col = 0; col < cols; col++) {
            const circle = document.createElement('div');
            circle.className = 'btn icon';
            circle.style.width = "30px";
            circle.style.height = "30px";



            // Tooltip for reaction details
            const tooltip = document.createElement('div');
            tooltip.className = 'tooltip';
            tooltip.textContent = getReactionDetails(row, col);
            circle.appendChild(tooltip);

            // Add click event to focus circle
            circle.addEventListener('click', () => {
                // Populate the side panel with reaction details
                sidePanel.innerHTML = `
                    <div class="reaction-header">Reaction Constructor</div>
                    <p>Details for well at row ${row + 1}, column ${col + 1}:</p>
                    <p>This is a placeholder for the reaction constructor.</p>
                `;
                sidePanel.classList.add('active');
            });

            grid.appendChild(circle);
        }
    }
}

function createCarousel() {
    const numPoints = parseInt(document.getElementById("carousel-reactions").value, 10); // Get the number of points
    const carousel = document.getElementById('carousel');
    const sidePanel = document.getElementById('side-panel');
    carousel.innerHTML = ''; // Clear existing carousel

    const { centerX, centerY, radius } = getCarouselDimensions();

    for (let i = 0; i < numPoints; i++) {
        // Calculate angle for this point
        const angle = (2 * Math.PI * i) / numPoints;

        // Calculate x and y positions for the point
        const x = centerX + radius * Math.cos(angle);
        const y = centerY + radius * Math.sin(angle);

        // Create the circle element
        const circle = document.createElement('div');
        circle.className = 'btn icon';
        circle.style.position = 'absolute';
        circle.style.left = `${x - centerX + radius - 42}px`;
        circle.style.top = `${y - centerY + radius - 42}px`;

        // Tooltip for reaction details
        const tooltip = document.createElement('div');
        tooltip.className = 'tooltip';
        tooltip.textContent = `Reaction at Point ${i + 1}`;
        circle.appendChild(tooltip);

        // Add click event to focus on the circle
        circle.addEventListener('click', () => {
            // Populate the side panel with reaction details
            sidePanel.innerHTML = `
            <div class="reaction-header">Reaction Constructor</div>
            <p>Details for Point ${i + 1}:</p>
            <p>This is a placeholder for the reaction constructor.</p>
        `;
            sidePanel.classList.add('active');
        });

        carousel.appendChild(circle);
    }
}


function getCarouselDimensions() {
    const carousel = document.getElementById('carousel');
    const rect = carousel.getBoundingClientRect(); // Get position and size of the container

    const radius = Math.min(rect.width, rect.height) / 2; // Radius is half the smaller dimension
    const centerX = rect.left + rect.width / 2; // X-coordinate of the center
    const centerY = rect.top + rect.height / 2; // Y-coordinate of the center

    return { centerX: centerX, centerY: centerY, radius: radius };
}


function updateReactorType() {
    const reactorType = document.getElementById("reactor-type").value; // Get selected reactor type

    // Get references to selectors and reactor containers
    const wellplateContainer = document.getElementById("wellplate-selector");
    const carouselContainer = document.getElementById("carousel-selector");
    const multiwellReactor = document.getElementById("wellplate-container");
    const carouselReactor = document.getElementById("carousel-container");

    if (reactorType === "well-plate") {
        // Show multiwell setup and reactor, hide carousel
        wellplateContainer.style.display = "block";
        multiwellReactor.style.display = "block";

        carouselContainer.style.display = "none";
        carouselReactor.style.display = "none";
    } else if (reactorType === "carousel") {
        // Show carousel setup and reactor, hide multiwell
        carouselContainer.style.display = "block";
        carouselReactor.style.display = "block";

        wellplateContainer.style.display = "none";
        multiwellReactor.style.display = "none";

        createCarousel(); // Ensure carousel is created/updated
    }
}

// Initialize the page with the correct reactor type settings
window.onload = function() {
    createGrid(); // Ensure the page is loaded with the correct display settings
};

