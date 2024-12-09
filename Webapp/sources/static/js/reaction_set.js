
// Mock database function to get reaction details
const getReactionDetails = (row, col) => {
    return `Reaction at row ${row + 1}, column ${col + 1}`;
};

const createGrid = () => {
    const cols = document.getElementById("column-size").value;
    const rows = document.getElementById("row-size").value;
    const grid = document.getElementById('grid');
    const sidePanel = document.getElementById('side-panel');
    grid.innerHTML = ''; // Clear existing grid
    grid.style.gridTemplateRows = `repeat(${rows}, 1fr)`;
    grid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
    console.log(cols, rows)

    for (let row = 0; row < rows; row++) {
        for (let col = 0; col < cols; col++) {
            console.log("s")
            console.log(row, col)
            const circle = document.createElement('div');
            circle.className = 'btn icon';

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

    const createCarousel = () => {
    const numPoints = parseInt(document.getElementById("carousel-reactions").value, 10); // Get the number of points
    const carousel = document.getElementById('carousel');
    const sidePanel = document.getElementById('side-panel');
    carousel.innerHTML = ''; // Clear existing carousel

    const radius = 150; // Radius of the circle in pixels
    const centerX = 200; // X-coordinate of the center
    const centerY = 200; // Y-coordinate of the center

    for (let i = 0; i < numPoints; i++) {
        // Calculate angle for this point
        const angle = (2 * Math.PI * i) / numPoints;

        // Calculate x and y positions for the point
        const x = centerX + radius * Math.cos(angle);
        const y = centerY + radius * Math.sin(angle);

        // Create the circle element
        const circle = document.createElement('div');
        circle.className = 'carousel-circle';
        circle.style.position = 'absolute';
        circle.style.left = `${x}px`;
        circle.style.top = `${y}px`;

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
};

    // // Update grid on dropdown change
    // gridSizeSelector.addEventListener('change', () => {
    //     const [rows, cols] = gridSizeSelector.value.split('x').map(Number);
    //     createGrid(rows, cols);
    // });
};

