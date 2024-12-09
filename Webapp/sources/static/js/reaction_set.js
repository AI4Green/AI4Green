const grid = document.getElementById('grid');
const gridSizeSelector = document.getElementById('grid-size');
const sidePanel = document.getElementById('side-panel');

// Mock database function to get reaction details
const getReactionDetails = (row, col) => {
    return `Reaction at row ${row + 1}, column ${col + 1}`;
};

const createGrid = (rows, cols) => {
    grid.innerHTML = ''; // Clear existing grid
    grid.style.gridTemplateRows = `repeat(${rows}, 1fr)`;
    grid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;

    for (let row = 0; row < rows; row++) {
        for (let col = 0; col < cols; col++) {
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
};

// Initial grid setup
const [initialRows, initialCols] = gridSizeSelector.value.split('x').map(Number);
createGrid(initialRows, initialCols);

// Update grid on dropdown change
gridSizeSelector.addEventListener('change', () => {
    const [rows, cols] = gridSizeSelector.value.split('x').map(Number);
    createGrid(rows, cols);
});
