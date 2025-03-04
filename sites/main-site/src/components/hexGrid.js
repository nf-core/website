/**
 * Hex Grid Layout
 * Arranges hexagonal elements in a honeycomb pattern
 * Inspired by hexb.in's implementation
 */

export function arrangeHexGrid(container, hexElements, options = {}) {
  const defaults = {
    spacing: 4, // Spacing between hexagons
    startingRow: 3, // Number of hexagons in the first row
    maxPerRow: null, // Maximum hexagons per row (calculated automatically if null)
    hexWidth: null, // Width of a hexagon (calculated from first element if null)
    hexHeight: null, // Height of a hexagon (calculated from first element if null)
    minHexPerRow: 2, // Minimum hexagons per row for small screens
    transitionDuration: '0.4s', // Duration for smooth transitions
  };

  const settings = { ...defaults, ...options };

  if (!container || !hexElements || hexElements.length === 0) return;

  // Calculate hex dimensions if not provided
  if (!settings.hexWidth || !settings.hexHeight) {
    const firstHex = hexElements[0];
    const rect = firstHex.getBoundingClientRect();
    settings.hexWidth = rect.width + settings.spacing;
    settings.hexHeight = rect.height + settings.spacing;
  }

  // Calculate the maximum number of hexagons per row based on container width
  if (!settings.maxPerRow) {
    const containerWidth = container.clientWidth;
    settings.maxPerRow = Math.max(Math.floor(containerWidth / settings.hexWidth), settings.minHexPerRow);

    // Adjust starting row for very small screens
    if (containerWidth < 576) { // Bootstrap's sm breakpoint
      settings.startingRow = Math.min(settings.startingRow, settings.minHexPerRow);
    }
  }

  // Calculate the horizontal offset for centering
  const hexHeight = settings.hexHeight;
  const hexWidth = settings.hexWidth;
  const verticalOverlap = hexHeight * 0.25; // Vertical overlap between rows

  // Arrange hexagons in a honeycomb pattern
  let currentRow = 0;
  let currentCol = 0;
  let hexPerRow = settings.startingRow;
  let increasing = true;
  let rowsCount = 0;

  // Calculate total rows needed
  const totalHexagons = hexElements.length;
  let remainingHexagons = totalHexagons;
  let tempHexPerRow = settings.startingRow;

  while (remainingHexagons > 0) {
    remainingHexagons -= tempHexPerRow;
    rowsCount++;

    if (tempHexPerRow >= settings.maxPerRow) {
      increasing = false;
    } else if (tempHexPerRow <= settings.startingRow && rowsCount > 1) {
      increasing = true;
    }

    tempHexPerRow = increasing ? tempHexPerRow + 1 : tempHexPerRow - 1;
    tempHexPerRow = Math.min(Math.max(tempHexPerRow, settings.startingRow), settings.maxPerRow);
  }

  // Set container height to accommodate all hexagons
  container.style.height = `${(rowsCount) * (hexHeight - verticalOverlap) + verticalOverlap}px`;
  container.style.position = 'relative';

  // Ensure all hexagons have transition property set for smooth movement
  hexElements.forEach(hex => {
    if (!hex.style.transition || !hex.style.transition.includes('left') || !hex.style.transition.includes('top')) {
      hex.style.transition = `left ${settings.transitionDuration} ease-out, top ${settings.transitionDuration} ease-out`;
    }
  });

  // Position each hexagon - use requestAnimationFrame to avoid layout thrashing
  requestAnimationFrame(() => {
    let hexIndex = 0;
    currentRow = 0;
    hexPerRow = settings.startingRow;
    increasing = true;

    while (hexIndex < hexElements.length) {
      // Calculate the horizontal offset for centering the row
      const rowWidth = hexPerRow * hexWidth;
      const offsetX = (container.clientWidth - rowWidth) / 2;

      // Position hexagons in the current row
      for (currentCol = 0; currentCol < hexPerRow && hexIndex < hexElements.length; currentCol++) {
        const hex = hexElements[hexIndex];

        // Calculate position
        const x = offsetX + currentCol * hexWidth;
        const y = currentRow * (hexHeight - verticalOverlap);

        // Apply position
        hex.style.position = 'absolute';
        hex.style.left = `${x}px`;
        hex.style.top = `${y}px`;

        hexIndex++;
      }

      // Update hexagons per row for next row
      if (hexPerRow >= settings.maxPerRow) {
        increasing = false;
      } else if (hexPerRow <= settings.startingRow && currentRow > 0) {
        increasing = true;
      }

      hexPerRow = increasing ? hexPerRow + 1 : hexPerRow - 1;
      hexPerRow = Math.min(Math.max(hexPerRow, settings.startingRow), settings.maxPerRow);
      currentRow++;
    }
  });
}

// Handle window resize with improved performance
export function setupHexGridResizing(container, hexElements, options = {}) {
  let resizeTimeout;
  let lastWidth = window.innerWidth;
  let lastHeight = window.innerHeight;
  let ticking = false;

  // More aggressive throttling for resize events
  function handleResize() {
    // Only process if dimensions actually changed significantly (at least 20px)
    const currentWidth = window.innerWidth;
    const currentHeight = window.innerHeight;

    const widthChanged = Math.abs(currentWidth - lastWidth) > 20;
    const heightChanged = Math.abs(currentHeight - lastHeight) > 20;

    if (!widthChanged && !heightChanged) return;

    lastWidth = currentWidth;
    lastHeight = currentHeight;

    // Clear any existing timeout
    clearTimeout(resizeTimeout);

    // Don't schedule multiple rAF callbacks
    if (!ticking) {
      ticking = true;

      // Use requestAnimationFrame for better performance
      requestAnimationFrame(() => {
        // Set a timeout to avoid too frequent updates
        resizeTimeout = setTimeout(() => {
          arrangeHexGrid(container, hexElements, options);
          ticking = false;
        }, 150); // Longer delay for smoother experience
      });
    }
  }

  // Debounced version for final layout after resize ends
  function handleResizeEnd() {
    clearTimeout(resizeTimeout);
    resizeTimeout = setTimeout(() => {
      arrangeHexGrid(container, hexElements, options);
    }, 250);
  }

  // Listen for both window resize and orientation change for mobile
  window.addEventListener('resize', handleResize, { passive: true });
  window.addEventListener('resize', handleResizeEnd, { passive: true });
  window.addEventListener('orientationchange', handleResizeEnd);

  // Initial arrangement
  arrangeHexGrid(container, hexElements, options);

  // Return cleanup function
  return () => {
    window.removeEventListener('resize', handleResize);
    window.removeEventListener('resize', handleResizeEnd);
    window.removeEventListener('orientationchange', handleResizeEnd);
  };
}

// Helper function to determine if a resize is significant enough to trigger rearrangement
export function shouldRearrange(oldWidth, newWidth, threshold = 0.05) {
  // Only rearrange if width changed by more than threshold percentage
  return Math.abs(oldWidth - newWidth) / oldWidth > threshold;
}
