# nf-core Website Components

This directory contains the Astro components used in the nf-core website.

## HexGrid Component

The `HexGrid.astro` component is used to arrange hexagonal elements in a honeycomb pattern. It replaces the previous JavaScript implementation with a more integrated Astro component.

### Usage

```astro
---
import HexGrid from "./components/HexGrid.astro";
import HexSticker from "./components/HexSticker.astro";
---

<HexGrid class="custom-class" id="my-hex-grid">
  <HexSticker imageUrl="/path/to/image.png" alt="Sticker description" />
  <HexSticker imageUrl="/path/to/image2.png" alt="Another sticker" />
  <!-- Add more hex stickers as needed -->
</HexGrid>
```

### Props

| Prop | Type | Description |
|------|------|-------------|
| `class` | `string` | Optional CSS class to add to the container |
| `id` | `string` | Optional ID for the container |

### Features

- Responsive layout that adjusts to different screen sizes
- Automatic arrangement of hexagonal elements in a honeycomb pattern
- Smooth transitions when resizing the window
- Bootstrap-compatible styling
- Optimized performance with throttled resize handling

### CSS Variables

The component uses the following CSS variables that can be customized:

- `--hex-size`: Controls the size of the hexagons (default varies by screen size)
- `--hex-margin`: Controls the margin around each hexagon (default: 0.5rem)

## HexSticker Component

The `HexSticker.astro` component represents an individual hexagonal sticker in the grid.

### Props

| Prop | Type | Description |
|------|------|-------------|
| `imageUrl` | `string` | URL to the sticker image |
| `alt` | `string` | Alt text for the image |
| `category` | `string` | Category of the sticker (optional) |
| `link` | `string` | Optional link for the sticker |
| `description` | `string` | Optional description for the sticker |

## Implementation Details

The HexGrid component uses client-side JavaScript to arrange the hexagonal elements in a honeycomb pattern. It calculates the optimal layout based on the container width and the number of elements.

The arrangement is done using absolute positioning, with each element positioned at specific coordinates to create the honeycomb pattern. The component also handles window resizing, ensuring that the layout is updated when the window size changes.

The component exposes three functions globally:

- `arrangeHexGrid`: Arranges the hexagonal elements in a honeycomb pattern
- `setupHexGridResizing`: Sets up resize handling for the hex grid
- `shouldRearrange`: Helper function to determine if a resize is significant enough to trigger rearrangement

These functions are used internally by the component but can also be used by other components if needed.
