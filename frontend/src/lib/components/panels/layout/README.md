# Panel Layout System - Carousel Documentation

The panel layout system provides responsive UI layouts with a sophisticated carousel system for navigation between panels. This document provides comprehensive coverage of the carousel system, which has been a critical component requiring careful handling.

## Architecture Overview

```
layout/
├── PanelCarousel.svelte   # Main carousel component (Swiper.js)
├── PanelGrid.svelte       # Grid layout for desktop
├── PanelOverlay.svelte    # Modal/overlay panels
├── index.js               # Layout exports
└── README.md              # This documentation
```

## Carousel System Deep Dive

### Design Philosophy

The carousel system uses a **focus-driven architecture** where `focusedPanel` is the single source of truth for all navigation state. This eliminates synchronization issues common in traditional index-based carousels.

```
State Flow: focusedPanel → carouselMode → currentIndex → UI Updates
```

### Core State Management

#### Primary State

```typescript
// Single source of truth
export const focusedPanel = writable<string | null>(null);

// Derived states (automatically computed)
export const currentCarouselMode = derived(focusedPanel, ($focused) => $focused !== null);

export const currentCarouselIndex = derived(
  [focusedPanel, currentPanels],
  ([$focused, $panels]) => {
    return $panels.findIndex((panel) => panel.id === $focused);
  }
);
```

#### Page-Specific Panel Orders

```typescript
const PANEL_ORDERS = {
  'spectral-simulation': ['spectrum', 'fragmentation', 'structure', 'peaks', 'console', 'export'],
  'how-it-works': [
    'system-overview',
    'workflow',
    'ml-engine',
    'fragmentation-science',
    'accuracy-validation',
    'technical-requirements',
    'visualization-features',
    'advanced-features',
  ],
  about: [
    'developer-hero',
    'developer-journey',
    'project-metrics',
    'tech-stack',
    'features-showcase',
    'contact-connect',
    'acknowledgments',
    'open-source',
  ],
  'chat-with-spectrum': ['chat-conversation', 'chat-spectrum', 'compound-stats'],
} as const;
```

### Navigation Functions

#### Core Navigation API

```typescript
// Enter/exit carousel mode
export function setCarouselMode(pageKey: PageKey, mode: boolean) {
  if (mode) {
    const firstPanelId = getFirstPanelIdForPage(pageKey);
    focusedPanel.set(firstPanelId);
  } else {
    focusedPanel.set(null);
  }
}

// Navigate panels
export function nextCarouselPanel(pageKey: PageKey) {
  const panels = getPanelIdsForPage(pageKey);
  const currentPanelId = get(focusedPanel);
  const currentIdx = panels.findIndex((id) => id === currentPanelId);
  const nextIdx = (currentIdx + 1) % panels.length; // Wraps around
  focusedPanel.set(panels[nextIdx]);
}

export function prevCarouselPanel(pageKey: PageKey) {
  const panels = getPanelIdsForPage(pageKey);
  const currentPanelId = get(focusedPanel);
  const currentIdx = panels.findIndex((id) => id === currentPanelId);
  const prevIdx = currentIdx === 0 ? panels.length - 1 : currentIdx - 1;
  focusedPanel.set(panels[prevIdx]);
}

// Direct navigation
export function focusPanel(panelId: string | null) {
  focusedPanel.set(panelId);
}
```

#### Navigation Utilities

```typescript
// Get panels for current page
function getPanelIdsForPage(pageKey: PageKey): string[] {
  return PANEL_ORDERS[pageKey] || [];
}

// Reset state (important for page transitions)
export function resetCarouselState() {
  focusedPanel.set(null);
}
```

## PanelCarousel Component

### Swiper.js Integration

The carousel uses Swiper.js for smooth animations and touch support:

```svelte
<script>
  import { Swiper, SwiperSlide } from 'swiper/svelte';
  import 'swiper/css';

  // Swiper configuration for optimal UX
  const swiperConfig = {
    spaceBetween: 20,
    slidesPerView: 1,
    speed: 600,
    allowTouchMove: true,
    preventClicks: false, // Allow button clicks
    preventClicksPropagation: false,
    slideToClickedSlide: false, // Prevent conflicts
    touchStartPreventDefault: false,
    keyboard: {
      enabled: true,
      onlyInViewport: true,
    },
    a11y: {
      enabled: true,
      nextSlideMessage: 'Next panel',
      prevSlideMessage: 'Previous panel',
      slideRole: 'tabpanel',
    },
  };
</script>
```

### State Synchronization

The critical challenge is keeping Swiper state synchronized with Svelte stores:

```typescript
// Subscribe to focus changes → update Swiper
unsub = focusedPanel.subscribe(($focused) => {
  if (!$focused || !initialized || !swiperEl?.swiper) return;

  const panelIndex = currentPanels.findIndex((panel) => panel.id === $focused);

  if (panelIndex !== -1 && swiperEl.swiper.activeIndex !== panelIndex) {
    // Update Swiper to match store state
    swiperEl.swiper.slideTo(panelIndex, 600);

    // Trigger plot resize after slide animation
    setTimeout(() => {
      resizePlotsInActiveSlide();
    }, 650); // Slightly after animation completes
  }
});

// Handle user swipe → update store
function handleSlideChange() {
  if (!swiperEl?.swiper || currentPanels.length === 0) return;

  const newIndex = swiperEl.swiper.activeIndex;
  const panelId = currentPanels[newIndex]?.id;

  if (panelId && panelId !== get(focusedPanel)) {
    focusedPanel.set(panelId); // This triggers all other updates
  }
}
```

### Plot Resizing Strategy

One of the most complex aspects is ensuring Plotly.js plots resize correctly:

#### Multi-Level Resize Strategy

```typescript
// 1. Fast resize: Only active slide
function resizePlotsInActiveSlide() {
  const activeSlide = swiperEl?.querySelector('.swiper-slide-active');
  if (!activeSlide) return;

  const plots = activeSlide.querySelectorAll('.js-plotly-plot');
  plots.forEach((plot) => {
    if (plot._plot) {
      Plotly.Plots.resize(plot);
    }
  });
}

// 2. Thorough resize: All slides (backup)
function resizeAllPlotsInCarousel() {
  if (!swiperEl) return;

  const allPlots = swiperEl.querySelectorAll('.swiper-slide .js-plotly-plot');
  allPlots.forEach((plot) => {
    if (plot._plot) {
      setTimeout(() => Plotly.Plots.resize(plot), 100);
    }
  });
}

// 3. Integration with global plot manager
async function triggerCarouselPlotResize() {
  const { triggerPlotResize } = await import('$lib/stores');
  await triggerPlotResize();
}
```

#### Timing Coordination

```typescript
// Resize after slide transitions
function handleSlideChange() {
  // Update store immediately
  updateFocusedPanel();

  // Resize plots after animation
  setTimeout(resizePlotsInActiveSlide, 300);
  setTimeout(resizeAllPlotsInCarousel, 500); // Backup
}

// Window resize handling
const handleResize = debounce(() => {
  requestAnimationFrame(async () => {
    await triggerCarouselPlotResize();
  });
}, 150);
```

### Keyboard Navigation

Comprehensive keyboard support with accessibility:

```typescript
let isNavigating = false; // Prevent rapid navigation

keyboardListener = (e: KeyboardEvent) => {
  // Only handle if carousel is active and no interactive element focused
  if (!initialized || !swiperEl?.swiper || isNavigating) return;

  // Skip if input field is focused
  if (
    e.target instanceof HTMLInputElement ||
    e.target instanceof HTMLTextAreaElement ||
    isInteractiveElement(e.target as HTMLElement)
  ) {
    return;
  }

  if (e.key === 'ArrowRight') {
    e.preventDefault();
    isNavigating = true;
    nextCarouselPanel($currentPage);

    // Reset navigation flag after animation
    setTimeout(() => {
      isNavigating = false;
    }, 650);
  } else if (e.key === 'ArrowLeft') {
    e.preventDefault();
    isNavigating = true;
    prevCarouselPanel($currentPage);

    setTimeout(() => {
      isNavigating = false;
    }, 650);
  }
};

// Interactive element detection
function isInteractiveElement(element: HTMLElement): boolean {
  const interactiveTags = ['BUTTON', 'TEXTAREA', 'INPUT', 'SELECT', 'A'];
  const interactiveClasses = ['send', 'export-btn', 'action-btn', 'cta-button'];
  const interactiveAttributes = ['contenteditable', 'tabindex'];

  return (
    interactiveTags.includes(element.tagName) ||
    interactiveClasses.some((cls) => element.classList.contains(cls)) ||
    interactiveAttributes.some((attr) => element.hasAttribute(attr)) ||
    element.closest('button, a, input, textarea, select') !== null
  );
}
```

## Common Issues and Solutions

### 1. Plot Rendering Issues

**Problem**: Plots not displaying correctly in carousel slides

**Causes**:

- Plot rendered when slide not visible
- Container dimensions not available
- Timing issues with Swiper initialization

**Solutions**:

```typescript
// Wait for slide to be visible before creating plot
const observer = new IntersectionObserver((entries) => {
  entries.forEach((entry) => {
    if (entry.isIntersecting) {
      const plotContainer = entry.target.querySelector('[data-plot-container]');
      if (plotContainer && !plotContainer._plot) {
        initializePlot(plotContainer);
      }
    }
  });
});

// Observe all slides
swiperEl.querySelectorAll('.swiper-slide').forEach((slide) => {
  observer.observe(slide);
});
```

### 2. Focus State Desynchronization

**Problem**: Swiper position doesn't match focused panel

**Causes**:

- Race conditions between user interaction and store updates
- Multiple components updating focus simultaneously
- Page navigation not clearing state properly

**Solutions**:

```typescript
// Single update path
export function setCarouselMode(pageKey: PageKey, mode: boolean) {
  if (!mode) {
    focusedPanel.set(null); // Only update the store
    return;
  }

  const firstPanelId = getFirstPanelIdForPage(pageKey);
  if (firstPanelId) {
    focusedPanel.set(firstPanelId); // Let everything else derive
  }
}

// Page transition cleanup
beforeUpdate(() => {
  if ($currentPage !== previousPage) {
    resetCarouselState();
    previousPage = $currentPage;
  }
});
```

### 3. Interactive Elements Not Working

**Problem**: Buttons/inputs in carousel slides don't respond to clicks

**Causes**:

- Swiper preventing click events
- CSS pointer-events disabled
- Z-index issues

**Solutions**:

```css
/* Ensure interactions work */
.swiper-slide button,
.swiper-slide textarea,
.swiper-slide input,
.swiper-slide select {
  pointer-events: auto !important;
  position: relative;
  z-index: 10;
  cursor: pointer;
}

/* Container styles */
.swiper-slide {
  overflow: visible;
}

.slide-content {
  height: 100%;
  display: flex;
  flex-direction: column;
}
```

```typescript
// Swiper configuration
Object.assign(swiperEl, {
  preventClicks: false,
  preventClicksPropagation: false,
  slideToClickedSlide: false,
  touchStartPreventDefault: false,
});
```

### 4. Memory Leaks

**Problem**: Application slows down over time due to accumulating event listeners

**Causes**:

- Event listeners not cleaned up
- Plotly instances not properly destroyed
- Store subscriptions not unsubscribed

**Solutions**:

```typescript
onDestroy(() => {
  // Clean up ALL listeners
  if (keyboardListener) {
    window.removeEventListener('keydown', keyboardListener);
  }

  window.removeEventListener('resize', handleResize);

  if (swiperEl) {
    swiperEl.removeEventListener('slidechange', handleSlideChange);
  }

  // Clean up store subscriptions
  unsub?.();

  // Clean up intersection observer
  if (observer) {
    observer.disconnect();
  }

  // Clean up plots
  const plots = swiperEl?.querySelectorAll('.js-plotly-plot');
  plots?.forEach((plot) => {
    if (plot._plot) {
      Plotly.purge(plot);
    }
  });
});
```

## Performance Optimizations

### Lazy Loading Strategy

Plot management and Swiper are loaded only when needed:

```typescript
// Lazy load plot utilities
async function initializePlotSupport() {
  const { triggerPlotResize } = await import('$lib/stores');
  return triggerPlotResize;
}

// Lazy load Swiper modules
import('swiper/modules').then(({ Navigation, Keyboard, A11y }) => {
  // Configure additional modules
});
```

### Efficient Updates

```typescript
// Debounced resize handling
const handleResize = debounce(() => {
  requestAnimationFrame(async () => {
    await triggerCarouselPlotResize();
  });
}, 150);

// Throttled scroll handling
const handleScroll = throttle(() => {
  updateVisiblePanels();
}, 50);
```

### Memory Management

```typescript
// Clean component state
function resetComponent() {
  initialized = false;
  isNavigating = false;
  currentPanels = [];

  // Clear timers
  if (resizeTimeout) {
    clearTimeout(resizeTimeout);
    resizeTimeout = null;
  }
}
```

## Debugging Guide

### State Debugging

Add these reactive statements for development:

```typescript
// Monitor carousel state
$: if (import.meta.env.DEV) {
  console.log('Carousel Debug:', {
    focusedPanel: $focusedPanel,
    carouselMode: $currentCarouselMode,
    currentIndex: $currentCarouselIndex,
    pageKey: $currentPage,
    panelCount: currentPanels.length,
  });
}

// Monitor Swiper state
$: if (import.meta.env.DEV && swiperEl?.swiper) {
  console.log('Swiper Debug:', {
    activeIndex: swiperEl.swiper.activeIndex,
    slidesLength: swiperEl.swiper.slides.length,
    initialized: swiperEl.swiper.initialized,
  });
}
```

### Plot Debugging

```typescript
// Monitor plot operations
function debugPlotResize() {
  const plots = document.querySelectorAll('.js-plotly-plot');
  console.log(`Found ${plots.length} plots for resize`);

  plots.forEach((plot, index) => {
    console.log(`Plot ${index}:`, {
      hasPlotlyInstance: !!plot._plot,
      isVisible: plot.offsetWidth > 0 && plot.offsetHeight > 0,
      dimensions: { width: plot.offsetWidth, height: plot.offsetHeight },
    });
  });
}
```

### Navigation Debugging

```typescript
// Monitor navigation events
export function nextCarouselPanel(pageKey: PageKey) {
  if (import.meta.env.DEV) {
    console.log(`Navigation: next panel on ${pageKey}`, {
      current: get(focusedPanel),
      available: getPanelIdsForPage(pageKey),
    });
  }

  // ... navigation logic
}
```

## Best Practices

### State Management

1. **Single Source of Truth**: Always use `focusedPanel` as the primary state
2. **Derived States**: Compute other values from `focusedPanel`
3. **Clean Transitions**: Reset state when changing pages
4. **Avoid Direct Manipulation**: Don't directly update Swiper state

### Performance

1. **Lazy Load**: Load heavy components only when needed
2. **Debounce**: Limit frequency of expensive operations
3. **Cleanup**: Always clean up resources in `onDestroy`
4. **Efficient Queries**: Use specific selectors for DOM queries

### Accessibility

1. **Keyboard Support**: Implement comprehensive keyboard navigation
2. **Focus Management**: Handle focus properly for interactive elements
3. **ARIA Labels**: Provide descriptive labels for navigation
4. **Screen Reader**: Test with screen readers

### Testing

1. **State Testing**: Test store integration thoroughly
2. **Navigation Testing**: Verify all navigation paths work
3. **Plot Testing**: Test plot rendering and resizing
4. **Cleanup Testing**: Verify proper resource cleanup

## Future Enhancements

- **Virtual Scrolling**: For carousels with many panels
- **Gesture Support**: Enhanced touch/swipe gestures
- **Transition Effects**: Custom transition animations
- **Panel Preloading**: Preload adjacent panels for smoother UX
- **Accessibility Improvements**: Enhanced screen reader support
- **Performance Monitoring**: Built-in performance metrics
