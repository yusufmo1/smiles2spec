# Store System - Modular Architecture

The SMILES2SPEC frontend uses a modular store architecture designed for performance, maintainability, and tree-shaking optimisation.

## Architecture Overview

```
stores/
├── index.ts            # Main exports with lazy loading utilities
├── appState.ts         # Global application state
├── carouselStore.ts    # Panel navigation and carousel
├── pageStore.ts        # Page routing and navigation
├── panelStore.ts       # Panel definitions and management
├── plotEffects.ts      # Lazy-loaded plot management
└── types.ts           # TypeScript type definitions
```

## Key Design Principles

### 🚀 **Performance First**

- **Lazy Loading**: Plot management code loads only when needed
- **Tree Shaking**: Unused store logic excluded from bundles
- **Code Splitting**: Heavy functionality in separate chunks

### 🧩 **Modular Design**

- **Single Responsibility**: Each store handles one concern
- **Clear Boundaries**: Well-defined interfaces between modules
- **Easy Testing**: Isolated modules simplify unit testing

### 📦 **Bundle Optimization**

- Plot effects separated into `plotEffects.ts` (lazy-loaded)
- Main bundle excludes heavy plot management code
- Better initial load performance

## Store Modules

### `index.ts` - Smart Entry Point

Main export file with intelligent lazy loading:

```typescript
// Direct imports (lightweight)
export { appState, currentPredictionData } from './appState';
export { focusedPanel, setCarouselMode } from './carouselStore';
export { currentPage, isLoading } from './pageStore';

// Lazy loading utilities (heavy functionality)
export async function initializePlots() {
  const { initializePlotSystem } = await import('./plotEffects');
  return initializePlotSystem();
}

export async function triggerPlotResize() {
  const { triggerPlotResize: trigger } = await import('./plotEffects');
  return trigger();
}
```

### `appState.ts` - Global State

Manages application-wide state including predictions, bulk processing, and console logging:

```typescript
import { appState } from '$lib/stores/appState';

// Core functionality
appState.setPredictionData(result);
appState.setBulkList(['CCO', 'CCC', 'CC(=O)O']);
appState.addConsoleEntry({
  type: 'success',
  message: 'Prediction completed',
});
```

**Features:**

- Prediction data management
- Bulk SMILES processing
- Console logging system
- Loading and error states

### `carouselStore.ts` - Panel Navigation

Handles panel focus, navigation, and carousel mode:

```typescript
import { focusedPanel, setCarouselMode, nextCarouselPanel } from '$lib/stores/carouselStore';

// Panel navigation
setCarouselMode('home', true); // Enter carousel mode
nextCarouselPanel('home'); // Navigate to next panel
focusedPanel.set('spectrum'); // Focus specific panel
```

**Features:**

- Panel focus management
- Carousel navigation controls
- Keyboard navigation support
- Page-specific panel ordering

### `pageStore.ts` - Page Management

Manages page routing, navigation history, and page-specific state:

```typescript
import { currentPage, isLoading, addToNavigationHistory } from '$lib/stores/pageStore';

// Page state
$: pageKey = $currentPage; // 'home' | 'about' | etc.
isLoading.set(true);
addToNavigationHistory('/chat-with-spectrum');
```

**Features:**

- SvelteKit page integration
- Navigation history tracking
- Global loading states
- Page-specific configuration

### `panelStore.ts` - Panel Definitions

Manages panel definitions, props, and categorization:

```typescript
import { panelStore, simulationPanels, infoPanels } from '$lib/stores/panelStore';

// Panel management
panelStore.updatePanelProps('spectrum', { data: newData });
$: currentPanels = $simulationPanels;
```

**Features:**

- Panel registration system
- Dynamic prop updates
- Category-based organization
- Component mapping

### `plotEffects.ts` - Lazy Plot Management

Heavy plot functionality loaded on-demand:

```typescript
// This module is lazy-loaded automatically
import { initializePlots } from '$lib/stores';

onMount(async () => {
  // Plot management loads only when needed
  const cleanup = await initializePlots();
  return cleanup;
});
```

**Features:**

- Lazy loading for performance
- Automatic plot resizing
- Carousel plot synchronization
- Memory cleanup

### `types.ts` - Type Definitions

Shared TypeScript interfaces and types:

```typescript
export type PageKey = 'home' | 'how-it-works' | 'about' | 'chat-with-spectrum';

export interface AppState {
  focusedPanel: string | null;
  carouselIndices: Record<string, number>;
  // ... more types
}
```

## Usage Patterns

### Basic Store Usage

```typescript
// Import from main index for convenience
import { appState, focusedPanel, currentPage } from '$lib/stores';

// Use reactive statements
$: predictionData = $appState.currentPredictionData;
$: isCarouselMode = $focusedPanel !== null;
$: pageKey = $currentPage;
```

### Performance-Optimized Usage

```typescript
// Import directly from modules for better tree-shaking
import { appState } from '$lib/stores/appState';
import { focusedPanel } from '$lib/stores/carouselStore';
import { currentPage } from '$lib/stores/pageStore';
```

### Lazy Plot Initialization

```typescript
import { onMount } from 'svelte';
import { initializePlots } from '$lib/stores';

let plotCleanup: (() => void) | undefined;

onMount(async () => {
  // Initialize plot system only when component mounts
  plotCleanup = await initializePlots();
});

onDestroy(() => {
  // Clean up plot effects
  plotCleanup?.();
});
```

## Migration Guide

### From Old Architecture

**Before (Monolithic):**

```typescript
import {
  focusedPanel,
  carouselIndex,
  isCarouselMode, // No longer needed
} from '$lib/stores';
```

**After (Modular):**

```typescript
import {
  focusedPanel,
  currentCarouselIndex, // Derived from focusedPanel
  currentCarouselMode, // Derived from focusedPanel
} from '$lib/stores';
```

## Performance Benefits

### Bundle Analysis

**Before Optimization:**

- Single 200+ line store file
- All plot management in main bundle
- Poor tree-shaking

**After Optimization:**

- Plot effects: Separate lazy-loaded chunk (`~1KB`)
- Main bundle: Reduced by excluding heavy plot code
- Better tree-shaking for unused functionality

### Lazy Loading

Plot management functionality is loaded only when:

- Component requires plot functionality
- `initializePlots()` is called
- Plot-heavy pages are accessed

## Best Practices

### Component Integration

```typescript
<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { appState, initializePlots } from '$lib/stores';

  let plotCleanup: (() => void) | undefined;

  // Reactive data
  $: predictionData = $appState.currentPredictionData;

  onMount(async () => {
    // Initialize heavy functionality only when needed
    if (needsPlots) {
      plotCleanup = await initializePlots();
    }
  });

  onDestroy(() => {
    plotCleanup?.();
  });
</script>
```

### Store Testing

```typescript
import { get } from 'svelte/store';
import { appState } from '$lib/stores/appState';

test('manages prediction data', () => {
  const mockData = { smiles: 'CCO', spectrum: {...} };
  appState.setPredictionData(mockData);

  const state = get(appState);
  expect(state.currentSmiles).toBe('CCO');
});
```

## Future Enhancements

- **Service Workers**: Cache store state for offline usage
- **Persistence**: Save state to localStorage
- **Undo/Redo**: Implement state history
- **Real-time Sync**: WebSocket integration for multi-user features
- **Performance Monitoring**: Track store performance metrics

## Debugging

### Store DevTools

```typescript
// Enable store debugging in development
if (import.meta.env.DEV) {
  appState.subscribe((state) => {
    console.log('AppState changed:', state);
  });
}
```

### Performance Monitoring

```typescript
// Monitor lazy loading performance
export async function initializePlots() {
  const start = performance.now();
  const { initializePlotSystem } = await import('./plotEffects');
  const end = performance.now();
  console.log(`Plot effects loaded in ${end - start}ms`);
  return initializePlotSystem();
}
```

The modular store architecture provides a solid foundation for scalable, performant state management while maintaining developer experience and code maintainability.
