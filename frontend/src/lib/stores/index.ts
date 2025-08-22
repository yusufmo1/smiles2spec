/**
 * Modular Store Index
 *
 * Main entry point for the store system with lazy loading and better tree-shaking.
 * This file exports the most commonly used stores directly and provides lazy
 * imports for heavy functionality like plot management.
 */

// ===== CORE EXPORTS =====

// Page management (lightweight)
export {
  currentPage,
  isLoading,
  appError,
  canNavigateBack,
  navigationHistory,
  addToNavigationHistory,
  clearNavigationHistory,
  getPreviousPage,
} from './pageStore';

// Carousel and navigation (core functionality)
export {
  focusedPanel,
  carouselIndices,
  panelPropsReady,
  carouselDataSync,
  currentCarouselMode,
  carouselModes,
  currentCarouselIndex,
  activePanelData,
  setCarouselMode,
  setCarouselIndex,
  nextCarouselPanel,
  prevCarouselPanel,
  resetCarouselState,
} from './carouselStore';

// Panel store (already modular)
export {
  panelStore,
  simulationPanels,
  infoPanels,
  aboutPanels,
  chatPanels,
  type Panel,
  type PanelStoreState,
} from './panelStore';

// App state (already modular)
export {
  appState,
  currentPredictionData,
  consoleText,
  type GlobalAppState,
  type ConsoleEntry,
} from './appState';

// Types
export type { PageKey } from './types';

// ===== LAZY LOADING UTILITIES =====

/**
 * Lazy-load plot effects when needed
 *
 * Call this function when you need plot management functionality.
 * Returns a cleanup function to stop plot effects.
 *
 * @example
 * ```typescript
 * import { initializePlots } from '$lib/stores';
 *
 * onMount(async () => {
 *   const cleanup = await initializePlots();
 *   return cleanup; // Cleanup on unmount
 * });
 * ```
 */
export async function initializePlots() {
  const { initializePlotSystem } = await import('./plotEffects');
  return initializePlotSystem();
}

/**
 * Manually trigger plot resize (lazy-loaded)
 *
 * Use this when you need to resize plots without changing carousel state.
 *
 * @example
 * ```typescript
 * import { triggerPlotResize } from '$lib/stores';
 *
 * function handleResize() {
 *   triggerPlotResize();
 * }
 * ```
 */
export async function triggerPlotResize() {
  const { triggerPlotResize: trigger } = await import('./plotEffects');
  return trigger();
}

// ===== COMPLETE EXPORTS (for convenience) =====

/**
 * Re-export all stores for convenience
 * Import from specific modules when you need better tree-shaking
 */
export * from './carouselStore';
export * from './pageStore';
export * from './panelStore';
export * from './appState';
export * from './types';
