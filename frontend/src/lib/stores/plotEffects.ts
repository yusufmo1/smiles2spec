/**
 * Plot Effects Module
 *
 * Lazy-loaded plot management side effects for better performance.
 * Only imported when plot functionality is needed, reducing bundle size
 * for pages that don't require plot management.
 */

import { browser } from '$app/environment';
import { focusedPanel } from './carouselStore';

/**
 * Initialize plot resize effects
 *
 * Sets up automatic plot resizing when carousel mode changes.
 * Uses lazy imports to avoid loading plot management code until needed.
 */
export function initializePlotEffects() {
  if (!browser) return;

  let resizeTimeout: ReturnType<typeof setTimeout>;

  // Subscribe to focusedPanel changes for debounced plot resizing
  const unsubscribe = focusedPanel.subscribe(async ($focused) => {
    if (!browser) return;

    // Clear any pending resize timeouts
    clearTimeout(resizeTimeout);

    if ($focused) {
      // Shorter delay + debounced resize prevents multiple calls
      resizeTimeout = setTimeout(async () => {
        try {
          // Dynamic import to avoid SSR issues and reduce initial bundle size
          const { plotManager } = await import('$lib/services/plotManager');
          // Debounced resize automatically handles rapid successive calls
          plotManager.resizeCarouselPlots();
        } catch (error) {
          console.warn('Plot manager not available:', error);
        }
      }, 200); // Reduced from 400ms to 200ms since resize is now debounced
    }
  });

  // Return cleanup function
  return unsubscribe;
}

/**
 * Manually trigger plot resize
 *
 * Useful for components that need to force a plot resize
 * without changing carousel state.
 */
export async function triggerPlotResize() {
  if (!browser) return;

  try {
    const { plotManager } = await import('$lib/services/plotManager');
    plotManager.resizeCarouselPlots();
  } catch (error) {
    console.warn('Plot manager not available:', error);
  }
}

/**
 * Initialize all plot-related effects
 *
 * Call this function when plot functionality is needed.
 * Returns cleanup function to stop all effects.
 */
export function initializePlotSystem() {
  const cleanupFunctions: Array<() => void> = [];

  // Initialize plot resize effects
  const plotEffectsCleanup = initializePlotEffects();
  if (plotEffectsCleanup) {
    cleanupFunctions.push(plotEffectsCleanup);
  }

  // Return combined cleanup function
  return () => {
    cleanupFunctions.forEach((cleanup) => cleanup());
  };
}
