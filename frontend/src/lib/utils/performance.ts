/**
 * Performance Utilities
 * 
 * Optimised utilities for debouncing, batching, and performance monitoring
 */

/**
 * Creates a debounced version of the provided function
 * @param func Function to debounce
 * @param wait Delay in milliseconds
 * @param immediate Execute on leading edge instead of trailing
 */
export function debounce<T extends (...args: any[]) => any>(
  func: T,
  wait: number,
  immediate = false
): (...args: Parameters<T>) => void {
  let timeout: ReturnType<typeof setTimeout> | null = null;

  return function executedFunction(...args: Parameters<T>) {
    const later = () => {
      timeout = null;
      if (!immediate) func(...args);
    };

    const callNow = immediate && !timeout;
    
    if (timeout) clearTimeout(timeout);
    timeout = setTimeout(later, wait);
    
    if (callNow) func(...args);
  };
}

/**
 * Batches multiple operations in a single animation frame
 * @param operations Array of functions to execute
 */
export function batchOperations(operations: Array<() => void>): void {
  requestAnimationFrame(() => {
    operations.forEach(op => {
      try {
        op();
      } catch (error) {
        console.warn('Batched operation failed:', error);
      }
    });
  });
}

/**
 * Performance-optimised state updater
 * Batches multiple store updates in single frame
 */
export class StateBatcher {
  private pendingUpdates = new Set<() => void>();
  private frameRequested = false;

  /**
   * Queue a state update to be executed in next frame
   */
  queue(updateFn: () => void): void {
    this.pendingUpdates.add(updateFn);
    
    if (!this.frameRequested) {
      this.frameRequested = true;
      requestAnimationFrame(() => {
        this.flush();
      });
    }
  }

  /**
   * Execute all pending updates immediately
   */
  flush(): void {
    this.pendingUpdates.forEach(updateFn => {
      try {
        updateFn();
      } catch (error) {
        console.warn('State update failed:', error);
      }
    });
    
    this.pendingUpdates.clear();
    this.frameRequested = false;
  }
}

// Export singleton state batcher
export const stateBatcher = new StateBatcher();

/**
 * Development-only logging
 * Automatically disabled in production builds
 */
export const devLog = {
  log: (...args: any[]) => {
    if (import.meta.env.DEV) {
      console.log(...args);
    }
  },
  warn: (...args: any[]) => {
    if (import.meta.env.DEV) {
      console.warn(...args);
    }
  },
  error: (...args: any[]) => {
    // Always log errors, even in production
    console.error(...args);
  }
};