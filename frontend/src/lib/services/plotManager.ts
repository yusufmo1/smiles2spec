import { browser } from '$app/environment';
import { debounce } from '$lib/utils/performance';

/**
 * PlotManager service to track and manage all plotly plots
 * Handles global operations like resize across all plots
 */
class PlotManager {
  private plots = new Set<HTMLElement>();
  
  // ✅ OPTIMIZED: Debounced resize to prevent multiple rapid calls
  private debouncedResizeAll = debounce(this.performResizeAll.bind(this), 100);
  private debouncedResizeCarousel = debounce(this.performCarouselResize.bind(this), 100);

  /**
   * Register a plot element to be managed
   */
  registerPlot(element: HTMLElement) {
    if (browser) {
      this.plots.add(element);
    }
  }

  /**
   * Unregister a plot element
   */
  unregisterPlot(element: HTMLElement) {
    this.plots.delete(element);
  }

  /**
   * Resize all registered plots (debounced)
   */
  async resizeAllPlots() {
    this.debouncedResizeAll();
  }

  /**
   * Resize all plots in carousel or overlay context (debounced)
   */
  async resizeCarouselPlots() {
    this.debouncedResizeCarousel();
  }

  /**
   * Force immediate resize of all plots (bypass debouncing)
   */
  async resizeAllPlotsImmediate() {
    return this.performResizeAll();
  }

  /**
   * Force immediate carousel resize (bypass debouncing)
   */
  async resizeCarouselPlotsImmediate() {
    return this.performCarouselResize();
  }

  /**
   * Internal: Perform actual resize of all registered plots
   */
  private async performResizeAll() {
    if (!browser) return;

    const { loadPlotly } = await import('./plotlyService');
    const Plotly = await loadPlotly();

    if (Plotly) {
      this.plots.forEach((plot) => {
        try {
          Plotly.Plots.resize(plot);
        } catch (error) {
          console.warn('Failed to resize plot:', error);
        }
      });
    }
  }

  /**
   * Internal: Perform actual carousel plot resize
   */
  private async performCarouselResize() {
    if (!browser) return;

    // Target plots specifically in carousel/overlay contexts
    const carouselPlots = document.querySelectorAll(
      '.swiper-slide .js-plotly-plot, .fullscreen-overlay .js-plotly-plot'
    );

    const { loadPlotly } = await import('./plotlyService');
    const Plotly = await loadPlotly();

    if (Plotly) {
      carouselPlots.forEach((plot) => {
        try {
          Plotly.Plots.resize(plot as HTMLElement);
        } catch (error) {
          console.warn('Failed to resize carousel plot:', error);
        }
      });
    }
  }
}

// Export a singleton instance
export const plotManager = new PlotManager();
