import { browser } from '$app/environment';

let PlotlyModule: any = null;
let loadingPromise: Promise<any> | null = null;

/**
 * Safely load Plotly.js with SSR compatibility and caching
 */
export async function loadPlotly(): Promise<any> {
  if (!browser) {
    console.warn('Plotly can only be loaded in the browser');
    return null;
  }

  // Return cached module if already loaded
  if (PlotlyModule) {
    return PlotlyModule;
  }

  // Return existing loading promise if already loading
  if (loadingPromise) {
    return loadingPromise;
  }

  // Start loading Plotly with retry logic and better error handling
  loadingPromise = import('plotly.js-dist-min')
    .then((module) => {
      PlotlyModule = module.default || module;
      console.log('Plotly.js loaded successfully');
      return PlotlyModule;
    })
    .catch((error) => {
      console.error('Failed to load Plotly.js:', error);
      loadingPromise = null; // Reset to allow retry
      return null; // Return null instead of throwing
    });

  return loadingPromise;
}

/**
 * Create a new plot with proper error handling
 */
export async function createPlot(
  element: HTMLElement,
  data: any,
  layout: any = {},
  config: any = {}
): Promise<void> {
  if (!browser || !element) {
    console.warn('createPlot: Browser or element not available');
    return;
  }

  try {
    const Plotly = await loadPlotly();
    if (!Plotly) {
      throw new Error('Plotly failed to load');
    }

    // Default config for consistent behavior
    const plotConfig = {
      responsive: true,
      displayModeBar: false,
      staticPlot: false,
      ...config,
    };

    // Default layout for SSG compatibility
    const plotLayout = {
      autosize: true,
      paper_bgcolor: 'transparent',
      plot_bgcolor: 'transparent',
      font: {
        family: 'SF Pro Display, -apple-system, sans-serif',
        size: 12,
        color: '#1d1d1f',
      },
      margin: { l: 60, r: 40, b: 50, t: 30 },
      ...layout,
    };

    await Plotly.newPlot(element, data, plotLayout, plotConfig);
    console.log('Plot created successfully');
  } catch (error) {
    console.error('Failed to create plot:', error);
    throw error;
  }
}

/**
 * Update an existing plot
 */
export async function updatePlot(element: HTMLElement, data: any, layout: any = {}): Promise<void> {
  if (!browser || !element) return;

  try {
    const Plotly = await loadPlotly();
    if (!Plotly) return;

    await Plotly.react(element, data, layout);
    console.log('Plot updated successfully');
  } catch (error) {
    console.error('Failed to update plot:', error);
  }
}

/**
 * Resize plot for responsive behavior
 */
export async function resizePlot(element: HTMLElement): Promise<void> {
  if (!browser || !element) return;

  try {
    const Plotly = await loadPlotly();
    if (!Plotly) return;

    await Plotly.Plots.resize(element);
  } catch (error) {
    console.error('Failed to resize plot:', error);
  }
}

/**
 * Download plot as image
 */
export async function downloadPlot(
  element: HTMLElement,
  filename: string = 'plot',
  format: 'png' | 'jpeg' | 'pdf' | 'svg' = 'png',
  width: number = 1200,
  height: number = 800
): Promise<void> {
  if (!browser || !element) return;

  try {
    const Plotly = await loadPlotly();
    if (!Plotly) return;

    const dataUrl = await Plotly.toImage(element, {
      format,
      width,
      height,
      scale: 2,
    });

    const link = document.createElement('a');
    link.href = dataUrl;
    link.download = `${filename}.${format}`;
    link.click();
  } catch (error) {
    console.error('Failed to download plot:', error);
    throw error;
  }
}
