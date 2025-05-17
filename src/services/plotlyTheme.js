/**
 * Shared Plotly.js configuration for consistent appearance across visualization components
 */

// Default shared theme for all Plotly charts
export const defaultTheme = {
  paper_bgcolor: 'transparent',
  plot_bgcolor: 'transparent',
  font: {
    family: '"SF Pro", system-ui, -apple-system, sans-serif',
    color: 'rgba(60, 60, 67, 0.85)',
    size: 12
  },
  margin: {
    l: 60,
    r: 40,
    b: 50,
    t: 25,
    pad: 10
  },
  xaxis: {
    gridcolor: 'rgba(60, 60, 67, 0.08)',
    zerolinecolor: 'rgba(60, 60, 67, 0.15)',
    tickfont: {
      family: '"SF Pro", system-ui, -apple-system, sans-serif',
      size: 11,
      color: 'rgba(60, 60, 67, 0.6)'
    },
    linecolor: 'rgba(60, 60, 67, 0.12)'
  },
  yaxis: {
    gridcolor: 'rgba(60, 60, 67, 0.08)',
    zerolinecolor: 'rgba(60, 60, 67, 0.15)',
    tickfont: {
      family: '"SF Pro", system-ui, -apple-system, sans-serif',
      size: 11,
      color: 'rgba(60, 60, 67, 0.6)'
    },
    linecolor: 'rgba(60, 60, 67, 0.12)'
  },
  colorway: ['var(--accent)', 'var(--accent-secondary)', '#9d6bff', '#ff9f0a'],
  hoverlabel: {
    bgcolor: 'rgba(255, 255, 255, 0.95)',
    font: {
      family: '"SF Pro", system-ui, -apple-system, sans-serif',
      color: 'rgba(60, 60, 67, 0.85)',
      size: 12
    },
    bordercolor: 'var(--accent)'
  }
};

// Default configuration options
export const defaultConfig = {
  responsive: true,
  displayModeBar: false,
  useResizeHandler: true,
  toImageButtonOptions: {
    format: 'png',
    filename: 'spectrum_plot',
    height: 500,
    width: 700,
    scale: 2
  }
};

// Helper function for light mode plots
export function createLightModePlot(plotElement, data, layout) {
  const lightModeLayout = {
    ...defaultTheme,
    ...layout
  };
  
  // Add gradients to bar charts
  if (data[0]?.type === 'bar') {
    data[0].marker = {
      ...data[0].marker,
      color: 'var(--accent)',
      opacity: 0.9,
      line: {
        width: 1,
        color: 'rgba(0, 0, 0, 0.03)'
      }
    };
  }
  
  Plotly.newPlot(plotElement, data, lightModeLayout, defaultConfig);
}

// Helper function to merge user options with defaults
export function createLayout(customOptions = {}) {
  return { ...defaultTheme, ...customOptions };
} 