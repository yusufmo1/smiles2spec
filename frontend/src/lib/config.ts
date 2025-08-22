import { browser } from '$app/environment';

/**
 * Get the API URL based on environment
 */
export function getApiUrl(): string {
  if (browser) {
    const hostname = window.location.hostname;

    // Development environment
    if (hostname === 'localhost' || hostname === '127.0.0.1') {
      return '/api'; // Proxied via Vite
    }

    // Production environment
    return 'https://api.spectralsimulation.com';
  }

  // Server-side default (build time)
  return 'http://localhost:5050';
}

interface PlotlyConfig {
  displayModeBar: boolean;
  responsive: boolean;
  staticPlot: boolean;
  config: {
    displaylogo: boolean;
    modeBarButtonsToRemove: string[];
  };
}

interface ChatConfig {
  maxMessages: number;
  streamTimeout: number;
  retryAttempts: number;
}

interface UIConfig {
  animationDuration: number;
  debounceDelay: number;
  panelHeight: string;
  breakpoints: {
    mobile: number;
    tablet: number;
    desktop: number;
  };
}

interface Config {
  apiUrl: string;
  isDevelopment: boolean;
  plotly: PlotlyConfig;
  chat: ChatConfig;
  ui: UIConfig;
}

// Configuration object
export const config: Config = {
  // Update apiUrl to use the dynamic function
  get apiUrl() {
    return getApiUrl();
  },
  isDevelopment: false, // Will be set at runtime

  // Plotly configuration
  plotly: {
    displayModeBar: false,
    responsive: true,
    staticPlot: false,
    config: {
      displaylogo: false,
      modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d'],
    },
  },

  // Chat configuration
  chat: {
    maxMessages: 100,
    streamTimeout: 30000,
    retryAttempts: 3,
  },

  // UI configuration
  ui: {
    animationDuration: 300,
    debounceDelay: 500,
    panelHeight: '80vh',
    breakpoints: {
      mobile: 768,
      tablet: 1024,
      desktop: 1200,
    },
  },
};

export const getBaseUrl = (): string => {
  if (typeof window === 'undefined') {
    // Server-side rendering
    return '';
  }

  // Development environment
  if (window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1') {
    return '';
  }

  // Production environment - ensure consistent base path
  return '';
};

// Add environment detection helper
export const isProduction = (): boolean => {
  if (typeof window === 'undefined') return false;
  return (
    window.location.hostname !== 'localhost' &&
    window.location.hostname !== '127.0.0.1' &&
    window.location.port !== '8080' &&
    window.location.port !== '5173' &&
    window.location.port !== '5174'
  );
};

interface EnvironmentInfo {
  environment: string;
  hostname?: string;
  port?: string;
  href?: string;
  pathname?: string;
  apiUrl?: string;
  baseUrl?: string;
}

// Debug helper for production
export const getEnvironmentInfo = (): EnvironmentInfo => {
  if (typeof window === 'undefined') return { environment: 'server' };

  return {
    environment: isProduction() ? 'production' : 'development',
    hostname: window.location.hostname,
    port: window.location.port,
    href: window.location.href,
    pathname: window.location.pathname,
    apiUrl: getApiUrl(),
    baseUrl: getBaseUrl(),
  };
};
