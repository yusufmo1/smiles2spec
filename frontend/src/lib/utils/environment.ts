import { browser } from '$app/environment';

export interface EnvironmentInfo {
  isDevelopment: boolean;
  isProduction: boolean;
  apiUrl: string;
  hostname: string;
  port: string;
}

export function getEnvironmentInfo(): EnvironmentInfo {
  if (!browser) {
    return {
      isDevelopment: true,
      isProduction: false,
      apiUrl: 'http://localhost:5050',
      hostname: 'localhost',
      port: '5050',
    };
  }

  const hostname = window.location.hostname;
  const isDevelopment = hostname === 'localhost' || hostname === '127.0.0.1';
  const isProduction = !isDevelopment;

  const apiUrl = isDevelopment
    ? '/api' // Proxied in development
    : 'https://api.spectralsimulation.com'; // Direct in production

  return {
    isDevelopment,
    isProduction,
    apiUrl,
    hostname,
    port: window.location.port,
  };
}

export function logEnvironmentInfo(): void {
  if (browser) {
    const info = getEnvironmentInfo();
    console.log('🌍 Environment Info:', info);
  }
}
