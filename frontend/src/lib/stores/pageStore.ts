/**
 * Page Management Store
 *
 * Handles page-specific state, navigation history, and SvelteKit page integration.
 * Separated for better modularity and to reduce coupling with carousel logic.
 */

import { writable, derived } from 'svelte/store';
import { page } from '$app/stores';
import type { Readable } from 'svelte/store';
import type { PageKey } from './types';

/**
 * SvelteKit-aware current page store with proper typing
 * Derives the current page from SvelteKit's page store
 */
export const currentPage: Readable<PageKey> = derived(page, ($page) => {
  if ($page.url.pathname === '/') return 'home';
  if ($page.url.pathname === '/spectral-simulation') return 'spectral-simulation';
  if ($page.url.pathname === '/how-it-works') return 'how-it-works';
  if ($page.url.pathname === '/about') return 'about';
  if ($page.url.pathname === '/chat-with-spectrum') return 'chat-with-spectrum';
  return 'home';
});

/**
 * Global app state stores
 */
export const isLoading = writable<boolean>(false);
export const appError = writable<string | null>(null);

/**
 * Navigation state management
 */
export const canNavigateBack = writable<boolean>(false);
export const navigationHistory = writable<string[]>([]);

/**
 * Update navigation history
 * @param path - The path to add to history
 */
export function addToNavigationHistory(path: string) {
  navigationHistory.update((history) => {
    const newHistory = [...history, path];
    // Keep history length reasonable
    return newHistory.slice(-10);
  });
  canNavigateBack.set(true);
}

/**
 * Clear navigation history
 */
export function clearNavigationHistory() {
  navigationHistory.set([]);
  canNavigateBack.set(false);
}

/**
 * Get the previous page from navigation history
 */
export function getPreviousPage(): string | null {
  let previousPage: string | null = null;
  navigationHistory.subscribe((history) => {
    previousPage = history.length > 1 ? history[history.length - 2] : null;
  })();
  return previousPage;
}
