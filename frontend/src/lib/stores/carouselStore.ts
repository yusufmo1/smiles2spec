/**
 * Carousel Navigation Store
 *
 * Handles panel focus, navigation, and carousel mode management.
 * Separated from main store for better modularity and tree-shaking.
 */

import { writable, derived, get } from 'svelte/store';
import type { PageKey } from './types';
import { simulationPanels, infoPanels, aboutPanels, chatPanels } from './panelStore';
import { stateBatcher } from '$lib/utils/performance';

/** ID of the panel that is currently in focus, or null */
export const focusedPanel = writable<string | null>(null);

/** Page-specific carousel indices for backward compatibility */
export const carouselIndices = writable<Record<PageKey, number>>({
  home: 0,
  'spectral-simulation': 0,
  'how-it-works': 0,
  about: 0,
  'chat-with-spectrum': 0,
});

/** Whether panel props are ready for rendering */
export const panelPropsReady = writable<boolean>(false);

/** Whether carousel data is synced */
export const carouselDataSync = writable<boolean>(true);

/**
 * Derive carousel mode from focusedPanel
 * Returns true when any panel is focused (carousel mode active)
 */
export const currentCarouselMode = derived(focusedPanel, ($focused) => $focused !== null);

/**
 * Derive carousel mode per page
 * Returns object with carousel state for each page
 */
export const carouselModes = derived([focusedPanel], ([$focused]) => ({
  home: $focused !== null,
  'spectral-simulation': $focused !== null,
  'how-it-works': $focused !== null,
  about: $focused !== null,
  'chat-with-spectrum': $focused !== null,
}));

/**
 * Derive current carousel index from focused panel
 * Finds the index of the focused panel within the current page's panel list
 */
export const currentCarouselIndex = derived(
  [focusedPanel, simulationPanels, infoPanels, aboutPanels, chatPanels],
  ([$focused, $sim, $info, $about, $chat]) => {
    if (!$focused) return 0;

    // Check all panel lists to find the focused panel
    const allPanels = [
      { panels: [], page: 'home' }, // Landing page has no panels
      { panels: $sim, page: 'spectral-simulation' },
      { panels: $info, page: 'how-it-works' },
      { panels: $about, page: 'about' },
      { panels: $chat, page: 'chat-with-spectrum' },
    ];

    for (const { panels } of allPanels) {
      const index = panels.findIndex((panel) => panel.id === $focused);
      if (index !== -1) return index;
    }

    return 0;
  }
);

/**
 * Panel ordering configuration for each page
 * Used for navigation and index calculations
 */
const PANEL_ORDERS = {
  home: [], // Landing page has no panels
  'spectral-simulation': ['spectrum', 'fragmentation', 'structure', 'peaks', 'console', 'export'],
  'how-it-works': [
    'system-overview',
    'workflow',
    'ml-engine',
    'fragmentation-science',
    'visualization-features',
    'accuracy-validation',
    'advanced-features',
    'technical-requirements',
  ],
  about: [
    'developer-hero',
    'developer-journey',
    'project-metrics',
    'tech-stack',
    'features-showcase',
    'open-source',
    'contact-connect',
    'acknowledgments',
  ],
  'chat-with-spectrum': ['chat-conversation', 'chat-spectrum', 'compound-stats'],
} as const;

/**
 * Get the first panel ID for a given page
 */
function getFirstPanelIdForPage(pageKey: PageKey): string {
  const panelOrder = PANEL_ORDERS[pageKey];
  return panelOrder?.[0] || 'spectrum';
}

/**
 * Get panel ID by index for a specific page
 */
function getPanelIdByIndex(pageKey: PageKey, index: number): string | null {
  const panelOrder = PANEL_ORDERS[pageKey];
  return panelOrder?.[index] || null;
}

/**
 * Get all panel IDs for a specific page
 */
function getPanelIdsForPage(pageKey: PageKey): readonly string[] {
  return PANEL_ORDERS[pageKey] || [];
}

/**
 * Set carousel mode for a specific page
 * @param pageKey - The page to set carousel mode for
 * @param mode - Whether to enable carousel mode
 */
export function setCarouselMode(pageKey: PageKey, mode: boolean) {
  if (mode) {
    const firstPanelId = getFirstPanelIdForPage(pageKey);
    focusedPanel.set(firstPanelId);
  } else {
    focusedPanel.set(null);
  }
}

/**
 * Set carousel index for a specific page
 * @param pageKey - The page to set index for
 * @param index - The panel index to focus
 */
export function setCarouselIndex(pageKey: PageKey, index: number) {
  const panelId = getPanelIdByIndex(pageKey, index);
  if (panelId) {
    focusedPanel.set(panelId);
  }
}

/**
 * Navigate to next panel in carousel
 * @param pageKey - Current page
 * @param maxIndex - Maximum index (unused, kept for compatibility)
 */
export function nextCarouselPanel(pageKey: PageKey, maxIndex?: number) {
  const panels = getPanelIdsForPage(pageKey);
  const currentPanelId = get(focusedPanel);
  const currentIdx = panels.findIndex((id) => id === currentPanelId);
  const nextIdx = (currentIdx + 1) % panels.length;
  focusedPanel.set(panels[nextIdx]);
}

/**
 * Navigate to previous panel in carousel
 * @param pageKey - Current page
 * @param maxIndex - Maximum index (unused, kept for compatibility)
 */
export function prevCarouselPanel(pageKey: PageKey, maxIndex?: number) {
  const panels = getPanelIdsForPage(pageKey);
  const currentPanelId = get(focusedPanel);
  const currentIdx = panels.findIndex((id) => id === currentPanelId);
  const prevIdx = currentIdx === 0 ? panels.length - 1 : currentIdx - 1;
  focusedPanel.set(panels[prevIdx]);
}

/**
 * Reset carousel state (exit carousel mode)
 * @param pageKey - Page to reset (unused, kept for compatibility)
 */
export function resetCarouselState(pageKey?: PageKey) {
  focusedPanel.set(null);
}

/**
 * Derived store for active panel data (backward compatibility)
 */
export const activePanelData = derived(
  [simulationPanels, carouselIndices],
  ([$panels, $indices]) => $panels[$indices.home] || null
);
