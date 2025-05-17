import { writable } from 'svelte/store';

/** id of the panel that is currently in focus, or null */
export const focusedPanel = writable(null); 