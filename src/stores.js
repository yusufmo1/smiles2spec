import { writable } from 'svelte/store';

/** DOM element reference of the panel that is currently in focus, or null */
export const focusedPanel = writable(null);

/** index of the card that is currently centred in the carousel */
export const carouselIndex = writable(0); 