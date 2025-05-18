import { writable } from 'svelte/store';

/** DOM element reference of the panel that is currently in focus, or null */
export const focusedPanel = writable(null);

/** index of the card that is currently centred in the carousel */
export const carouselIndex = writable(0);

/** Whether the app is currently in carousel mode */
export const isCarouselMode = writable(false); 