<script lang="ts">
  import { createEventDispatcher, onMount, onDestroy } from 'svelte';
  import { browser } from '$app/environment';
  import {
    focusedPanel,
    setCarouselIndex,
    currentCarouselIndex,
    simulationPanels,
    infoPanels,
    aboutPanels,
    chatPanels,
    triggerPlotResize,
  } from '$lib/stores/index';
  import { fade, scale } from 'svelte/transition';
  import { cubicOut } from 'svelte/easing';
  import PanelCarousel from './PanelCarousel.svelte';
  import type { Panel as PanelType } from '$lib/stores/panelStore';
  import type { PageKey } from '$lib/stores/index';
  import { handlePanelKeyNavigation, FocusManager } from '../utils/navigation';

  export let currentPage: string = 'home';

  const dispatch = createEventDispatcher<{
    close: void;
  }>();

  // ✅ SIMPLIFIED: Get current page carousel state
  $: pageKey = (currentPage as PageKey) || 'home';
  $: overlayVisible = $focusedPanel !== null;
  $: carouselIndex = $currentCarouselIndex;

  // Use reactive derived stores instead of getPanelsForPage
  $: panelsToShow =
    currentPage === 'home' || currentPage === 'spectral-simulation'
      ? $simulationPanels
      : currentPage === 'how-it-works'
        ? $infoPanels
        : currentPage === 'about'
          ? $aboutPanels
          : currentPage === 'chat-with-spectrum'
            ? $chatPanels
            : $simulationPanels;

  $: maxIndex = panelsToShow.length - 1;

  const focusManager = new FocusManager();

  // ✅ SIMPLIFIED: Only modify focusedPanel
  const close = () => {
    focusedPanel.set(null);
    dispatch('close');
  };

  // ✅ SIMPLIFIED: Use setCarouselIndex helper
  const next = () => {
    const newIndex = Math.min(carouselIndex + 1, maxIndex);
    setCarouselIndex(pageKey, newIndex);
  };

  const prev = () => {
    const newIndex = Math.max(carouselIndex - 1, 0);
    setCarouselIndex(pageKey, newIndex);
  };

  const goToSlide = (index: number) => {
    setCarouselIndex(pageKey, index);
  };

  let isNavigating = false;

  const handleKey = (e: KeyboardEvent) => {
    if (!$focusedPanel) return;
    if (isNavigating) return; // Prevent multiple rapid keydowns

    // Ensure this is a direct overlay keyboard event, not bubbled
    if (e.target !== document.body && !e.currentTarget) return;

    // Add more specific navigation key handling to prevent double-stepping
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight' || e.key === 'ArrowUp' || e.key === 'ArrowDown') {
      e.preventDefault();
      e.stopImmediatePropagation(); // Stop any other handlers from firing
      
      isNavigating = true;
      
      // Handle navigation with debounced execution
      requestAnimationFrame(() => {
        handlePanelKeyNavigation(
          e,
          carouselIndex,
          panelsToShow.length,
          (index) => {
            setCarouselIndex(pageKey, index);
          },
          close
        );
        
        // Reset navigation flag after animation completes
        setTimeout(() => {
          isNavigating = false;
        }, 650); // Slightly longer than carousel transition
      });
    } else if (e.key === 'Escape') {
      e.preventDefault();
      close();
    }
  };

  const handleDotKeydown = (e: KeyboardEvent, index: number) => {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      goToSlide(index);
    }
  };

  // UPDATED: Resize ALL plots in carousel, not just active slide
  function handleTransitionEnd() {
    setTimeout(async () => {
      // Resize ALL plots when overlay opens, not just active slide
      await triggerPlotResize();
    }, 100);
  }

  // ✅ SIMPLIFIED: Body scroll prevention
  $: if (browser && $focusedPanel) {
    // Store previous focus
    focusManager.saveFocus();

    // Prevent body scroll
    document.body.style.overflow = 'hidden';

    // Add keyboard listener with capture to handle before other handlers
    document.addEventListener('keydown', handleKey, { capture: true });

    // NEW: Add immediate resize when overlay opens
    // Give DOM time to render, then resize all plots
    setTimeout(async () => {
      await triggerPlotResize();
    }, 500); // Longer delay for overlay transition
  } else if (browser) {
    // Restore body scroll
    document.body.style.overflow = '';

    // Remove keyboard listener with same capture option
    document.removeEventListener('keydown', handleKey, { capture: true });

    // Restore previous focus
    focusManager.restoreFocus();
  }

  onDestroy(() => {
    if (browser) {
      document.body.style.overflow = '';
      document.removeEventListener('keydown', handleKey, { capture: true });
      // Reset navigation flag if component is destroyed during navigation
      isNavigating = false;
    }
  });
</script>

{#if $focusedPanel}
  <div
    class="fullscreen-overlay"
    role="dialog"
    tabindex="-1"
    on:click|self={close}
    on:keydown={(e) => e.key === 'Escape' && close()}
    aria-modal="true"
    transition:fade={{ duration: 400 }}
    on:introend={handleTransitionEnd}
  >
    <div class="overlay-title">
      <h1 class="display text-gradient prosto-one-regular">
        {#if currentPage === 'home' || currentPage === 'spectral-simulation'}
          SPECTRAL SIMULATION
        {:else if currentPage === 'how-it-works'}
          HOW IT WORKS
        {:else if currentPage === 'about'}
          ABOUT
        {:else if currentPage === 'chat-with-spectrum'}
          CHAT WITH SPECTRUM
        {/if}
      </h1>
    </div>

    <div
      class="overlay-content"
      transition:scale={{ duration: 400, easing: cubicOut, start: 0.85 }}
    >
      <PanelCarousel currentPanels={panelsToShow} hideNavigation={true} />

      <button
        class="arrow left"
        on:click={prev}
        disabled={carouselIndex === 0}
        aria-label="Previous panel">◀</button
      >
      <button
        class="arrow right"
        on:click={next}
        disabled={carouselIndex === maxIndex}
        aria-label="Next panel">▶</button
      >

      <button class="control-button close" on:click={close} aria-label="Close">✕</button>

      <div class="pagination-indicator" role="tablist" aria-label="Panel navigation">
        {#each panelsToShow as _, i}
          <button
            class="indicator-dot"
            class:active={i === carouselIndex}
            on:click={() => goToSlide(i)}
            on:keydown={(e) => handleDotKeydown(e, i)}
            role="tab"
            aria-selected={i === carouselIndex}
            aria-label={`Go to panel ${i + 1}`}
          ></button>
        {/each}
      </div>

      <div class="keyboard-shortcuts">
        <div class="shortcut"><kbd>Esc</kbd> Close</div>
        <div class="shortcut"><kbd>←</kbd><kbd>→</kbd> Navigate</div>
      </div>
    </div>

    <div class="overlay-watermark">Spectral Simulation Tool</div>
  </div>
{/if}

<style>
  .fullscreen-overlay {
    position: fixed;
    inset: 0; /* fills the viewport by itself */
    z-index: 10000;
    background: rgba(247, 247, 248, 0.95);
    backdrop-filter: blur(25px) saturate(180%);
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding: 0;
    /* Modern browsers: dynamic units that ignore browser chrome */
    width: 100dvw; /* fallback-safe: ignored if not supported */
    height: 100dvh;
    animation: reveal 0.5s cubic-bezier(0.16, 1, 0.3, 1) forwards;
    box-sizing: border-box;
    /* GPU acceleration for smooth overlay transitions */
    will-change: opacity, transform;
    transform: translate3d(0, 0, 0);
  }

  /* Fallback for older browsers that don't support dynamic viewport units */
  :global(.fullscreen-overlay) {
    width: 100%;
    height: 100%;
  }

  /* Enhanced vignette effect */
  .fullscreen-overlay::after {
    content: '';
    position: absolute;
    inset: 0;
    background: radial-gradient(
      circle at center,
      transparent 30%,
      rgba(0, 0, 0, 0.05) 60%,
      rgba(0, 0, 0, 0.2) 100%
    );
    pointer-events: none;
  }

  @keyframes reveal {
    from {
      opacity: 0;
      transform: translate3d(0, 40px, 0) scale(0.95);
    }
    to {
      opacity: 1;
      transform: translate3d(0, 0, 0) scale(1);
    }
  }

  .overlay-title {
    position: absolute;
    top: 2rem;
    text-align: center;
    z-index: 10001;
  }

  .overlay-title h1 {
    font:
      600 clamp(2.5rem, 4vw, 3.75rem) / 1 'Prosto One',
      sans-serif;
    margin: 0;
    letter-spacing: 0.05em;
  }

  /* ✅ FIXED: Use original structure */
  .overlay-content {
    position: relative;
    width: 100%; /* Use 100% instead of 100vw to respect parent boundaries */
    height: 80vh;
    margin-top: 5rem;
    display: flex;
    align-items: center;
    justify-content: center;
    box-sizing: border-box;
  }

  .arrow {
    position: absolute;
    top: 50%;
    transform: translate3d(0, -50%, 0);
    background: var(--accent-soft);
    border: none;
    border-radius: 50%;
    width: 60px;
    height: 60px;
    font-size: 1.2rem;
    color: var(--accent);
    cursor: pointer;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
    z-index: 10;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all var(--transition-smooth);
    /* GPU acceleration for smooth button interactions */
    will-change: transform, background-color;
  }

  .arrow:hover:not(:disabled) {
    transform: translate3d(0, -50%, 0) scale(1.1);
    background: var(--accent);
    color: white;
    box-shadow: 0 6px 20px rgba(120, 121, 255, 0.4);
  }

  .arrow:disabled {
    opacity: 0.3;
    cursor: not-allowed;
  }

  .arrow.left {
    left: 5%;
  }
  .arrow.right {
    right: 5%;
  }

  .control-button {
    position: absolute;
    top: -3rem;
    right: 2rem;
    width: 50px;
    height: 50px;
    border: none;
    border-radius: 50%;
    background: var(--accent-soft);
    color: var(--accent);
    cursor: pointer;
    font-size: 1.25rem;
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    transition: all var(--transition-smooth);
  }

  .control-button:hover {
    transform: scale(1.1);
    background: #ff453a;
    color: white;
  }

  .pagination-indicator {
    position: absolute;
    bottom: 2rem;
    left: 50%;
    transform: translateX(-50%);
    display: flex;
    gap: 8px;
    z-index: 10;
  }

  .indicator-dot {
    width: 12px;
    height: 12px;
    border-radius: 50%;
    background: var(--accent-soft);
    border: none;
    cursor: pointer;
    transition: all 0.3s ease;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .indicator-dot:hover {
    background: var(--accent);
    transform: scale(1.1);
  }

  .indicator-dot:focus {
    outline: 2px solid var(--accent);
    outline-offset: 2px;
  }

  .indicator-dot.active {
    background: var(--accent);
    transform: scale(1.25);
  }

  .keyboard-shortcuts {
    position: fixed;
    bottom: 1.5rem;
    left: 50%;
    transform: translateX(-50%);
    display: flex;
    gap: 1.5rem;
    color: white;
    font-size: 0.85rem;
    z-index: 10002;
    background: var(--accent);
    padding: 0.5rem 1.25rem;
    border-radius: 50px;
    box-shadow: 0 4px 12px rgba(120, 121, 255, 0.4);
  }

  .shortcut {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    font-weight: 500;
  }

  kbd {
    display: inline-block;
    padding: 0.2rem 0.5rem;
    background: rgba(255, 255, 255, 0.25);
    border: 1px solid rgba(255, 255, 255, 0.3);
    border-radius: 4px;
    box-shadow: 0 2px 0 rgba(0, 0, 0, 0.1);
    font-family: 'SF Mono', monospace;
    font-size: 0.8rem;
    color: white;
  }

  .overlay-watermark {
    position: absolute;
    bottom: 1rem;
    left: 1rem;
    font-size: 0.8rem;
    opacity: 0.4;
    color: var(--text-tertiary);
  }

  @media (max-width: 1280px) {
    .arrow.left {
      left: 5%;
    }
    .arrow.right {
      right: 5%;
    }
  }

  @media (max-width: 768px) {
    .overlay-content {
      height: 75vh;
    }

    .arrow.left {
      left: 5%;
    }
    .arrow.right {
      right: 5%;
    }

    .keyboard-shortcuts {
      flex-direction: column;
      align-items: center;
      gap: 0.75rem;
    }
  }
</style>
