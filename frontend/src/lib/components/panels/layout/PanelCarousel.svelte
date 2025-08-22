<script lang="ts">
  import { focusedPanel, triggerPlotResize } from '$lib/stores/index';
  import Panel from '../Panel.svelte';
  import { onMount, onDestroy } from 'svelte';
  import { browser } from '$app/environment';
  import type { Panel as PanelType } from '$lib/stores/panelStore';
  import type { PageKey } from '$lib/stores/index';
  import { DynamicIcon } from '$lib/components/icons';
  import { swiperService } from '$lib/services/swiperService';

  export let currentPanels: PanelType[] = [];
  export let hideNavigation: boolean = false;

  let swiperEl: any;
  let initialized = false;
  let unsub: (() => void) | undefined;
  let isCarouselInitialized = false;
  let keyboardListener: ((e: KeyboardEvent) => void) | undefined;
  let isNavigating = false;

  // Enhanced resize handling for Plotly charts
  function resizePlotsInActiveSlide() {
    const activeSlide = swiperEl?.querySelector('.swiper-slide-active');
    if (activeSlide) {
      const plots = activeSlide.querySelectorAll('.js-plotly-plot');
      plots.forEach((plot: any) => {
        if ((window as any).Plotly) {
          (window as any).Plotly.Plots.resize(plot);
        }
      });
    }
  }

  // NEW: Resize ALL plots, not just active slide
  function resizeAllPlotsInCarousel() {
    if (!swiperEl) return;

    // Find ALL slides, not just active one
    const allSlides = swiperEl.querySelectorAll('.swiper-slide');
    allSlides.forEach((slide: Element) => {
      const plots = slide.querySelectorAll('.js-plotly-plot');
      plots.forEach((plot: Element) => {
        if ((window as any).Plotly) {
          (window as any).Plotly.Plots.resize(plot);
        }
      });
    });
  }

  // ✅ UPDATE: Modified to set focusedPanel directly
  function handleSlideChange() {
    if (swiperEl?.swiper && currentPanels.length > 0) {
      const newIndex = swiperEl.swiper.activeIndex;
      const panelId = currentPanels[newIndex]?.id;

      if (panelId) {
        focusedPanel.set(panelId); // ← This updates everything else
      }

      // Resize Plotly charts in active slide after slide change
      setTimeout(() => {
        resizePlotsInActiveSlide();
      }, 300);
    }
  }

  // Set up coverflow effect options
  const coverflowEffect = {
    rotate: 0,
    stretch: 120,
    depth: 250,
    modifier: 1.05,
    slideShadows: false,
  };

  onMount(async () => {
    if (!browser || !swiperEl) return;

    try {
      // Use singleton service instead of heavy import
      await swiperService.ensureRegistered();

      // Configure swiper after component is mounted
      swiperEl.addEventListener('slidechange', handleSlideChange);

      // Use service configuration with coverflow effect
      swiperService.configureSwiper(swiperEl, {
        coverflowEffect: coverflowEffect,
      });

      await swiperEl.initialize();
      initialized = true;

      // Add keyboard navigation for inline carousel mode
      keyboardListener = (e: KeyboardEvent) => {
        if (!initialized || !swiperEl?.swiper || isNavigating) return;

        // Only handle if no modal/overlay is open and we're not in an input field
        // Also check if any overlay elements exist in DOM to prevent conflicts
        if (
          $focusedPanel ||
          document.querySelector('.fullscreen-overlay') ||
          e.target instanceof HTMLInputElement ||
          e.target instanceof HTMLTextAreaElement
        ) {
          return;
        }

        if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') {
          e.preventDefault();
          e.stopPropagation(); // Prevent event bubbling
          isNavigating = true;

          if (e.key === 'ArrowRight') {
            swiperEl.swiper.slideNext();
          } else if (e.key === 'ArrowLeft') {
            swiperEl.swiper.slidePrev();
          }

          // Reset navigation flag after animation
          setTimeout(() => {
            isNavigating = false;
          }, 600);
        }
      };

      window.addEventListener('keydown', keyboardListener);

      // Subscribe to focusedPanel changes instead of indices
      unsub = focusedPanel.subscribe(($focused) => {
        if (!$focused || !initialized || !swiperEl?.swiper) return;

        // Find index of focused panel
        const panelIndex = currentPanels.findIndex((panel) => panel.id === $focused);
        if (panelIndex !== -1 && swiperEl.swiper.activeIndex !== panelIndex) {
          swiperEl.swiper.slideTo(panelIndex, 600);
        }
      });

      // Initial resize of Plotly charts - now resize ALL plots
      setTimeout(() => {
        resizeAllPlotsInCarousel();
      }, 500);

      // Handle window resize events to resize all plots
      const handleResize = () => {
        requestAnimationFrame(async () => {
          await triggerPlotResize();
        });
      };

      window.addEventListener('resize', handleResize);

      onDestroy(() => {
        window.removeEventListener('resize', handleResize);
        swiperEl?.removeEventListener('slidechange', handleSlideChange);
        if (keyboardListener) {
          window.removeEventListener('keydown', keyboardListener);
        }
        unsub?.();
      });
    } catch (error) {
      console.warn('Failed to initialize Swiper:', error);
    }
  });

  // NEW: Detect when carousel becomes visible and resize all plots
  $: if (browser && swiperEl && initialized && !isCarouselInitialized) {
    isCarouselInitialized = true;
    // Delay to ensure DOM is ready
    setTimeout(async () => {
      resizeAllPlotsInCarousel();
      await triggerPlotResize();
    }, 300);
  }

  // Additional cleanup for keyboard listener
  onDestroy(() => {
    if (browser && keyboardListener) {
      window.removeEventListener('keydown', keyboardListener);
      isNavigating = false;
    }
  });
</script>

<div class="carousel-container">
  {#if browser}
    <swiper-container bind:this={swiperEl} init="false">
      {#each currentPanels as panel (panel.id)}
        <swiper-slide>
          <div class="slide-content">
            <Panel {...panel} clickable={true} isCarousel={true} />
          </div>
        </swiper-slide>
      {/each}
    </swiper-container>

    <!-- Navigation arrows for inline carousel -->
    {#if !hideNavigation && currentPanels.length > 1 && initialized}
      <button
        class="nav-arrow prev"
        on:click={() => swiperEl?.swiper?.slidePrev()}
        disabled={!swiperEl?.swiper || swiperEl.swiper.isBeginning}
        aria-label="Previous panel"
      >
        <DynamicIcon name="ArrowLeftIcon" size={20} />
      </button>
      <button
        class="nav-arrow next"
        on:click={() => swiperEl?.swiper?.slideNext()}
        disabled={!swiperEl?.swiper || swiperEl.swiper.isEnd}
        aria-label="Next panel"
      >
        <DynamicIcon name="ArrowRightIcon" size={20} />
      </button>
    {/if}

    <!-- Keyboard hints -->
    {#if !hideNavigation && currentPanels.length > 1 && !$focusedPanel}
      <div class="keyboard-hints">
        <div class="hint-content">
          <kbd>←</kbd> <kbd>→</kbd> to navigate
        </div>
      </div>
    {/if}
  {:else}
    <!-- SSR fallback -->
    <div class="ssr-fallback">
      {#each currentPanels as panel, index}
        {#if index === 0}
          <div class="slide-content">
            <Panel {...panel} clickable={true} isCarousel={true} />
          </div>
        {/if}
      {/each}
    </div>
  {/if}
</div>

<style>
  .carousel-container {
    width: 100%;
    height: 100%;
  }

  /* ✅ FIXED: Use original dimensions */
  .slide-content {
    width: 100%;
    height: 80vh; /* Fixed height like original */
    /* Ensure interactions can bubble up */
    position: relative;
    z-index: 1;
  }

  .ssr-fallback {
    width: 100%;
    height: 80vh;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  /* GPU-accelerated slide styling */
  :global(.swiper-slide) {
    width: 70vw;
    max-width: 1100px;
    flex-shrink: 0;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    /* Ensure clicks can pass through */
    pointer-events: auto;
    /* GPU acceleration for smooth transitions */
    will-change: transform;
    transform: translate3d(0, 0, 0) scale(0.85);
    transition: transform 0.35s cubic-bezier(0.4, 0, 0.2, 1);
  }
  :global(.swiper-slide-active) {
    transform: translate3d(0, 0, 0) scale(1);
  }
  :global(.swiper-slide-prev),
  :global(.swiper-slide-next) {
    transform: translate3d(0, 0, 0) scale(0.9);
  }

  /* ✅ FIXED: Responsive tweak */
  @media (max-width: 768px) {
    :global(.swiper-slide) {
      width: 90vw;
    }

    .slide-content {
      height: 75vh; /* Smaller on mobile */
    }
  }

  /* IMPORTANT: Ensure interactive elements work in carousel */
  :global(.slide-content) {
    pointer-events: auto;
  }

  /* Ensure all interactive elements can receive clicks */
  :global(.slide-content button),
  :global(.slide-content textarea),
  :global(.slide-content input),
  :global(.slide-content select),
  :global(.slide-content a),
  :global(.slide-content [role='button']),
  :global(.slide-content .clickable),
  :global(.slide-content .send),
  :global(.slide-content .export-btn),
  :global(.slide-content th[role='button']),
  :global(.slide-content .indicator-dot) {
    pointer-events: auto !important;
    position: relative;
    z-index: 10;
    cursor: pointer;
  }

  :global(.swiper-slide:not(.swiper-slide-active)) {
    filter: brightness(0.8) blur(1px) grayscale(0.15);
    transition: filter 0.25s cubic-bezier(0.4, 0, 0.2, 1);
    /* GPU acceleration for filter effects */
    will-change: filter;
  }

  /* Override swiper's click prevention */
  :global(.swiper-container) {
    touch-action: pan-y pinch-zoom;
  }

  /* Navigation arrows for inline carousel */
  .nav-arrow {
    position: absolute;
    top: 50%;
    transform: translateY(-50%);
    z-index: 15;
    background: var(--accent-soft);
    border: none;
    border-radius: 50%;
    width: 50px;
    height: 50px;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    transition: all var(--transition-smooth);
    color: var(--accent);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    backdrop-filter: blur(10px);
  }

  .nav-arrow:hover:not(:disabled) {
    background: var(--accent);
    color: white;
    transform: translateY(-50%) scale(1.1);
    box-shadow: 0 6px 20px rgba(120, 121, 255, 0.4);
  }

  .nav-arrow:disabled {
    opacity: 0.3;
    cursor: not-allowed;
  }

  .nav-arrow.prev {
    left: 10px;
  }

  .nav-arrow.next {
    right: 10px;
  }

  /* Keyboard hints */
  .keyboard-hints {
    position: absolute;
    bottom: 1rem;
    left: 50%;
    transform: translateX(-50%);
    z-index: 15;
    pointer-events: none;
  }

  .hint-content {
    background: var(--accent);
    color: white;
    padding: 0.5rem 1rem;
    border-radius: 50px;
    font-size: 0.85rem;
    font-weight: 500;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    box-shadow: 0 4px 12px rgba(120, 121, 255, 0.4);
    backdrop-filter: blur(10px);
  }

  .keyboard-hints kbd {
    display: inline-block;
    padding: 0.2rem 0.4rem;
    background: rgba(255, 255, 255, 0.25);
    border: 1px solid rgba(255, 255, 255, 0.3);
    border-radius: 4px;
    font-family: 'SF Mono', monospace;
    font-size: 0.75rem;
    color: white;
    box-shadow: 0 1px 0 rgba(0, 0, 0, 0.1);
  }

  /* Hide navigation elements on mobile */
  @media (max-width: 768px) {
    .nav-arrow {
      width: 40px;
      height: 40px;
    }

    .nav-arrow.prev {
      left: 5px;
    }

    .nav-arrow.next {
      right: 5px;
    }

    .keyboard-hints {
      display: none; /* Hide keyboard hints on mobile */
    }
  }
</style>
