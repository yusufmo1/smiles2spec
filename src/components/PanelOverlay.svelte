<script>
    import { focusedPanel, carouselIndex, isCarouselMode } from '../stores.js';
    import { fade, scale, fly } from 'svelte/transition';
    import { cubicOut, quartOut } from 'svelte/easing';
    import { onMount, onDestroy, tick } from 'svelte';
    import Plotly from 'plotly.js-dist-min';
    
    // Import Swiper and required modules
    import { register } from 'swiper/element/bundle';
    import 'swiper/css';
    import 'swiper/css/effect-coverflow';
    import 'swiper/css/navigation';
  
    // Register Swiper custom elements
    register();
  
    let swiperEl;
    let slideEl;  // Reference to the slide container
    let initialized = false;
    let panels = []; // Array to store panel elements
    
    // Helpers
    const close = () => {
      focusedPanel.set(null);
      isCarouselMode.set(false);
    };
    const next = () => {
      if (swiperEl && swiperEl.swiper) {
        swiperEl.swiper.slideNext();
        carouselIndex.set(swiperEl.swiper.activeIndex);
      }
    };
    const prev = () => {
      if (swiperEl && swiperEl.swiper) {
        swiperEl.swiper.slidePrev();
        carouselIndex.set(swiperEl.swiper.activeIndex);
      }
    };
  
    /* keyboard shortcuts */
    const handleKey = (e) => {
      if (!$focusedPanel) return;
      if (e.key === 'Escape') close();
      if (e.key === 'ArrowRight') next();
      if (e.key === 'ArrowLeft') prev();
    };
  
    // Global resize handler for Plotly plots
    const resizePlotlyCharts = () => {
      const activePlot = document.querySelector('.swiper-slide-active .js-plotly-plot');
      if (activePlot) {
        Plotly.Plots.resize(activePlot);
      }
    };
  
    // Handle resize to ensure proper centering
    const handleResize = () => {
      if (swiperEl && swiperEl.swiper) {
        swiperEl.swiper.update();
      }
      
      // Resize Plotly plots using the global handler
      resizePlotlyCharts();
    };
  
    // Initialize swiper once the component is mounted
    onMount(async () => {
      if (swiperEl) {
        // Configure Swiper using a single object
        const opts = {
          slidesPerView: 'auto',
          centeredSlides: true,
          effect: 'coverflow',
          grabCursor: true,
          initialSlide: $carouselIndex,
          watchSlidesProgress: true,
          loop: false,
          coverflowEffect: {
            rotate: 0,
            stretch: 120,
            depth: 250,
            modifier: 1,
            slideShadows: false
          }
        };
        
        // Set all options at once
        Object.entries(opts).forEach(([k, v]) => 
          swiperEl.setAttribute(k, typeof v === 'object' ? JSON.stringify(v) : v.toString())
        );
        
        // Initialize with delay to ensure DOM is fully ready
        setTimeout(async () => {
          if (!swiperEl.initialized) {
            await swiperEl.initialize();
          }
          
          if (swiperEl.swiper) {
            // Set responsive breakpoints
            swiperEl.swiper.params.breakpoints = {
              0: { coverflowEffect: {stretch: 60, depth: 120}},
              768: { coverflowEffect: {stretch: 120, depth: 250}},
              1280: { coverflowEffect: {stretch: 180, depth: 300}}
            };
            
            // Force update and go to the active slide
            swiperEl.swiper.update();
            swiperEl.swiper.slideTo($carouselIndex, 0);
            initialized = true;
          }
        }, 150);
        
        // Event listeners for swiper
        swiperEl.addEventListener('slidechange', () => {
          if (swiperEl.swiper) {
            // Update the active index
            carouselIndex.set(swiperEl.swiper.activeIndex);
            // Resize the plot in the active slide
            resizePlotlyCharts();
          }
        });
        
        // Add global resize listener for Plotly charts
        window.addEventListener('resize', resizePlotlyCharts);
        window.addEventListener('resize', handleResize);
      }
    });
    
    // Clean up
    onDestroy(() => {
      window.removeEventListener('resize', handleResize);
      window.removeEventListener('resize', resizePlotlyCharts);
      swiperEl?.swiper?.destroy(true, true);
    });
    
    // Handle transition end
    const handleTransitionEnd = () => {
      setTimeout(() => {
        if (swiperEl && swiperEl.swiper) {
          swiperEl.swiper.update();
          
          // Resize any plots
          resizePlotlyCharts();
        }
      }, 300);
    };

    // Function to handle the opening of the overlay
    const handleOverlayOpen = async () => {
      if ($focusedPanel) {
        // Collect all panels for the carousel
        const allPanels = Array.from(document.querySelectorAll('.panel'));
        panels = [];
        
        // Create slots for each panel
        for (let i = 0; i < allPanels.length; i++) {
          const slide = document.createElement('swiper-slide');
          const content = document.createElement('div');
          content.className = 'swiper-slide-content glass-card';
          slide.appendChild(content);
          swiperEl.appendChild(slide);
          
          // Clone the panel and its children deeply
          const panelClone = allPanels[i].cloneNode(true);
          
          // Add to our collection
          panels.push({ slide, content, panelClone });
        }
        
        // After DOM is updated
        await tick();
        
        // Move all panel clones to their slides
        for (let i = 0; i < panels.length; i++) {
          panels[i].content.appendChild(panels[i].panelClone);
        }
        
        // For each panel with a Plotly plot, we need to re-create the plot
        for (let i = 0; i < panels.length; i++) {
          const plotDiv = panels[i].panelClone.querySelector('.js-plotly-plot');
          if (plotDiv) {
            // Find the original plot div
            const originalPlot = allPanels[i].querySelector('.js-plotly-plot');
            if (originalPlot && originalPlot.data) {
              // Re-create the plot using the same data and layout
              Plotly.newPlot(plotDiv, originalPlot.data, 
                { ...originalPlot.layout, autosize: true });
            }
          }
        }
        
        // Initialize with the correct panel
        if (swiperEl && swiperEl.swiper) {
          swiperEl.swiper.slideTo($carouselIndex, 0, false);
        }
        
        // After all plots are created, resize the active one
        resizePlotlyCharts();
      }
    };
    
    // Function to handle the closing of the overlay
    const handleOverlayClose = () => {
      // Clear the panel clones
      panels = [];
    };
  </script>
  
  <svelte:window on:keydown={handleKey} />
  
  {#if $focusedPanel}
    <!-- Full-screen overlay backdrop -->
    <div 
      class="fullscreen-overlay" 
      role="dialog" 
      on:click|self={close} 
      on:keydown={(e) => e.key === 'Escape' && close()}
      aria-modal="true"
      transition:fade={{ duration: 400 }}
      on:introstart={handleOverlayOpen}
      on:outroend={handleOverlayClose}
    >
      <!-- App title remains visible -->
      <div class="overlay-title" transition:fly={{ y: -50, duration: 400, easing: quartOut }}>
        <h1 class="display text-gradient prosto-one-regular">SPECTRAL SIMULATION</h1>
      </div>
      
      <!-- Content container -->
      <div 
        class="overlay-content" 
        transition:scale={{ duration: 400, easing: cubicOut, start: 0.85 }}
        on:introend={handleTransitionEnd}
      >
        <!-- CAROUSEL VIEW -->
        <div class="swiper-container-wrapper">
          <swiper-container bind:this={swiperEl} class="swiper-container">
            <!-- Slides will be added dynamically in handleOverlayOpen -->
          </swiper-container>
        </div>
        
        <!-- Navigation arrows -->
        <button class="arrow left" on:click={prev} aria-label="Previous panel">◀</button>
        <button class="arrow right" on:click={next} aria-label="Next panel">▶</button>
        
        <!-- Close button -->
        <button class="control-button close" on:click={close} aria-label="Close">✕</button>
        
        <!-- Pagination indicators -->
        <div class="pagination-indicator">
          {#each panels as _, i}
            <div class="indicator-dot" class:active={i === $carouselIndex}></div>
          {/each}
        </div>
        
        <!-- Keyboard shortcuts indicator -->
        <div class="keyboard-shortcuts">
          <div class="shortcut"><kbd>Esc</kbd> Close</div>
          <div class="shortcut"><kbd>←</kbd><kbd>→</kbd> Navigate</div>
        </div>
      </div>
      
      <!-- Subtle branding watermark -->
      <div class="overlay-watermark">
        Spectral Simulation Tool
      </div>
    </div>
  {/if}
  
  <style>
    /* Ensure consistent box-sizing for all elements */
    :global(.fullscreen-overlay *) {
      box-sizing: border-box;
    }
  
    .fullscreen-overlay {
      position: fixed;
      inset: 0;
      z-index: 10000;
      background: rgba(247, 247, 248, 0.95);
      backdrop-filter: blur(25px) saturate(180%);
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      padding: 0;
      width: 100vw;
      height: 100vh;
      animation: reveal 0.5s cubic-bezier(0.16, 1, 0.3, 1) forwards;
      box-sizing: border-box;
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
        transform: translateY(40px) scale(0.95);
      }
      to { 
        opacity: 1;
        transform: translateY(0) scale(1);
      }
    }
    
    .overlay-title {
      position: absolute;
      top: 2rem;
      text-align: center;
      z-index: 10001;
    }
    
    .overlay-title h1 {
      font: 600 clamp(2.5rem, 4vw, 3.75rem)/1 'Prosto One', sans-serif;
      margin: 0;
      letter-spacing: 0.05em;
    }
    
    .overlay-content {
      position: relative;
      width: 100vw;
      height: 80vh;
      margin-top: 5rem;
      display: flex;
      align-items: center;
      justify-content: center;
      box-sizing: border-box;
    }
    
    /* Swiper container wrapper */
    .swiper-container-wrapper {
      width: 100%;
      height: 100%;
      position: relative;
      overflow: visible;
      display: flex;
      justify-content: center;
      align-items: center;
    }
    
    /* Swiper container */
    :global(.swiper-container) {
      width: 90% !important;
      height: 100% !important;
      position: relative;
      overflow: visible !important;
    }
    
    /* Override Swiper's internal styles for the coverflow effect */
    :global(.swiper-wrapper) {
      align-items: center;
      display: flex;
      height: 100%;
    }
    
    /* Individual slides */
    :global(.swiper-slide) {
      height: 100% !important;
      max-height: 90%;
      display: flex;
      align-items: center;
      justify-content: center;
    }
    
    /* Centre card pops, side cards fade */
    :global(.swiper-slide:not(.swiper-slide-active)) {
      filter: brightness(.8) blur(1px) grayscale(.15);
      transition: filter .25s var(--transition-smooth);
    }
    
    /* Individual slide content styling */
    :global(.swiper-slide-content) {
      display: flex;
      flex-direction: column;
      height: 100%;
      width: 100%;
      overflow: hidden !important;
      border-radius: var(--enforce-pill);
      box-shadow: 
        0 25px 60px rgba(0, 0, 0, 0.15),
        0 0 0 1px rgba(255, 255, 255, 0.5);
      transition: all var(--transition-smooth);
    }
    
    :global(.swiper-slide-content .panel) {
      display: flex;
      flex-direction: column;
      flex: 1;
    }
    
    :global(.swiper-slide-content .panel-content) {
      flex: 1;
      display: flex;
      flex-direction: column;
    }
    
    :global(.swiper-slide-content .panel-content > *) {
      flex: 1 1 auto;
    }
    
    /* Ensure plots and charts fill the space */
    :global(.swiper-slide-content .panel-content svg),
    :global(.swiper-slide-content .panel-content canvas),
    :global(.swiper-slide-content .panel-content img) {
      flex: 1;
      max-width: 100%;
      object-fit: contain;
    }
    
    /* Pagination indicator styles */
    .pagination-indicator {
      position: absolute;
      bottom: -3rem;
      left: 50%;
      transform: translateX(-50%);
      display: flex;
      gap: 8px;
      z-index: 10;
    }
    
    .indicator-dot {
      width: 8px;
      height: 8px;
      border-radius: 50%;
      background: var(--accent-soft);
      transition: all 0.3s ease;
    }
    
    .indicator-dot.active {
      background: var(--accent);
      transform: scale(1.25);
    }
    
    /* Navigation arrows */
    .arrow {
      position: absolute;
      top: 50%;
      transform: translateY(-50%);
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
    }
    
    .arrow:hover {
      transform: translateY(-50%) scale(1.1);
      background: var(--accent);
      color: white;
      box-shadow: 0 6px 20px rgba(120, 121, 255, 0.4);
    }
    
    /* Position arrows at fixed positions relative to the center */
    .arrow.left { left: 5%; }
    .arrow.right { right: 5%; }
  
    /* Control button */
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
    
    /* Keyboard shortcuts */
    .keyboard-shortcuts {
      position: absolute;
      bottom: -5rem;
      left: 50%;
      transform: translateX(-50%);
      display: flex;
      gap: 1.5rem;
      color: var(--text-secondary);
      font-size: 0.85rem;
    }
    
    .shortcut {
      display: flex;
      align-items: center;
      gap: 0.5rem;
    }
    
    kbd {
      display: inline-block;
      padding: 0.2rem 0.5rem;
      background: rgba(255, 255, 255, 0.8);
      border: 1px solid rgba(0, 0, 0, 0.1);
      border-radius: 4px;
      box-shadow: 0 2px 0 rgba(0, 0, 0, 0.1);
      font-family: 'SF Mono', monospace;
      font-size: 0.8rem;
    }
    
    /* Watermark */
    .overlay-watermark {
      position: absolute;
      bottom: 1rem;
      left: 1rem;
      font-size: 0.8rem;
      opacity: 0.4;
      color: var(--text-tertiary);
    }
    
    /* Responsive adjustments */
    @media (max-width: 1280px) {
      .arrow.left { left: 5%; }
      .arrow.right { right: 5%; }
    }
    
    @media (max-width: 768px) {
      .overlay-content {
        height: 75vh;
      }
      
      .arrow.left { left: 5%; }
      .arrow.right { right: 5%; }
      
      .keyboard-shortcuts {
        flex-direction: column;
        align-items: center;
        gap: 0.75rem;
      }
      
      :global(.swiper-slide) {
        height: 80% !important;
      }
    }
    
    /* Accessibility - reduced motion */
    @media (prefers-reduced-motion: reduce) {
      :global(swiper-container) {
        --swiper-transition-duration: 0ms !important;
      }
      
      .fullscreen-overlay {
        animation: none;
      }
    }
  </style>