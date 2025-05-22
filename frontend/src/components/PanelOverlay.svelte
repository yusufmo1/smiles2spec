<script>
  import { focusedPanel, carouselIndex, isCarouselMode } from '../stores.js';
  import { fade, scale } from 'svelte/transition';
  import PanelCarousel from './PanelCarousel.svelte';

  const close = () => {
    focusedPanel.set(null);
    isCarouselMode.set(false);
  };

  const next = () => {
    carouselIndex.update(i => Math.min(i + 1, 5));
  };

  const prev = () => {
    carouselIndex.update(i => Math.max(i - 1, 0));
  };

  function handleKey(e) {
    if (!$focusedPanel) return;
    if (e.key === 'Escape') close();
    if (e.key === 'ArrowRight') next();
    if (e.key === 'ArrowLeft') prev();
  }

  function handleTransitionEnd() {
    setTimeout(() => {
      const plots = document.querySelectorAll('.swiper-slide-active .js-plotly-plot');
      plots.forEach(plot => {
        if (window.Plotly) {
          window.Plotly.Plots.resize(plot);
        }
      });
    }, 100);
  }
</script>

<svelte:window on:keydown={handleKey} />

{#if $focusedPanel}
  <div class="fullscreen-overlay" on:click|self={close} transition:fade={{ duration: 400 }} on:introend={handleTransitionEnd}>
    <div class="overlay-title">
      <h1 class="display text-gradient prosto-one-regular">SPECTRAL SIMULATION</h1>
    </div>

    <div class="overlay-content" transition:scale={{ duration: 400, start: 0.85 }}>
      <PanelCarousel />

      <button class="arrow left" on:click={prev} disabled={$carouselIndex === 0}>◀</button>
      <button class="arrow right" on:click={next} disabled={$carouselIndex === 5}>▶</button>

      <button class="control-button close" on:click={close}>✕</button>

      <div class="pagination-indicator">
        {#each Array(6) as _, i}
          <div class="indicator-dot" class:active={i === $carouselIndex}></div>
        {/each}
      </div>
    </div>
  </div>
{/if}

<style>
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
  }

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
    z-index: 10;
    transition: all 0.2s;
  }

  .arrow:hover:not(:disabled) {
    background: var(--accent);
    color: white;
    transform: translateY(-50%) scale(1.1);
  }

  .arrow:disabled {
    opacity: 0.3;
    cursor: not-allowed;
  }

  .arrow.left { left: 5%; }
  .arrow.right { right: 5%; }

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
    transition: all 0.2s;
  }

  .control-button:hover {
    background: #ff453a;
    color: white;
    transform: scale(1.1);
  }

  .pagination-indicator {
    position: absolute;
    bottom: -3rem;
    left: 50%;
    transform: translateX(-50%);
    display: flex;
    gap: 8px;
  }

  .indicator-dot {
    width: 8px;
    height: 8px;
    border-radius: 50%;
    background: var(--accent-soft);
    transition: all 0.3s;
  }

  .indicator-dot.active {
    background: var(--accent);
    transform: scale(1.25);
  }
</style>
