<script lang="ts">
  import { focusedPanel } from '$lib/stores/index';
  import PanelRenderer from './PanelRenderer.svelte';
  import { isInteractiveElement } from './utils/navigation';
  import { devLog } from '$lib/utils/performance';

  export let id: string;
  export let title = '';
  export let component: any;
  export let props = {};
  export let clickable = true;
  export let isCarousel = false;
  export let compact = false;

  // Log when props change
  $: {
    devLog.log(`Panel ${id} props updated:`, props);
  }

  let panelEl: HTMLDivElement;

  function handleClick(event: MouseEvent) {
    if (!clickable) return;

    // In carousel mode, allow interactive elements to work but prevent panel expansion
    if (isCarousel) {
      // Check if click is on an interactive element
      let target = event.target as HTMLElement;
      while (target && target !== panelEl) {
        if (isInteractiveElement(target)) {
          // Let the interactive element handle the click normally
          return;
        }
        target = target.parentElement as HTMLElement;
      }
      // If clicked on non-interactive area in carousel, do nothing
      return;
    }

    // Enhanced interactive element detection for grid mode
    let target = event.target as HTMLElement;
    while (target && target !== panelEl) {
      if (isInteractiveElement(target)) {
        event.stopPropagation();
        return;
      }
      target = target.parentElement as HTMLElement;
    }

    // ✅ SIMPLIFIED: Only set focusedPanel - everything else derives from this
    focusedPanel.set(id);
  }

  function handleKeyDown(event: KeyboardEvent) {
    if (!clickable || isCarousel) return; // Don't handle keys in carousel mode

    // Don't capture space or enter if the event target is a form control
    const target = event.target as HTMLElement;
    if (
      target.tagName === 'TEXTAREA' ||
      target.tagName === 'INPUT' ||
      target.getAttribute('contenteditable') === 'true'
    ) {
      return;
    }

    if (event.key === 'Enter' || event.key === ' ') {
      handleClick(event as any);
      event.preventDefault();
    }
  }
</script>

<div
  bind:this={panelEl}
  class="glass-card panel hover-lift"
  class:clickable
  class:carousel={isCarousel}
  class:compact
  data-panel-id={id}
  on:click={handleClick}
  on:keydown={handleKeyDown}
  role={clickable && !isCarousel ? 'button' : 'region'}
  aria-label={title ? `${isCarousel ? 'View' : 'Expand'} ${title} panel` : 'Panel'}
>
  {#if title}
    <div class="title-pill">
      <span class="title-dot"></span>
      <h2 class="prosto-one-regular">{title}</h2>
    </div>
  {/if}

  <div class="panel-content" class:compact>
    <PanelRenderer {component} {props} {isCarousel} />
  </div>
</div>

<style>
  /* Import only needed animations */
  .hover-lift {
    transition:
      transform var(--transition-smooth),
      box-shadow var(--transition-smooth);
  }

  .hover-lift:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 25px rgba(0, 0, 0, 0.08);
  }
  .panel {
    display: flex;
    flex-direction: column;
    height: 440px;
    width: 100%;
    overflow: hidden !important;
    border-radius: var(--enforce-pill);
    position: relative;
  }

  /* ✅ FIXED: Simplified carousel mode */
  .panel.carousel {
    height: 100%;
    max-height: 100%;
    cursor: default;
  }

  .panel.clickable:not(.carousel) {
    cursor: pointer;
  }

  /* Disable hover lift in carousel mode */
  .panel.carousel {
    transition: none !important;
    transform: none !important;
  }

  .panel:focus {
    outline: 2px solid var(--accent);
    outline-offset: 2px;
  }

  .panel:focus:not(:focus-visible) {
    outline: none;
  }

  /* New title pill styling */
  .title-pill {
    width: calc(100% - 2rem);
    margin: 1rem auto 1.5rem;
    background: var(--surface-glass);
    border: 1px solid rgba(255, 255, 255, 0.9);
    border-radius: var(--enforce-pill);
    padding: 0.75rem 1.5rem;
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 1rem;
    backdrop-filter: blur(15px);
    box-shadow:
      0 4px 12px rgba(0, 0, 0, 0.08),
      0 1px 3px rgba(255, 255, 255, 0.6) inset;
    flex-shrink: 0;
    position: relative;
    z-index: 1;
  }

  .title-dot {
    width: 10px;
    height: 10px;
    border-radius: 50%;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    box-shadow: 0 0 12px var(--accent-soft);
    flex-shrink: 0;
  }

  h2 {
    margin: 0;
    font-size: 14px;
    font-weight: 600;
    background: linear-gradient(135deg, var(--text-primary), var(--text-secondary));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    letter-spacing: 0.8px;
    text-align: center;
  }

  .panel-content {
    flex: 1;
    padding: 0 1.75rem 1.75rem;
    overflow: auto;
    position: relative;
    z-index: 1;
    box-sizing: border-box;
    min-height: 0; /* Allow flex child to shrink */
  }

  .panel-content.compact {
    padding: 0 1.25rem 1.25rem;
  }

  @media (max-width: 640px) {
    .panel-content.compact {
      padding: 0 1rem 1rem;
    }

    .title-pill {
      width: calc(100% - 1.5rem);
      margin: 0.75rem auto 1rem;
      padding: 0.65rem 1.25rem;
    }

    h2 {
      font-size: 13px;
    }
  }

  /* Removing the inner glow effect to achieve uniform transparency */
</style>
