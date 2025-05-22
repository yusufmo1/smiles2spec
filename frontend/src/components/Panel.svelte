<script>
  import { focusedPanel, carouselIndex } from '../stores.js';
  import PanelRenderer from './PanelRenderer.svelte';

  export let id;
  export let title = '';
  export let component;
  export let props = {};
  export let clickable = true;
  export let isCarousel = false;

  let panelEl;

  function handleClick(event) {
    if (!clickable || $focusedPanel) return;

    let target = event.target;
    while (target && target !== panelEl) {
      if (isInteractiveElement(target)) {
        event.stopPropagation();
        return;
      }
      target = target.parentElement;
    }

    focusedPanel.set(id);

    const panelOrder = ['spectrum', 'fragments', 'structure', 'peaks', 'console', 'chat'];
    const index = panelOrder.indexOf(id);
    if (index !== -1) {
      carouselIndex.set(index);
    }
  }

  function isInteractiveElement(element) {
    return element.tagName === 'BUTTON' ||
           element.tagName === 'TEXTAREA' ||
           element.tagName === 'INPUT' ||
           element.tagName === 'SELECT' ||
           element.getAttribute('contenteditable') === 'true' ||
           element.getAttribute('role') === 'button' ||
           element.classList.contains('clickable') ||
           element.closest('button, textarea, input, select, [contenteditable="true"], [role="button"]');
  }
</script>

<div
  bind:this={panelEl}
  class="glass-card panel"
  class:clickable
  class:carousel={isCarousel}
  data-panel-id={id}
  on:click={handleClick}
  tabindex={clickable ? "0" : "-1"}
  role={clickable ? "button" : "region"}
  aria-label={title ? `${isCarousel ? 'View' : 'Expand'} ${title} panel` : 'Panel'}
>
  {#if title}
    <header>
      <div class="title-wrapper">
        <span class="title-dot"></span>
        <h2 class="prosto-one-regular">{title}</h2>
      </div>
    </header>
  {/if}

  <div class="panel-content">
    <PanelRenderer {component} {props} {isCarousel} />
  </div>
</div>

<style>
  .panel {
    display: flex;
    flex-direction: column;
    height: 440px;
    width: 100%;
    overflow: hidden !important;
    border-radius: var(--enforce-pill);
    position: relative;
    transition: transform 0.2s, box-shadow 0.2s;
  }

  .panel.carousel {
    height: 100%;
    max-height: 90vh;
  }

  .panel.clickable {
    cursor: pointer;
  }

  .panel.clickable:hover {
    transform: translateY(-2px);
    box-shadow: 0 12px 30px rgba(0, 0, 0, 0.1);
  }

  header {
    font-size: 13px;
    font-weight: 600;
    color: var(--text-secondary);
    padding: 1.5rem 1.75rem;
    border-bottom: 1px solid var(--surface-stroke);
    letter-spacing: 0.8px;
    display: flex;
    justify-content: center;
    align-items: center;
    backdrop-filter: blur(10px);
    position: relative;
    z-index: 1;
    height: 65px;
    box-sizing: border-box;
  }

  .title-wrapper {
    display: flex;
    align-items: center;
  }

  .title-dot {
    width: 10px;
    height: 10px;
    border-radius: 50%;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    margin-right: 1rem;
    box-shadow: 0 0 12px var(--accent-soft);
  }

  h2 {
    margin: 0;
    font-size: 14px;
    font-weight: 600;
    background: linear-gradient(135deg, var(--text-primary), var(--text-secondary));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
  }

  .panel-content {
    flex: 1;
    padding: 1.75rem;
    overflow: auto;
    height: calc(340px - 65px);
    position: relative;
    z-index: 1;
    box-sizing: border-box;
  }

  .panel-content::after {
    content: '';
    position: absolute;
    inset: 0;
    box-shadow: 0 0 25px rgba(255, 255, 255, 0.6) inset;
    pointer-events: none;
    z-index: -1;
  }
</style>
