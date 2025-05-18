<script>
  import { focusedPanel, carouselIndex } from '../stores.js';
  import { onMount } from 'svelte';

  export let title = '';
  let htmlString;          // will hold the rendered markup
  let panelEl;             // bind:this element

  // after mount grab our outerHTML so we can clone into overlay
  onMount(() => {
    htmlString = panelEl.outerHTML;
  });

  function handleClick(event) {
    // block if some other panel is already focused
    if ($focusedPanel) return;

    // Check if the click target is a button or a child of a button
    let target = event.target;
    while (target && target !== panelEl) {
      if (target.tagName === 'BUTTON' || 
          target.getAttribute('role') === 'button' ||
          target.classList.contains('button') ||
          target.closest('button')) {
        // If clicking on a button, do not enter carousel mode
        event.stopPropagation();
        return;
      }
      target = target.parentElement;
    }

    /* Send the actual element instead of outerHTML */
    focusedPanel.set(panelEl);
  }
  
  function handleKeyDown(event) {
    if (event.key === 'Enter' || event.key === ' ') {
      handleClick(event);
      event.preventDefault();
    }
  }
</script>

<div 
  bind:this={panelEl} 
  class="glass-card panel" 
  on:click={handleClick}
  on:keydown={handleKeyDown}
  tabindex="0"
  role="button"
  aria-label={title ? `Expand ${title} panel` : 'Expand panel'}
>
  {#if title}
    <header>
      <div class="title-wrapper">
        <span class="title-dot"></span>
        <h2>{title}</h2>
      </div>
    </header>
  {/if}
  <div class="panel-content">
    <slot />
  </div>
</div>

<style>
  .panel {
    display: flex;
    flex-direction: column;
    height: 100%;
    width: 100%;
    overflow: hidden !important;
    min-height: 340px;
    border-radius: var(--enforce-pill);
    position: relative;
    cursor: pointer;
    transition: transform 0.2s, box-shadow 0.2s;
  }
  
  .panel:hover {
    transform: translateY(-2px);
    box-shadow: 0 12px 30px rgba(0, 0, 0, 0.1);
  }
  
  .panel:focus {
    outline: 2px solid var(--accent);
    outline-offset: 2px;
  }
  
  .panel:focus:not(:focus-visible) {
    outline: none;
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
    height: 100%;
    position: relative;
    z-index: 1;
  }
  
  /* Subtle inner glow effect */
  .panel-content::after {
    content: '';
    position: absolute;
    inset: 0;
    box-shadow: 0 0 25px rgba(255, 255, 255, 0.6) inset;
    pointer-events: none;
    z-index: -1;
  }
</style> 