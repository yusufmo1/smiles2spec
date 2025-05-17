<script>
  import { focusedPanel } from '../stores.js';
  import { fade, scale } from 'svelte/transition';
  import { spring } from 'svelte/motion';
  import { cubicOut } from 'svelte/easing';

  // spring keeps it buttery-smooth when window resizes during focus
  const xy = spring({ x: 0, y: 0 }, { stiffness: 0.15, damping: 0.4 });

  /** close overlay on backdrop click or Esc */
  function close() {
    focusedPanel.set(null);
  }

  // close on Esc
  const handleKey = (e) => $focusedPanel && e.key === 'Escape' && close();
</script>

<svelte:window on:keydown={handleKey} />

{#if $focusedPanel}
  <!-- backdrop -->
  <div
    class="backdrop"
    on:click|self={close}
    transition:fade={{ duration: 160 }}
  />

  <!-- the zoomed panel itself -->
  <div
    class="overlay"
    transition:scale={{ duration: 240, easing: cubicOut }}
    style="transform-origin: center;"
  >
    <!-- slot contains the actual Panel markup that we clone -->
    {@html $focusedPanel}
    <button class="close" aria-label="Close" on:click={close}>✕</button>
  </div>
{/if}

<style>
  .backdrop {
    position: fixed;
    inset: 0;
    backdrop-filter: blur(12px);
    background: rgba(247, 247, 248, 0.55);
    z-index: 1000;
  }
  .overlay {
    position: fixed;
    inset: 0;
    margin: auto;
    width: min(90vw, 960px);
    height: min(90vh, 680px);
    z-index: 1001;
    overflow: hidden;
    border-radius: var(--enforce-pill);
    box-shadow: 0 25px 55px rgba(0,0,0,0.12);
  }
  .close {
    position: absolute;
    top: 14px;
    right: 20px;
    background: var(--accent-soft);
    border: none;
    border-radius: var(--enforce-pill);
    width: 32px;
    height: 32px;
    font-size: 1rem;
    cursor: pointer;
  }
</style> 