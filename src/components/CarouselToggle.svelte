<script>
  import { createEventDispatcher } from 'svelte';
  export let value = false;
  const dispatch = createEventDispatcher();
  function toggle() { 
    value = !value; 
    dispatch('change', value); 
  }
</script>

<button class="toggle" on:click={toggle} aria-label="Toggle carousel view">
  <span class:active={!value}>Grid</span>
  <span class:active={value}>Carousel</span>
  <div class="thumb" class:right={value}></div>
</button>

<style>
  .toggle {
    position: relative;
    display: flex;
    background: var(--toggle-bg, rgba(255, 255, 255, 0.6));
    border: 1px solid var(--surface-stroke);
    border-radius: var(--radius-pill);
    padding: 0.25rem;
    width: 160px;
    height: 36px;
    cursor: pointer;
    box-shadow: var(--shadow-sm);
    transition: var(--transition-smooth);
  }

  .toggle:hover {
    background: var(--surface-glass-hover);
  }

  span {
    flex: 1;
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 2;
    font-size: 0.875rem;
    font-weight: 500;
    color: var(--text-secondary);
    transition: var(--transition-smooth);
    padding: 0 0.5rem;
  }

  span.active {
    color: var(--text-primary);
    font-weight: 600;
  }

  .thumb {
    position: absolute;
    top: 3px;
    left: 3px;
    width: calc(50% - 6px);
    height: calc(100% - 6px);
    background: var(--toggle-thumb, var(--accent));
    opacity: 0.15;
    border-radius: calc(var(--radius-pill) - 4px);
    transition: transform var(--transition-elastic);
    transform: translateX(0);
  }

  .thumb.right {
    transform: translateX(100%);
  }
</style> 