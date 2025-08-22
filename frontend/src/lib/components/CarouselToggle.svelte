<script lang="ts">
  import { createEventDispatcher } from 'svelte';
  import { DynamicIcon } from '$lib/components/icons';

  export let value = false;

  const dispatch = createEventDispatcher<{
    change: boolean;
  }>();

  function toggle() {
    value = !value;
    dispatch('change', value);
  }
</script>

<div class="toggle-wrapper">
  <button
    class="dual-toggle"
    on:click={toggle}
    title={value ? 'Switch to grid view' : 'Switch to carousel view'}
  >
    <div class="option grid" class:active={!value} aria-label="Grid view">
      <DynamicIcon name="GridIcon" size={18} />
    </div>
    <div class="option carousel" class:active={value} aria-label="Carousel view">
      <DynamicIcon name="CarouselIcon" size={18} />
    </div>
    <div class="indicator" class:right={value}></div>
  </button>
</div>

<style>
  .toggle-wrapper {
    display: flex;
    align-items: center;
  }

  .dual-toggle {
    position: relative;
    display: flex;
    align-items: center;
    width: 84px;
    height: 38px;
    border-radius: 19px;
    background: rgba(255, 255, 255, 0.6);
    border: 1px solid rgba(0, 0, 0, 0.1);
    padding: 2px;
    cursor: pointer;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
    transition: all 0.2s ease;
  }

  .dual-toggle:hover {
    background: rgba(255, 255, 255, 0.8);
  }

  .dual-toggle:focus {
    outline: 2px solid #7879ff;
    outline-offset: 2px;
  }

  .option {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 40px;
    height: 32px;
    border-radius: 16px;
    position: relative;
    z-index: 2;
    transition: all 0.2s ease;
    color: #666;
  }

  .option.active {
    color: #333;
  }

  .indicator {
    position: absolute;
    top: 3px;
    left: 3px;
    width: 40px;
    height: 32px;
    background: #7879ff;
    opacity: 0.15;
    border-radius: 16px;
    transition: transform 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    transform: translateX(0);
    z-index: 1;
  }

  .indicator.right {
    transform: translateX(38px);
  }

  /* Responsive adjustments */
  @media (max-width: 640px) {
    .dual-toggle {
      width: 76px;
      height: 34px;
    }

    .option {
      width: 36px;
      height: 28px;
    }

    .indicator {
      width: 36px;
      height: 28px;
    }

    .indicator.right {
      transform: translateX(34px);
    }
  }
</style>
