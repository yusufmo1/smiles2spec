<script lang="ts">
  export let title: string = '';
  export let description: string = '';
  export let layout: 'single' | 'grid-2x2' | 'grid-1x4' | 'flex' = 'single';
  export let compact: boolean = false;
</script>

<div class="inner-container" class:compact>
  {#if title || description}
    <div class="inner-header">
      {#if title}
        <h3 class="inner-title">{title}</h3>
      {/if}
      {#if description}
        <p class="inner-description" class:compact>{description}</p>
      {/if}
    </div>
  {/if}

  <div
    class="inner-content"
    class:grid-2x2={layout === 'grid-2x2'}
    class:grid-1x4={layout === 'grid-1x4'}
    class:flex={layout === 'flex'}
    class:compact
  >
    <slot />
  </div>
</div>

<style>
  .inner-container {
    height: 100%;
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
  }

  /* Compact mode adjustments */
  .inner-container.compact {
    gap: 1rem; /* Reduced from 1.5rem */
  }

  .inner-header {
    text-align: center;
  }

  .inner-title {
    font-size: 1.25rem;
    font-weight: 600;
    color: var(--accent);
    margin: 0 0 0.5rem 0;
  }

  .inner-description {
    font-size: 1rem;
    color: var(--text-secondary);
    line-height: 1.6;
    margin: 0;
  }

  .inner-description.compact {
    margin: 0;
    font-size: 0.95rem; /* Slightly smaller */
    line-height: 1.5; /* Tighter */
  }

  .inner-content {
    flex: 1;
    display: flex;
    flex-direction: column;
  }

  .inner-content.grid-2x2 {
    display: grid;
    grid-template-columns: 1fr 1fr;
    grid-template-rows: 1fr 1fr;
    gap: 1.5rem;
  }

  .inner-content.compact.grid-2x2 {
    gap: 1rem; /* Reduced from 1.5rem */
  }

  .inner-content.grid-1x4 {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1rem;
  }

  .inner-content.compact.grid-1x4 {
    gap: 0.75rem; /* Reduced from 1rem */
  }

  .inner-content.flex {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
  }

  @media (max-width: 1024px) {
    .inner-content.grid-2x2 {
      grid-template-columns: 1fr;
      grid-template-rows: repeat(4, 1fr);
    }

    .inner-content.grid-1x4 {
      grid-template-columns: 1fr 1fr;
      grid-template-rows: 1fr 1fr;
    }
  }

  @media (max-width: 640px) {
    .inner-content.grid-1x4 {
      grid-template-columns: 1fr;
      grid-template-rows: repeat(4, 1fr);
    }
  }
</style>
