<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let items: Array<{
    id: string;
    title: string;
    value: string | number;
    icon?: any;
    color?: string;
    size?: 'small' | 'medium' | 'large';
  }> = [];

  let hoveredItem: string | null = null;
</script>

<div class="bento-grid">
  {#each items as item}
    <div
      class="bento-item {item.size || 'medium'}"
      class:hovered={hoveredItem === item.id}
      on:mouseenter={() => (hoveredItem = item.id)}
      on:mouseleave={() => (hoveredItem = null)}
      style="--item-color: {item.color || 'var(--accent)'}"
      role="button"
      tabindex="0"
      on:keydown={(e) =>
        e.key === 'Enter' && (hoveredItem = hoveredItem === item.id ? null : item.id)}
    >
      <div class="bento-content">
        {#if item.icon}
          <div class="bento-icon">
            <DynamicIcon name={item.icon as any} size={24} />
          </div>
        {/if}
        <div class="bento-value">{item.value}</div>
        <div class="bento-title">{item.title}</div>
      </div>
      <div class="bento-glow"></div>
    </div>
  {/each}
</div>

<style>
  .bento-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
    gap: 1rem;
    height: 100%;
    padding: 1rem;
  }

  .bento-item {
    position: relative;
    background: rgba(255, 255, 255, 0.6);
    border-radius: 24px;
    padding: 1.5rem;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    transition: all 0.3s ease;
    overflow: hidden;
    min-height: 120px;
  }

  .bento-item.small {
    grid-column: span 1;
  }

  .bento-item.medium {
    grid-column: span 2;
  }

  .bento-item.large {
    grid-column: span 3;
    grid-row: span 2;
  }

  .bento-item:hover {
    transform: translateY(-4px);
    background: rgba(255, 255, 255, 0.8);
    box-shadow: 0 12px 32px rgba(0, 0, 0, 0.1);
  }

  .bento-content {
    position: relative;
    z-index: 2;
    text-align: center;
  }

  .bento-icon {
    color: var(--item-color);
    margin-bottom: 0.75rem;
    opacity: 0.8;
  }

  .bento-value {
    font-size: 2.5rem;
    font-weight: 700;
    color: var(--item-color);
    line-height: 1;
    margin-bottom: 0.25rem;
  }

  .bento-title {
    font-size: 0.85rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.5px;
  }

  .bento-glow {
    position: absolute;
    inset: 0;
    background: radial-gradient(circle at center, var(--item-color), transparent 70%);
    opacity: 0;
    transition: opacity 0.3s ease;
    z-index: 1;
  }

  .bento-item:hover .bento-glow {
    opacity: 0.1;
  }
</style>
