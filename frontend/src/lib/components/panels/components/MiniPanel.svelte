<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let icon: string = '';
  export let iconComponent: any = null;
  export let title: string = '';
  export let description: string = '';
  export let variant: 'default' | 'highlight' | 'accent' = 'default';
  export let compact: boolean = false;
</script>

<div
  class="mini-panel glass-card hover-lift"
  class:highlight={variant === 'highlight'}
  class:accent={variant === 'accent'}
  class:compact
>
  {#if iconComponent}
    <div class="mini-icon" class:compact>
      <DynamicIcon name={iconComponent as any} size={compact ? 32 : 40} color="var(--accent)" />
    </div>
  {:else if icon}
    <div class="mini-icon" class:compact>{icon}</div>
  {/if}
  <div class="mini-content">
    <h4 class="mini-title" class:compact>{title}</h4>
    <p class="mini-description" class:compact>{description}</p>
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
  .mini-panel {
    padding: 1.5rem;
    border-radius: var(--radius-lg);
    background: rgba(255, 255, 255, 0.6);
    border: 1px solid rgba(255, 255, 255, 0.8);
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    height: 100%;
  }

  /* Compact mode adjustments */
  .mini-panel.compact {
    padding: 1.25rem; /* Reduced from 1.5rem */
  }

  .mini-panel:hover {
    background: rgba(255, 255, 255, 0.8);
  }

  .mini-panel.highlight {
    background: var(--accent-soft);
    border-color: rgba(120, 121, 255, 0.3);
  }

  .mini-panel.accent {
    background: linear-gradient(135deg, var(--accent-soft), rgba(120, 121, 255, 0.1));
    border-color: var(--accent);
  }

  .mini-icon {
    font-size: 2.5rem;
    margin-bottom: 1rem;
    line-height: 1;
  }

  .mini-icon.compact {
    font-size: 2rem; /* Reduced from 2.5rem */
    margin-bottom: 0.75rem; /* Reduced from 1rem */
  }

  .mini-title {
    font-size: 1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.5rem 0;
  }

  .mini-title.compact {
    font-size: 0.95rem; /* Reduced from 1rem */
    margin: 0 0 0.4rem 0; /* Reduced bottom margin */
  }

  .mini-description {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
    margin: 0;
  }

  .mini-description.compact {
    font-size: 0.85rem; /* Reduced from 0.9rem */
    line-height: 1.4; /* Tighter from 1.5 */
  }

  @media (max-width: 768px) {
    .mini-panel {
      padding: 1.25rem;
    }

    .mini-panel.compact {
      padding: 1rem; /* Even more compact on mobile */
    }

    .mini-icon {
      font-size: 2rem;
      margin-bottom: 0.75rem;
    }

    .mini-icon.compact {
      font-size: 1.75rem;
      margin-bottom: 0.5rem;
    }

    .mini-title {
      font-size: 0.9rem;
    }

    .mini-description {
      font-size: 0.85rem;
    }
  }
</style>
