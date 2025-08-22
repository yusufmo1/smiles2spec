<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import type { IconName } from '$lib/components/icons';

  export let icon: IconName = 'LoadingIcon';
  export let message: string = 'Awaiting data';
  export let hint: string = '';
  export let variant: 'default' | 'error' | 'empty' = 'default';
</script>

<div class="placeholder-state" class:error={variant === 'error'} class:empty={variant === 'empty'}>
  <div class="placeholder-icon">
    <DynamicIcon
      name={icon}
      size={48}
      color={variant === 'error' ? 'var(--error)' : 'var(--text-secondary)'}
    />
  </div>
  <p class="placeholder-message">{message}</p>
  {#if hint}
    <p class="placeholder-hint">{hint}</p>
  {/if}
</div>

<style>
  .placeholder-state {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    height: 100%;
    width: 100%;
    padding: 2rem;
    text-align: center;
    min-height: 200px;
  }

  .placeholder-icon {
    margin-bottom: 1rem;
    opacity: 0.6;
    animation: pulse 2s ease-in-out infinite;
  }

  .placeholder-message {
    font-size: 1rem;
    color: var(--text-secondary);
    margin: 0 0 0.5rem 0;
    font-weight: 500;
  }

  .placeholder-hint {
    font-size: 0.875rem;
    color: var(--text-tertiary);
    margin: 0;
    max-width: 300px;
  }

  .placeholder-state.error .placeholder-message {
    color: var(--error);
  }

  .placeholder-state.empty .placeholder-icon {
    opacity: 0.4;
    animation: none;
  }

  @keyframes pulse {
    0%,
    100% {
      opacity: 0.6;
    }
    50% {
      opacity: 0.3;
    }
  }

  @media (max-width: 640px) {
    .placeholder-state {
      padding: 1.5rem;
    }

    .placeholder-message {
      font-size: 0.9rem;
    }

    .placeholder-hint {
      font-size: 0.8rem;
    }
  }
</style>
