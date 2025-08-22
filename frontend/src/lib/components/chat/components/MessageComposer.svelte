<script lang="ts">
  import { createEventDispatcher } from 'svelte';

  const dispatch = createEventDispatcher();

  export let disabled = false;
  export let loading = false;

  let message = '';

  function handleSubmit() {
    if (!message.trim() || disabled || loading) return;

    dispatch('send', { content: message.trim() });
    message = '';
  }

  function handleKeyDown(event: KeyboardEvent) {
    if (event.key === 'Enter' && !event.shiftKey) {
      event.preventDefault();
      handleSubmit();
    }
  }
</script>

<div class="message-composer">
  <div class="input-container">
    <textarea
      bind:value={message}
      on:keydown={handleKeyDown}
      placeholder="Ask about this spectrum..."
      {disabled}
      rows="1"
    ></textarea>
    <button
      type="button"
      on:click={handleSubmit}
      disabled={disabled || loading || !message.trim()}
      class="send-button"
    >
      {#if loading}
        ⏳
      {:else}
        ➤
      {/if}
    </button>
  </div>
</div>

<style>
  .message-composer {
    padding: 1rem;
    border-top: 1px solid var(--surface-stroke);
  }

  .input-container {
    display: flex;
    align-items: flex-end;
    gap: 0.5rem;
    background: var(--surface-glass);
    border-radius: var(--enforce-pill);
    padding: 0.5rem;
  }

  textarea {
    flex: 1;
    border: none;
    background: transparent;
    resize: none;
    outline: none;
    font-family: inherit;
    font-size: 0.9rem;
    line-height: 1.4;
    padding: 0.5rem;
    max-height: 120px;
    min-height: 20px;
  }

  .send-button {
    background: var(--accent);
    color: white;
    border: none;
    border-radius: 50%;
    width: 32px;
    height: 32px;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    font-size: 1rem;
    transition: all var(--transition-smooth);
  }

  .send-button:hover:not(:disabled) {
    background: var(--accent-secondary);
    transform: scale(1.05);
  }

  .send-button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
    transform: none;
  }
</style>
