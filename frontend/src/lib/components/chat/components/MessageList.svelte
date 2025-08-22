<script lang="ts">
  export let messages: any[] = [];
  export let isCarousel = false;
</script>

<div class="message-list" class:carousel={isCarousel}>
  {#each messages as message (message.id)}
    <div
      class="message"
      class:user={message.type === 'user'}
      class:assistant={message.type === 'assistant'}
    >
      <div class="message-content">
        <div class="message-header">
          <span class="message-name">{message.name}</span>
        </div>
        <div class="message-text">{message.message}</div>
      </div>
    </div>
  {/each}
</div>

<style>
  .message-list {
    flex: 1;
    overflow-y: auto;
    padding: 1rem;
    display: flex;
    flex-direction: column;
    gap: 1rem;
  }

  .message {
    display: flex;
    max-width: 80%;
  }

  .message.user {
    align-self: flex-end;
  }

  .message.assistant {
    align-self: flex-start;
  }

  .message-content {
    background: var(--surface-glass);
    border-radius: var(--enforce-pill);
    padding: 0.75rem 1rem;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
  }

  .message.user .message-content {
    background: var(--accent);
    color: white;
  }

  .message-header {
    font-size: 0.8rem;
    font-weight: 600;
    margin-bottom: 0.25rem;
    opacity: 0.8;
  }

  .message-text {
    font-size: 0.9rem;
    line-height: 1.4;
    white-space: pre-wrap;
  }
</style>
