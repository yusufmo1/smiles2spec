<script lang="ts">
  import { onMount } from 'svelte';
  import { chatStore, messageHistory, canSendMessage } from './stores/chatStore';
  import { chatService } from './services/chatService';
  import MessageList from './components/MessageList.svelte';
  import MessageComposer from './components/MessageComposer.svelte';
  import PlaceholderState from './components/PlaceholderState.svelte';

  export let hasSmilesPrediction = false;
  export let currentSmiles: string | null = null;
  export let isCarousel = false;

  let prevSmiles: string | null = null;
  let messageListRef: any;

  // Update initial message when SMILES changes
  $: if (currentSmiles && hasSmilesPrediction && currentSmiles !== prevSmiles) {
    chatStore.updateInitialMessage(currentSmiles);
    prevSmiles = currentSmiles;
  }

  // Debug logging
  $: if (typeof window !== 'undefined') {
    console.log('Chat component - Current messages:', $chatStore.messages);
  }

  async function handleSend(event: CustomEvent) {
    const { content } = event.detail;

    if (!chatService.validateMessage(content)) {
      console.warn('Invalid message content:', content);
      return;
    }

    const sanitizedContent = chatService.sanitizeMessage(content);

    try {
      // Add user message
      chatStore.addUserMessage(sanitizedContent);
      chatStore.setLoading(true);

      // Start assistant message
      chatStore.startAssistantMessage();

      // Send message with streaming
      await chatService.sendMessage($messageHistory, {
        smiles: currentSmiles || undefined,
        onChunk: (chunk: string) => {
          console.log('Received streaming chunk:', chunk);
          chatStore.appendToStreamingMessage(chunk);
        },
      });
    } catch (error) {
      console.error('Chat error:', error);
      const errorMessage = chatService.formatErrorMessage(error);
      chatStore.handleError(errorMessage);
    } finally {
      chatStore.finishStreaming();
    }
  }

  // Reset chat when prediction status changes
  $: if (!hasSmilesPrediction) {
    chatStore.reset();
  }

  onMount(() => {
    console.log('Chat component mounted');
  });
</script>

{#if hasSmilesPrediction}
  <div class="chat-window" class:carousel={isCarousel}>
    <MessageList bind:this={messageListRef} messages={$chatStore.messages} {isCarousel} />
    <MessageComposer
      disabled={!$canSendMessage}
      loading={$chatStore.isLoading}
      on:send={handleSend}
    />
  </div>
{:else}
  <PlaceholderState />
{/if}

<style>
  .chat-window {
    height: 100%;
    width: 100%;
    display: flex;
    flex-direction: column;
    background: transparent; /* Let panel handle background */
    border-radius: 0; /* Let panel handle border radius */
    overflow: hidden;
    min-height: 0; /* Remove min-height constraint */
  }

  /* Carousel-specific styling */
  .chat-window.carousel {
    height: 100%;
    max-height: none; /* Remove height restrictions in carousel */
    min-height: 0; /* Allow shrinking in carousel */
    /* Add back some styling for carousel mode */
    background: rgba(255, 255, 255, 0.75);
    border-radius: var(--enforce-pill);
  }

  /* Accessibility improvements */
  @media (prefers-reduced-motion: reduce) {
    .chat-window :global(*) {
      animation-duration: 0.01ms !important;
      animation-iteration-count: 1 !important;
      transition-duration: 0.01ms !important;
    }
  }

  /* High contrast mode support */
  @media (prefers-contrast: high) {
    .chat-window.carousel {
      border: 2px solid currentColor;
      background: Canvas;
    }
  }
</style>
