<script>
  import { onMount } from 'svelte';
  import { writable } from 'svelte/store';
  import { chatWithSpectra } from '../services/api.js';
  
  // Chat message store
  const messages = writable([
    { 
      id: 1, 
      avatar: "/assets/images/spectra-avatar.png",
      name: "Spectra", 
      message: "Hello! I'm Spectra, your spectral analysis assistant. Ask me about molecular structures, SMILES notation, or help interpreting mass spectra.", 
      timestamp: new Date().toISOString(),
      type: "assistant"
    }
  ]);
  
  let userMessage = '';
  let chatElement;
  let loading = false;
  
  // Scroll to bottom when messages update
  $: if (chatElement && $messages) {
    setTimeout(() => {
      chatElement.scrollTo({ top: chatElement.scrollHeight, behavior: 'smooth' });
    }, 50);
  }
  
  async function handleSend() {
    if (!userMessage.trim() || loading) return;
    
    // Add user message
    const userMsg = {
      id: $messages.length + 1,
      avatar: "/assets/images/user-avatar.png",
      name: "You",
      message: userMessage,
      timestamp: new Date().toISOString(),
      type: "user"
    };
    
    messages.update(msgs => [...msgs, userMsg]);
    
    // Clear input
    const tempMessage = userMessage;
    userMessage = '';
    loading = true;
    
    try {
      // Show "thinking" message
      messages.update(msgs => [
        ...msgs,
        {
          id: $messages.length + 1,
          avatar: "/assets/images/spectra-avatar.png",
          name: "Spectra",
          message: "Thinking...",
          timestamp: new Date().toISOString(),
          type: "assistant",
          thinking: true
        }
      ]);
      
      // Use our API service to get Spectra's response
      const result = await chatWithSpectra(
        buildMessageHistory($messages.filter(m => !m.thinking))
      );
      
      // Remove the thinking message
      messages.update(msgs => msgs.filter(m => !m.thinking));
      
      // Add the bot's response
      const botMsg = {
        id: $messages.length + 1,
        avatar: "/assets/images/spectra-avatar.png",
        name: "Spectra",
        message: result.message,
        timestamp: new Date().toISOString(),
        type: "assistant"
      };
      
      messages.update(msgs => [...msgs, botMsg]);
    } catch (error) {
      console.error('Error:', error);
      messages.update(msgs => [
        ...msgs.filter(m => !m.thinking),
        {
          id: $messages.length + 1,
          avatar: "/assets/images/spectra-avatar.png",
          name: "Spectra",
          message: "Sorry, I encountered an error. Please try again.",
          timestamp: new Date().toISOString(),
          type: "assistant"
        }
      ]);
    } finally {
      loading = false;
    }
  }
  
  function buildMessageHistory(messages) {
    return messages.map(msg => ({
      role: msg.type,
      content: msg.message
    }));
  }
  
  function handleKeyDown(event) {
    if (event.key === 'Enter' && !event.shiftKey) {
      event.preventDefault();
      handleSend();
    }
  }
</script>

<div class="chat-container">
  <div class="chat-messages" bind:this={chatElement}>
    {#each $messages as message}
      <div class="message-bubble {message.type === 'user' ? 'user-message' : 'spectra-message'}">
        <div class="message-avatar">
          <img src={message.avatar} alt={message.name} />
        </div>
        <div class="message-content">
          <div class="message-header">
            <span class="message-name">{message.name}</span>
            <span class="message-time">{new Date(message.timestamp).toLocaleTimeString()}</span>
          </div>
          <div class="message-text">{message.message}</div>
        </div>
      </div>
    {/each}
  </div>
  
  <div class="input-area">
    <textarea 
      bind:value={userMessage}
      on:keydown={handleKeyDown}
      placeholder="Ask Spectra a question..."
      disabled={loading}
      rows="2"
    ></textarea>
    <button 
      class="pill-button send-button"
      on:click={handleSend}
      disabled={!userMessage.trim() || loading}
    >
      {#if loading}
        <span class="loading-indicator"></span>
      {:else}
        <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
          <path d="M2.01 21L23 12 2.01 3 2 10l15 2-15 2z"/>
        </svg>
      {/if}
    </button>
  </div>
</div>

<style>
  .chat-container {
    display: flex;
    flex-direction: column;
    height: 100%;
    width: 100%;
  }
  
  .chat-messages {
    flex: 1;
    overflow-y: auto;
    padding-right: 0.5rem;
    display: flex;
    flex-direction: column;
    gap: 0.75rem;
    margin-bottom: 1rem;
    min-height: 200px;
    max-height: 300px;
  }
  
  .message-bubble {
    display: flex;
    gap: 0.75rem;
    max-width: 80%;
  }
  
  .spectra-message {
    align-self: flex-start;
  }
  
  .user-message {
    align-self: flex-end;
    flex-direction: row-reverse;
  }
  
  .message-avatar {
    flex-shrink: 0;
    width: 28px;
    height: 28px;
    border-radius: 50%;
    overflow: hidden;
    background: var(--accent-soft);
  }
  
  .message-avatar img {
    width: 100%;
    height: 100%;
    object-fit: cover;
  }
  
  .message-content {
    background: rgba(255, 255, 255, 0.8);
    border-radius: var(--enforce-pill);
    padding: 0.5rem 0.75rem;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
  }
  
  .user-message .message-content {
    background: var(--accent-soft);
  }
  
  .message-header {
    display: flex;
    justify-content: space-between;
    margin-bottom: 0.25rem;
  }
  
  .message-name {
    font-weight: 600;
    color: var(--text-primary);
    font-size: 0.8rem;
  }
  
  .message-time {
    font-size: 0.7rem;
    color: var(--text-tertiary);
  }
  
  .message-text {
    color: var(--text-secondary);
    white-space: pre-wrap;
    word-break: break-word;
    font-size: 0.9rem;
  }
  
  .input-area {
    display: flex;
    gap: 0.75rem;
    margin-top: auto;
    padding-top: 0.75rem;
    border-top: 1px solid var(--surface-stroke, rgba(0,0,0,0.08));
  }
  
  textarea {
    flex: 1;
    border: 1px solid var(--surface-stroke, rgba(0,0,0,0.08));
    border-radius: var(--enforce-pill);
    padding: 0.5rem 0.75rem;
    resize: none;
    background: rgba(255, 255, 255, 0.8);
    color: var(--text-primary);
    font-family: inherit;
    line-height: 1.5;
    font-size: 0.9rem;
    height: 2.5rem;
  }
  
  textarea:focus {
    outline: 2px solid var(--accent-soft);
    border-color: var(--accent);
  }
  
  .send-button {
    align-self: flex-end;
    padding: 0.5rem;
    background: var(--accent-soft);
    border: none;
    color: var(--accent);
    border-radius: var(--enforce-pill);
    cursor: pointer;
    display: flex;
    justify-content: center;
    align-items: center;
    transition: all 0.2s;
    height: 2.5rem;
    width: 2.5rem;
  }
  
  .send-button:hover:not(:disabled) {
    background: var(--accent);
    color: white;
  }
  
  .send-button:disabled {
    opacity: 0.4;
    cursor: not-allowed;
  }
  
  .loading-indicator {
    width: 16px;
    height: 16px;
    border: 2px solid rgba(255, 255, 255, 0.3);
    border-top: 2px solid white;
    border-radius: 50%;
    animation: spin 1s linear infinite;
  }
  
  @keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
  }
  
  /* Custom scrollbar */
  .chat-messages::-webkit-scrollbar {
    width: 4px;
  }
  
  .chat-messages::-webkit-scrollbar-thumb {
    background-color: rgba(0, 0, 0, 0.15);
    border-radius: 2px;
  }
  
  .chat-messages::-webkit-scrollbar-track {
    background-color: rgba(0, 0, 0, 0.05);
  }
  
  /* Responsive adjustments */
  @media (max-width: 768px) {
    .message-bubble {
      max-width: 90%;
    }
    
    .input-area {
      flex-direction: column;
    }
    
    .send-button {
      align-self: stretch;
    }
  }
</style> 