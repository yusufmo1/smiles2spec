<script>
  import { onMount } from 'svelte';
  import { writable } from 'svelte/store';
  import { chatWithSpectrum } from '../services/api.js';
  
  // Get access to the hasFirstPrediction state
  export let hasSmilesPrediction = false;
  
  // Chat message store
  const messages = writable([
    { 
      id: 1, 
      avatar: "/assets/images/spectra-avatar.png",
      name: "Spectrum", 
      message: "Hello! I'm Spectrum, your spectral analysis assistant. Ask me about molecular structures, SMILES notation, or help interpreting mass spectra.", 
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
          name: "Spectrum",
          message: "Thinking...",
          timestamp: new Date().toISOString(),
          type: "assistant",
          thinking: true
        }
      ]);
      
      // Use our API service to get Spectrum's response
      const result = await chatWithSpectrum(
        buildMessageHistory($messages.filter(m => !m.thinking))
      );
      
      // Remove the thinking message
      messages.update(msgs => msgs.filter(m => !m.thinking));
      
      // Add the bot's response
      const botMsg = {
        id: $messages.length + 1,
        avatar: "/assets/images/spectra-avatar.png",
        name: "Spectrum",
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
          name: "Spectrum",
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
  {#if hasSmilesPrediction}
    <div class="chat-messages" bind:this={chatElement}>
      {#each $messages as message}
        <div class="message-bubble {message.type === 'user' ? 'user-message' : 'spectrum-message'}">
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
        placeholder="Ask Spectrum a question..."
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
  {:else}
    <div class="placeholder">
      <span class="placeholder-icon">💬</span>
      <p>Chat available after SMILES processing</p>
    </div>
  {/if}
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
  
  .spectrum-message {
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
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  }
  
  .message-avatar img {
    width: 100%;
    height: 100%;
    object-fit: cover;
  }
  
  .message-content {
    background: var(--glass-card-surface);
    border-radius: var(--enforce-pill);
    padding: 0.85rem 1.25rem;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
    position: relative;
    backdrop-filter: blur(8px);
    border: 1px solid rgba(255, 255, 255, 0.2);
  }
  
  .user-message .message-content {
    background: linear-gradient(135deg, var(--accent-soft) 0%, rgba(120, 121, 255, 0.1) 100%);
    border: 1px solid rgba(120, 121, 255, 0.2);
  }
  
  .message-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 0.5rem;
  }
  
  .message-name {
    font-weight: 600;
    font-size: 0.85rem;
    color: var(--accent);
  }
  
  .user-message .message-name {
    background: linear-gradient(135deg, var(--accent), var(--accent-secondary));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
  }
  
  .message-time {
    font-size: 0.75rem;
    color: var(--text-secondary);
  }
  
  .message-text {
    font-size: 0.9rem;
    line-height: 1.5;
    color: var(--text-primary);
    word-wrap: break-word;
  }
  
  .input-area {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    border-top: 1px solid var(--surface-stroke);
    padding-top: 0.75rem;
  }
  
  textarea {
    flex: 1;
    border: 1px solid var(--surface-stroke);
    border-radius: var(--enforce-pill);
    padding: 0.85rem 1.25rem;
    resize: none;
    background: var(--glass-card-surface);
    color: var(--text-primary);
    font-family: inherit;
    font-size: 0.9rem;
    line-height: 1.5;
    transition: all 0.2s var(--transition-smooth);
    height: auto;
    min-height: 42px;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.03);
  }
  
  textarea:focus {
    outline: none;
    border-color: var(--accent);
    box-shadow: 0 3px 12px rgba(120, 121, 255, 0.15);
  }
  
  .send-button {
    width: 44px;
    height: 44px;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-shrink: 0;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    color: white;
    border-radius: 50%;
    transition: all 0.2s var(--transition-smooth);
    border: none;
    box-shadow: 0 4px 12px rgba(120, 121, 255, 0.3);
  }
  
  .send-button:hover {
    transform: translateY(-2px);
    box-shadow: 0 6px 18px rgba(120, 121, 255, 0.4);
  }
  
  .send-button:active {
    transform: translateY(0px);
  }
  
  .send-button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
  }
  
  .loading-indicator {
    width: 20px;
    height: 20px;
    border: 2px solid rgba(255, 255, 255, 0.3);
    border-top-color: white;
    border-radius: 50%;
    animation: spin 0.8s linear infinite;
  }
  
  @keyframes spin {
    to { transform: rotate(360deg); }
  }
  
  /* Placeholder styles */
  .placeholder {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    height: 100%;
    width: 100%;
    background: rgba(0, 0, 0, 0.05);
    border-radius: 8px;
    color: rgba(0, 0, 0, 0.4);
  }
  
  .placeholder-icon {
    font-size: 3rem;
    margin-bottom: 1rem;
  }
  
  .placeholder p {
    font-size: 0.9rem;
    font-weight: 500;
  }
</style> 