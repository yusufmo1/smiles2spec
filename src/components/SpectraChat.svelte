<script>
  import { Chat, MessageList, Message, Input } from '@ai-sdk/svelte';
  import { chatMessages, latestSpectrum, latestStructure, latestSmiles } from '../stores.js';
  import { chat } from '../services/api.js';
  import { get } from 'svelte/store';
  import Plotly from 'plotly.js-dist-min';
  
  let isProcessing = false;
  
  async function handleSend(userText) {
    if (!userText.trim() || isProcessing) return;
    
    isProcessing = true;
    try {
      // Add user message to the chat
      chatMessages.update(m => [...m, { role: 'user', content: userText }]);
      
      // Set a temporary message to show typing indicator
      chatMessages.update(m => [...m, { role: 'assistant', content: '', typing: true }]);
      
      // Get chat response
      const res = await chat(get(chatMessages).filter(m => !m.typing));
      
      // Update assistant message with streamed response
      let content = '';
      
      // Replace the temporary typing message with actual content
      chatMessages.update(m => {
        const messages = [...m];
        // Remove the typing indicator
        messages.pop();
        // Add real assistant message
        messages.push({ role: 'assistant', content: '' });
        return messages;
      });
      
      // Stream chunks into the messages store
      for await (const chunk of res) {
        content += chunk;
        chatMessages.update(m => {
          const messages = [...m];
          const last = messages.at(-1);
          if (last?.role === 'assistant') {
            last.content = content;
          }
          return messages;
        });
      }
    } catch (error) {
      console.error('Chat error:', error);
      // Show error in chat
      chatMessages.update(m => {
        const messages = [...m];
        // Remove typing indicator if it exists
        if (messages.at(-1)?.typing) {
          messages.pop();
        }
        messages.push({ 
          role: 'assistant', 
          content: `Sorry, I encountered an error: ${error.message}. Please try again.` 
        });
        return messages;
      });
    } finally {
      isProcessing = false;
    }
  }
  
  async function attachSpectrum() {
    try {
      const plotElement = document.querySelector('.plot-container');
      if (!plotElement) {
        throw new Error("No spectrum plot found");
      }
      
      const img = await Plotly.toImage(plotElement, { 
        format: 'png',
        width: 800,
        height: 500
      });
      
      const smiles = get(latestSmiles);
      const message = `Here is the spectrum ${smiles ? `for ${smiles}` : "I mentioned"}.`;
      
      chatMessages.update(m => [
        ...m, 
        { 
          role: 'user', 
          content: message,
          files: [{ name: 'spectrum.png', data: img }] 
        }
      ]);
      
      // Immediately follow up with assistant reply
      handleSend(message);
    } catch (error) {
      console.error('Error attaching spectrum:', error);
      alert(`Could not attach spectrum: ${error.message}`);
    }
  }
</script>

<div class="spectra-chat-container">
  <Chat class="glass-card chat-shell">
    <MessageList messages={$chatMessages}>
      <Message let:message class="chat-bubble {message.role}">
        {#if message.typing}
          <div class="typing-indicator">
            <span></span><span></span><span></span>
          </div>
        {:else if message.files && message.files.length > 0}
          <div class="message-content">
            <p>{message.content}</p>
            <div class="file-attachments">
              {#each message.files as file}
                <img src={file.data} alt="Attached spectrum" class="attached-image" />
              {/each}
            </div>
          </div>
        {:else}
          <div class="message-content">
            {message.content}
          </div>
        {/if}
      </Message>
    </MessageList>
    
    <div class="chat-actions">
      {#if $latestSpectrum}
        <button class="pill-button" on:click={attachSpectrum}>
          <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
            <path d="M12 4V20M20 12H4" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
          </svg>
          Attach spectrum
        </button>
      {/if}
    </div>
    
    <Input 
      placeholder="Ask anything about mass spectra..." 
      on:submit={e => handleSend(e.detail)}
      disabled={isProcessing}
    />
  </Chat>
</div>

<style>
  .spectra-chat-container {
    height: 100%;
    display: flex;
    flex-direction: column;
  }
  
  .chat-shell {
    display: flex;
    flex-direction: column;
    height: 100%;
    overflow: hidden;
    padding: 1rem;
  }
  
  .chat-actions {
    display: flex;
    justify-content: center;
    margin: 0.5rem 0;
  }
  
  .pill-button {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    background: var(--accent-soft);
    color: var(--accent);
    border: none;
    font-size: 0.8rem;
    font-weight: 500;
    padding: 0.5rem 1rem;
    border-radius: var(--enforce-pill, 999px);
    cursor: pointer;
    transition: all 0.2s;
  }
  
  .pill-button:hover {
    background: var(--accent);
    color: white;
  }
  
  .attached-image {
    max-width: 100%;
    border-radius: 4px;
    margin-top: 0.5rem;
  }
  
  /* Typing indicator animation */
  .typing-indicator {
    display: flex;
    align-items: center;
    gap: 4px;
    padding: 4px 8px;
  }
  
  .typing-indicator span {
    width: 6px;
    height: 6px;
    background-color: var(--accent);
    border-radius: 50%;
    display: inline-block;
    animation: bounce 1.5s infinite ease-in-out;
  }
  
  .typing-indicator span:nth-child(1) { animation-delay: 0s; }
  .typing-indicator span:nth-child(2) { animation-delay: 0.2s; }
  .typing-indicator span:nth-child(3) { animation-delay: 0.4s; }
  
  @keyframes bounce {
    0%, 60%, 100% { transform: translateY(0); }
    30% { transform: translateY(-6px); }
  }
  
  /* Message styling */
  .chat-bubble {
    margin: 4px 0;
    padding: 8px 12px;
    border-radius: 12px;
    max-width: 80%;
    overflow-wrap: break-word;
  }
  
  .chat-bubble.user {
    align-self: flex-end;
    background: var(--accent-soft);
    color: var(--accent-text);
  }
  
  .chat-bubble.assistant {
    align-self: flex-start;
    background: var(--soft-bg);
    color: var(--text);
  }
  
  .file-attachments {
    margin-top: 8px;
  }
</style> 