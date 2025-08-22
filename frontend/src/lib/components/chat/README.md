# Chat System

The chat system provides AI-powered conversations about mass spectrometry and molecular analysis, with context-aware responses based on current spectrum data.

## Architecture Overview

```
chat/
├── Chat.svelte               # Main chat interface
├── index.js                  # Exports
├── components/               # Sub-components
│   ├── MessageComposer.svelte    # Message input interface
│   ├── MessageList.svelte        # Message display list
│   └── PlaceholderState.svelte   # Empty state UI
├── services/                 # Business logic
│   └── chatService.ts        # Chat API integration
└── stores/                   # State management
    └── chatStore.ts          # Chat state management
```

## Core Component

### Chat.svelte

The main chat interface that orchestrates the conversation experience:

```svelte
<script>
  import { Chat } from '$lib/components/chat';
</script>

<Chat
  smiles="CCO"
  spectrumData={currentSpectrum}
  placeholder="Ask about this spectrum..."
  on:messageAdded={handleNewMessage}
  on:error={handleError}
/>
```

**Features:**

- Real-time messaging with AI
- Streaming response support
- Spectrum context integration
- Message history persistence
- Error handling and retries
- Markdown rendering
- Loading states

**Props:**

- `smiles` (string) - Current molecule for context
- `spectrumData` (object) - Spectrum data for analysis
- `placeholder` (string) - Input placeholder text
- `maxMessages` (number) - Message history limit
- `enableRetry` (boolean) - Allow message retry on failure

**Events:**

- `messageAdded` - New message added to conversation
- `responseReceived` - AI response completed
- `error` - Error occurred during communication
- `streamChunk` - Real-time streaming chunk received

## Sub-Components

### MessageComposer.svelte

Message input interface with enhanced functionality:

```svelte
<MessageComposer
  placeholder="Type your message..."
  disabled={isLoading}
  on:submit={handleSubmit}
  on:clear={handleClear}
/>
```

**Features:**

- Auto-resizing textarea
- Send button state management
- Keyboard shortcuts (Ctrl+Enter)
- Input validation
- Character count (optional)
- Markdown preview (future)

**Keyboard Shortcuts:**

- `Ctrl/Cmd + Enter` - Send message
- `Escape` - Clear input
- `↑/↓` - Navigate message history

### MessageList.svelte

Displays conversation history with proper formatting:

```svelte
<MessageList
  messages={$chatStore.messages}
  isLoading={$chatStore.isLoading}
  showTimestamps={true}
  on:retry={retryMessage}
/>
```

**Features:**

- Automatic scrolling to latest message
- Message grouping by sender
- Timestamp display
- Retry failed messages
- Copy message content
- Message status indicators
- Markdown rendering with syntax highlighting

**Message Types:**

- **User messages**: User input with timestamp
- **AI responses**: Streaming responses with typing indicators
- **System messages**: Status updates and errors
- **Context messages**: Automatic spectrum context

### PlaceholderState.svelte

Empty state when no messages exist:

```svelte
<PlaceholderState
  title="Start a conversation"
  description="Ask questions about mass spectrometry..."
  suggestions={exampleQuestions}
  on:suggestionClick={handleSuggestion}
/>
```

**Features:**

- Welcome message
- Suggested questions
- Getting started tips
- Contextual help based on current data

## Services

### chatService.ts

Core chat API integration and business logic:

```typescript
import { chatService } from '$lib/components/chat/services';

// Send message with spectrum context
const response = await chatService.sendMessage({
  message: 'Explain this spectrum',
  smiles: 'CCO',
  spectrumData: currentSpectrum,
  stream: true,
  onChunk: (chunk) => console.log(chunk),
});

// Get message history
const history = await chatService.getMessageHistory();

// Clear conversation
await chatService.clearHistory();
```

**Core Methods:**

#### `sendMessage(options: ChatMessageOptions): Promise<ChatResponse>`

Sends a message to the AI service with optional context.

**Options:**

- `message` (string) - User message text
- `smiles` (string, optional) - Current molecule SMILES
- `spectrumData` (object, optional) - Current spectrum data
- `stream` (boolean) - Enable streaming responses
- `onChunk` (function) - Streaming chunk handler
- `temperature` (number) - Response creativity (0-1)

#### `getMessageHistory(): Promise<ChatMessage[]>`

Retrieves stored conversation history.

#### `clearHistory(): Promise<void>`

Clears all stored messages.

#### `retryMessage(messageId: string): Promise<ChatResponse>`

Retries a failed message with the same parameters.

#### `generateSuggestions(context?: SpectrumContext): string[]`

Generates contextual question suggestions.

### API Integration

The service integrates with the backend chat API:

```typescript
// Streaming request
async function sendStreamingMessage(options: ChatMessageOptions) {
  const response = await fetch('/api/chat', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      Accept: 'text/event-stream',
    },
    body: JSON.stringify({
      messages: buildMessageHistory(options.message),
      smiles: options.smiles,
      stream: true,
    }),
  });

  // Handle Server-Sent Events
  const reader = response.body?.getReader();
  const decoder = new TextDecoder();

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    const chunk = decoder.decode(value);
    const lines = chunk.split('\n');

    for (const line of lines) {
      if (line.startsWith('data: ')) {
        const data = line.slice(6);
        if (data === '[DONE]') {
          return;
        }

        try {
          const parsed = JSON.parse(data);
          options.onChunk?.(parsed.chunk);
        } catch (e) {
          console.warn('Failed to parse chunk:', data);
        }
      }
    }
  }
}
```

## State Management

### chatStore.ts

Manages chat state, message history, and UI state:

```typescript
import { chatStore } from '$lib/components/chat/stores';

// Current state
$: messages = $chatStore.messages;
$: isLoading = $chatStore.isLoading;
$: currentInput = $chatStore.currentInput;

// Actions
chatStore.addMessage({
  id: generateId(),
  role: 'user',
  content: 'Hello!',
  timestamp: Date.now(),
});

chatStore.setLoading(true);
chatStore.updateInput('New message...');
chatStore.clearMessages();
```

**State Properties:**

- `messages` - Array of chat messages
- `isLoading` - AI response loading state
- `currentInput` - Current input text
- `error` - Current error state
- `streamingMessageId` - ID of message being streamed
- `retryCount` - Number of retry attempts
- `lastActivity` - Timestamp of last activity

**Message Structure:**

```typescript
interface ChatMessage {
  id: string;
  role: 'user' | 'assistant' | 'system';
  content: string;
  timestamp: number;
  context?: {
    smiles?: string;
    spectrumData?: SpectrumData;
  };
  status?: 'sending' | 'sent' | 'failed' | 'streaming' | 'complete';
  metadata?: {
    tokens?: number;
    model?: string;
    processingTime?: number;
  };
}
```

**Actions:**

- `addMessage()` - Add new message to conversation
- `updateMessage()` - Update existing message content
- `setLoading()` - Set loading state
- `setError()` - Set error state
- `updateInput()` - Update input field text
- `clearMessages()` - Clear all messages
- `removeMessage()` - Remove specific message

## Context Integration

### Spectrum Context

The chat system automatically includes spectrum context in conversations:

```typescript
// Automatic context building
function buildSpectrumContext(smiles: string, spectrumData: any) {
  return {
    molecule: {
      smiles,
      molecularWeight: spectrumData.molecularWeight,
      formula: spectrumData.formula,
    },
    spectrum: {
      peakCount: spectrumData.peaks.length,
      basePeak: spectrumData.peaks[0],
      molecularIon: spectrumData.molecularIon,
      fragments: spectrumData.majorFragments,
    },
  };
}

// Context injection in messages
const contextualMessage = {
  role: 'user',
  content: userMessage,
  context: buildSpectrumContext(currentSmiles, currentSpectrum),
};
```

### Contextual Suggestions

The system provides contextual question suggestions:

```typescript
// Context-aware suggestions
function generateSuggestions(context: SpectrumContext): string[] {
  const suggestions = [
    'What is the molecular ion peak?',
    'Explain the fragmentation pattern',
    'What are the major fragments?',
    'How accurate is this prediction?',
  ];

  if (context.spectrum.peakCount > 20) {
    suggestions.push('Why are there so many peaks?');
  }

  if (context.molecule.molecularWeight > 500) {
    suggestions.push('Is this a large molecule?');
  }

  return suggestions;
}
```

## Streaming Responses

### Real-time Updates

The chat system supports streaming responses for better UX:

```typescript
// Streaming message handler
function handleStreamingResponse(messageId: string) {
  let accumulatedContent = '';

  return {
    onChunk: (chunk: string) => {
      accumulatedContent += chunk;
      chatStore.updateMessage(messageId, {
        content: accumulatedContent,
        status: 'streaming',
      });
    },
    onComplete: () => {
      chatStore.updateMessage(messageId, {
        status: 'complete',
      });
    },
    onError: (error: Error) => {
      chatStore.updateMessage(messageId, {
        status: 'failed',
        error: error.message,
      });
    },
  };
}
```

### Typing Indicators

Visual feedback during response generation:

```svelte
<!-- In MessageList.svelte -->
{#if $chatStore.isLoading}
  <div class="typing-indicator">
    <div class="typing-dots">
      <span></span>
      <span></span>
      <span></span>
    </div>
    <span class="typing-text">AI is thinking...</span>
  </div>
{/if}
```

## Error Handling

### Retry Logic

Robust error handling with retry capabilities:

```typescript
async function sendMessageWithRetry(options: ChatMessageOptions, maxRetries = 3) {
  let retries = 0;

  while (retries < maxRetries) {
    try {
      return await chatService.sendMessage(options);
    } catch (error) {
      retries++;

      if (error.type === 'NETWORK_ERROR' && retries < maxRetries) {
        // Exponential backoff
        await delay(Math.pow(2, retries) * 1000);
        continue;
      }

      throw error;
    }
  }
}
```

### Error Types

```typescript
interface ChatError {
  type: 'NETWORK_ERROR' | 'API_ERROR' | 'VALIDATION_ERROR' | 'RATE_LIMIT';
  message: string;
  retryable: boolean;
  details?: any;
}

// Error handling in UI
function handleChatError(error: ChatError) {
  switch (error.type) {
    case 'NETWORK_ERROR':
      showRetryOption(error);
      break;
    case 'RATE_LIMIT':
      showRateLimitMessage(error);
      break;
    case 'API_ERROR':
      showApiErrorMessage(error);
      break;
    default:
      showGenericError(error);
  }
}
```

## Performance Optimization

### Message Virtualization

For long conversations, implement message virtualization:

```typescript
// Virtual scrolling for large message lists
import { VirtualList } from '@sveltejs/virtual-list';

// Only render visible messages
const ITEM_HEIGHT = 80;
const visibleMessages = messages.slice(startIndex, endIndex);
```

### Debounced Input

Prevent excessive API calls during typing:

```typescript
const debouncedSend = debounce(async (message: string) => {
  await chatService.sendMessage({ message });
}, 1000);

// Handle input changes
function handleInputChange(value: string) {
  chatStore.updateInput(value);

  // Optional: Send typing indicators
  if (value.length > 3) {
    debouncedSend(value);
  }
}
```

### Message Persistence

Efficient storage of conversation history:

```typescript
// Local storage with size limits
class ChatPersistence {
  private static readonly MAX_MESSAGES = 100;
  private static readonly STORAGE_KEY = 'chat_history';

  static saveMessages(messages: ChatMessage[]) {
    const limitedMessages = messages.slice(-this.MAX_MESSAGES);
    localStorage.setItem(this.STORAGE_KEY, JSON.stringify(limitedMessages));
  }

  static loadMessages(): ChatMessage[] {
    const stored = localStorage.getItem(this.STORAGE_KEY);
    return stored ? JSON.parse(stored) : [];
  }

  static clearMessages() {
    localStorage.removeItem(this.STORAGE_KEY);
  }
}
```

## Accessibility

### Keyboard Navigation

Complete keyboard support for chat interface:

```typescript
// Message list navigation
const messageListKeyHandler = (e: KeyboardEvent) => {
  switch (e.key) {
    case 'ArrowUp':
      navigateToPreviousMessage();
      break;
    case 'ArrowDown':
      navigateToNextMessage();
      break;
    case 'Home':
      scrollToFirstMessage();
      break;
    case 'End':
      scrollToLastMessage();
      break;
  }
};

// Input field shortcuts
const inputKeyHandler = (e: KeyboardEvent) => {
  if (e.ctrlKey || e.metaKey) {
    switch (e.key) {
      case 'Enter':
        sendMessage();
        break;
      case 'k':
        clearMessages();
        break;
    }
  }
};
```

### Screen Reader Support

Proper ARIA labels and live regions:

```svelte
<!-- Message list with proper semantics -->
<ul role="log" aria-live="polite" aria-label="Chat conversation">
  {#each messages as message}
    <li role="article" aria-label="{message.role} message">
      <div class="message-content" id="message-{message.id}">
        {message.content}
      </div>
      <time datetime={message.timestamp} class="sr-only">
        {formatTimestamp(message.timestamp)}
      </time>
    </li>
  {/each}
</ul>

<!-- Input with proper labeling -->
<label for="chat-input" class="sr-only"> Type your message about the spectrum </label>
<textarea
  id="chat-input"
  aria-describedby="chat-help"
  aria-expanded={showSuggestions}
  bind:value={inputValue}
>
</textarea>
```

### Focus Management

Proper focus handling for dynamic content:

```typescript
// Focus new messages for screen readers
function announceNewMessage(message: ChatMessage) {
  const messageElement = document.getElementById(`message-${message.id}`);
  if (messageElement) {
    messageElement.focus();
    messageElement.scrollIntoView({ behavior: 'smooth' });
  }
}

// Manage focus during state changes
function handleLoadingStateChange(isLoading: boolean) {
  if (isLoading) {
    // Announce loading state
    announceToScreenReader('AI is generating response');
  } else {
    // Focus on new response
    announceNewMessage(lastMessage);
  }
}
```

## Usage Examples

### Basic Chat Implementation

```svelte
<script>
  import { Chat } from '$lib/components/chat';
  import { appState } from '$lib/stores';

  $: currentSpectrum = $appState.currentPredictionData;
  $: currentSmiles = $appState.currentSmiles;

  function handleNewMessage(event) {
    const { message } = event.detail;
    console.log('New message:', message);
  }

  function handleError(event) {
    const { error } = event.detail;
    console.error('Chat error:', error);
  }
</script>

<Chat
  smiles={currentSmiles}
  spectrumData={currentSpectrum}
  placeholder="Ask about this molecule..."
  on:messageAdded={handleNewMessage}
  on:error={handleError}
/>
```

### Advanced Configuration

```svelte
<Chat
  smiles={currentSmiles}
  spectrumData={currentSpectrum}
  maxMessages={50}
  enableRetry={true}
  streamingEnabled={true}
  showTimestamps={true}
  placeholder="What would you like to know about this spectrum?"
  suggestions={[
    'Explain the fragmentation pattern',
    'What is the base peak?',
    'How was this spectrum predicted?',
  ]}
  on:messageAdded={handleNewMessage}
  on:responseReceived={handleResponse}
  on:streamChunk={handleStreamChunk}
  on:error={handleError}
/>
```

### Custom Message Handling

```svelte
<script>
  import { chatStore } from '$lib/components/chat/stores';
  import { chatService } from '$lib/components/chat/services';

  async function sendCustomMessage(content: string) {
    // Add user message
    const userMessage = {
      id: generateId(),
      role: 'user',
      content,
      timestamp: Date.now(),
    };

    chatStore.addMessage(userMessage);
    chatStore.setLoading(true);

    try {
      // Send with context
      const response = await chatService.sendMessage({
        message: content,
        smiles: $appState.currentSmiles,
        spectrumData: $appState.currentPredictionData,
        stream: true,
        onChunk: (chunk) => {
          // Handle streaming chunks
          updateStreamingMessage(chunk);
        },
      });

      // Add complete response
      chatStore.addMessage({
        id: generateId(),
        role: 'assistant',
        content: response.content,
        timestamp: Date.now(),
        metadata: {
          tokens: response.tokens,
          processingTime: response.processingTime,
        },
      });
    } catch (error) {
      chatStore.setError(error.message);
    } finally {
      chatStore.setLoading(false);
    }
  }
</script>
```

## Testing

### Component Testing

```javascript
import { render, fireEvent, waitFor } from '@testing-library/svelte';
import Chat from './Chat.svelte';

describe('Chat Component', () => {
  test('sends message on submit', async () => {
    const { getByTestId, getByText } = render(Chat, {
      props: { smiles: 'CCO' },
    });

    const input = getByTestId('chat-input');
    const sendButton = getByTestId('send-button');

    await fireEvent.input(input, { target: { value: 'Test message' } });
    await fireEvent.click(sendButton);

    await waitFor(() => {
      expect(getByText('Test message')).toBeInTheDocument();
    });
  });

  test('handles streaming responses', async () => {
    const mockStream = vi.fn();
    const { component } = render(Chat);

    component.$on('streamChunk', mockStream);

    // Simulate streaming response
    await fireEvent.custom(component, 'streamChunk', {
      detail: { chunk: 'Streaming...' },
    });

    expect(mockStream).toHaveBeenCalled();
  });
});
```

### Service Testing

```javascript
import { chatService } from './services/chatService';

describe('Chat Service', () => {
  test('sends message with context', async () => {
    const response = await chatService.sendMessage({
      message: 'Test',
      smiles: 'CCO',
      spectrumData: mockSpectrumData,
    });

    expect(response).toHaveProperty('content');
    expect(response.content).toBeTruthy();
  });

  test('handles streaming responses', async () => {
    const chunks = [];

    await chatService.sendMessage({
      message: 'Test streaming',
      stream: true,
      onChunk: (chunk) => chunks.push(chunk),
    });

    expect(chunks.length).toBeGreaterThan(0);
  });
});
```

## Best Practices

### Message Handling

1. **Immediate UI updates** for user messages
2. **Progressive enhancement** for streaming responses
3. **Graceful degradation** for failed messages
4. **Context preservation** across conversations

### Performance

1. **Virtual scrolling** for long conversations
2. **Debounced inputs** to prevent spam
3. **Efficient re-renders** using message IDs
4. **Lazy loading** of message history

### User Experience

1. **Clear loading states** during responses
2. **Retry mechanisms** for failed messages
3. **Contextual suggestions** to guide users
4. **Accessibility compliance** for all interactions

### Error Handling

1. **Graceful error messages** for users
2. **Detailed logging** for developers
3. **Retry logic** for transient errors
4. **Offline support** for network issues

## Future Enhancements

- **Voice Input**: Speech-to-text integration
- **Message Search**: Search through conversation history
- **Export Conversations**: Save conversations to file
- **Message Threading**: Organize related messages
- **Collaborative Chat**: Multi-user conversations
- **Custom Prompts**: User-defined question templates
- **Integration APIs**: Connect with external chemistry tools
