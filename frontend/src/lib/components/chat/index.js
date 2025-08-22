// Export main chat component
export { default as Chat } from './Chat.svelte';

// Export chat stores
export { chatStore, messageHistory, canSendMessage, isStreaming } from './stores/chatStore';

// Export chat service - temporarily disabled due to build issues
// export { chatService } from './services/chatService';
