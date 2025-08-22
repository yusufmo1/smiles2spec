import { writable, derived } from 'svelte/store';

interface ChatMessage {
  id: number;
  avatar: string;
  name: string;
  message: string;
  timestamp: string;
  type: 'user' | 'assistant';
}

interface ChatState {
  messages: ChatMessage[];
  isLoading: boolean;
  streamingMessageId: number | null;
  currentSmiles: string | null;
}

interface ApiChatMessage {
  role: 'user' | 'assistant';
  content: string;
}

const initialState: ChatState = {
  messages: [
    {
      id: 1,
      avatar: '/assets/images/spectra-avatar.png',
      name: 'Spectrum',
      message:
        "Hello! I'm Spectrum, your spectral-analysis assistant – ask me anything about SMILES or these spectra.",
      timestamp: new Date().toISOString(),
      type: 'assistant',
    },
  ],
  isLoading: false,
  streamingMessageId: null,
  currentSmiles: null,
};

function createChatStore() {
  const { subscribe, set, update } = writable(initialState);

  return {
    subscribe,

    // Message management
    addUserMessage: (content: string) =>
      update((state) => ({
        ...state,
        messages: [
          ...state.messages,
          {
            id: Date.now(),
            avatar: '/assets/images/user-avatar.png',
            name: 'You',
            message: content,
            timestamp: new Date().toISOString(),
            type: 'user',
          },
        ],
      })),

    startAssistantMessage: () =>
      update((state) => {
        const messageId = Date.now() + 1;
        return {
          ...state,
          streamingMessageId: messageId,
          messages: [
            ...state.messages,
            {
              id: messageId,
              avatar: '/assets/images/spectra-avatar.png',
              name: 'Spectrum',
              message: '',
              timestamp: new Date().toISOString(),
              type: 'assistant',
            },
          ],
        };
      }),

    appendToStreamingMessage: (chunk: string) =>
      update((state) => {
        if (!state.streamingMessageId) return state;

        return {
          ...state,
          messages: state.messages.map((msg) =>
            msg.id === state.streamingMessageId ? { ...msg, message: msg.message + chunk } : msg
          ),
        };
      }),

    updateStreamingMessage: (content: string) =>
      update((state) => {
        if (!state.streamingMessageId) return state;

        return {
          ...state,
          messages: state.messages.map((msg) =>
            msg.id === state.streamingMessageId ? { ...msg, message: content } : msg
          ),
        };
      }),

    finishStreaming: () =>
      update((state) => ({
        ...state,
        streamingMessageId: null,
        isLoading: false,
      })),

    setLoading: (loading: boolean) => update((state) => ({ ...state, isLoading: loading })),

    updateInitialMessage: (smiles: string) =>
      update((state) => {
        const newMessage = `Hello! I'm Spectrum, your spectral-analysis assistant. Ask me anything about the current molecule (${smiles}) or its mass spectrum!`;
        return {
          ...state,
          currentSmiles: smiles,
          messages: state.messages.map((msg, index) =>
            index === 0 ? { ...msg, message: newMessage } : msg
          ),
        };
      }),

    handleError: (errorMessage: string = 'Sorry – something went wrong. Try again?') =>
      update((state) => {
        if (!state.streamingMessageId) return state;

        return {
          ...state,
          messages: state.messages.map((msg) =>
            msg.id === state.streamingMessageId ? { ...msg, message: errorMessage } : msg
          ),
          streamingMessageId: null,
          isLoading: false,
        };
      }),

    reset: () =>
      set({
        messages: [
          {
            id: 1,
            avatar: '/assets/images/spectra-avatar.png',
            name: 'Spectrum',
            message:
              "Hello! I'm Spectrum, your spectral-analysis assistant – ask me anything about SMILES or these spectra.",
            timestamp: new Date().toISOString(),
            type: 'assistant',
          },
        ],
        isLoading: false,
        streamingMessageId: null,
        currentSmiles: null,
      }),
  };
}

export const chatStore = createChatStore();

// Derived stores
export const canSendMessage = derived(chatStore, ($store) => !$store.isLoading);

export const messageHistory = derived(chatStore, ($store): ApiChatMessage[] =>
  $store.messages.map(({ type, message }) => ({
    role: type,
    content: message,
  }))
);

export const isStreaming = derived(chatStore, ($store) => $store.streamingMessageId !== null);
