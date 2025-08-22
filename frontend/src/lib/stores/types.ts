// Core page and navigation types
export type PageKey =
  | 'home'
  | 'spectral-simulation'
  | 'how-it-works'
  | 'about'
  | 'chat-with-spectrum';

// Panel and carousel types
export interface PanelData {
  id: string;
  title: string;
  content?: any;
  type?: 'info' | 'simulation' | 'about';
}

export interface CarouselState {
  indices: Record<PageKey, number>;
  modes: Record<PageKey, boolean>;
}

// Application state types
export interface AppState {
  currentPage: PageKey;
  isLoading: boolean;
  error: string | null;
}

// Navigation types
export interface NavigationState {
  canGoBack: boolean;
  history: string[];
  currentRoute: string;
}

// Spectrum and prediction types (for future use)
export interface SpectrumData {
  x: number[];
  y: number[];
  name?: string;
  peaks?: Array<{ mz: number; intensity: number }>;
}

export interface PredictionState {
  isLoading: boolean;
  error: string | null;
  data: SpectrumData | null;
  hasFirstPrediction: boolean;
}

// Chat system types (for future use)
export interface ChatMessage {
  id: string;
  role: 'user' | 'assistant';
  content: string;
  timestamp: Date;
  smiles?: string;
}

export interface ChatState {
  messages: ChatMessage[];
  isLoading: boolean;
  canSendMessage: boolean;
  error: string | null;
}

// SMILES input types (for future use)
export interface SmilesInputState {
  value: string;
  isValid: boolean;
  error: string | null;
  history: string[];
}

// Store action types
export type StoreAction<T = any> = {
  type: string;
  payload?: T;
};
