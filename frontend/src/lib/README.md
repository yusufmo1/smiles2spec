# Frontend Library (`/lib`)

The core library for the SMILES2SPEC frontend application, containing all reusable components, services, and utilities.

## Structure

```
lib/
├── components/      # Reusable UI components
├── services/        # API and business logic
├── stores/          # State management
├── styles/          # Global styles and themes
├── utils/           # Helper functions
└── config.js        # Frontend configuration
```

## Components (`/components`)

### Organization

```
components/
├── panels/          # Panel-based UI system
├── icons/           # SVG icon components
├── smiles-input/    # SMILES input system
├── chat/            # Chat interface
└── [core components] # Header, Navbar, etc.
```

### Core Components

#### SmilesInput

Main input interface for molecular structures:

- Multi-line input support
- Drag-and-drop file upload
- AI-powered generation
- Bulk processing detection

#### SpectrumPlot

Interactive mass spectrum visualization:

- Plotly.js integration
- Responsive resizing
- Peak filtering
- Hover tooltips

#### Chat

AI chat interface for spectrum analysis:

- Streaming responses
- Context awareness
- Message history
- Markdown rendering

### Panel System

Modular UI panels for different features:

```svelte
<script>
  import { PanelRenderer } from '$lib/components/panels';

  const panels = [
    { component: SpectrumPlot, props: { data } },
    { component: StructurePanel, props: { smiles } },
  ];
</script>

<PanelRenderer {panels} />
```

### Icons

Comprehensive icon library:

- 45+ custom SVG icons
- Consistent styling
- TypeScript support
- Dynamic icon loading system
- Size/color props

## Services (`/services`)

### api.ts

Complete API client with TypeScript types:

```typescript
// Predict spectrum
const result = await predictSpectrum('CCO');

// Chat with streaming
await chatWithSpectrum(messages, {
  smiles: 'CCO',
  onChunk: (chunk) => console.log(chunk),
});

// Generate SMILES
const { smiles } = await generateRandomSmiles(3, 'aromatic compounds');
```

### plotlyService.ts

Plotly.js wrapper for consistent plots:

- Theme management
- Responsive handling
- Memory cleanup
- Type safety

### plotManager.ts

Global plot lifecycle management:

- Track active plots
- Handle resizing
- Cleanup on unmount
- Performance optimisation

## Stores (`/stores`) - Modular Architecture

The store system is designed with modularity and performance in mind, featuring lazy-loaded components and clear separation of concerns.

### index.ts

Main store entry point with smart exports:

```typescript
import {
  appState,
  focusedPanel,
  currentPage,
  initializePlots, // Lazy-loaded plot functionality
  triggerPlotResize,
} from '$lib/stores';

// Lazy load plot management when needed
onMount(async () => {
  const plotCleanup = await initializePlots();
  return () => plotCleanup?.();
});
```

### appState.ts

Global application state:

```typescript
import { appState } from '$lib/stores/appState';

// Set prediction data
appState.setPredictionData(result);

// Manage bulk processing
appState.setBulkList(smilesList);
appState.setBulkIndex(0);

// Console logging
appState.addConsoleEntry({
  type: 'success',
  message: 'Prediction complete',
});
```

### carouselStore.ts

Panel navigation and carousel functionality:

```typescript
import { focusedPanel, setCarouselMode, nextCarouselPanel } from '$lib/stores/carouselStore';

// Navigate panels
setCarouselMode('home', true);
nextCarouselPanel('home');
```

### panelStore.ts

Panel definitions and management:

- Panel registration and props
- Dynamic panel loading
- Category-based organization

### pageStore.ts

Page routing and navigation state:

- Current page tracking
- Navigation history
- Loading states

### plotEffects.ts (Lazy-Loaded)

Plot management with performance optimisation:

- Lazy loading to reduce initial bundle size
- Automatic plot resizing on carousel changes
- Memory cleanup and resource management

### Architecture Benefits

**Modular Design**: Each store has a single responsibility
**Performance**: Plot management is lazy-loaded only when needed
**Tree-Shaking**: Unused store code is excluded from bundles
**Maintainability**: Clear separation makes debugging easier
**Type Safety**: Full TypeScript support across all modules

## Styles (`/styles`)

### theme.css

CSS custom properties for theming:

```css
:root {
  --accent-primary: #7879ff;
  --surface-primary: #ffffff;
  --text-primary: #1c1c1e;
  /* ... more tokens */
}
```

### tokens.css

Design system tokens:

- Colors
- Spacing
- Typography
- Shadows
- Animations

## Utils (`/utils`)

### debounce.ts

Debounce function for performance:

```typescript
const debouncedSearch = debounce(search, 300);
```

### environment.ts

Environment detection and configuration

### iconMapping.ts

Dynamic icon resolution system

## Configuration

### config.js

Frontend configuration:

```javascript
export const config = {
  apiUrl: import.meta.env.VITE_API_URL,
  maxFileSize: 10 * 1024 * 1024, // 10MB
  supportedFormats: ['.csv', '.txt'],
  // ... more config
};
```

## Development Guidelines

### Component Creation

1. Create component file in appropriate directory
2. Add TypeScript props interface
3. Include JSDoc documentation
4. Export from index file

Example:

```svelte
<script lang="ts">
  /**
   * MyComponent - Brief description
   * @component
   */
  export let value: string;
  export let onChange: (value: string) => void;
</script>
```

### Service Integration

1. Define TypeScript interfaces
2. Implement error handling
3. Add loading states
4. Include retry logic

### State Management

1. Use stores for global state
2. Component state for local UI
3. Derived stores for computed values
4. Clear action methods

## Best Practices

### Performance

- Lazy load heavy components
- Debounce user inputs
- Cleanup resources on destroy
- Use `{#key}` for force updates

### Accessibility

- Semantic HTML elements
- ARIA labels where needed
- Keyboard navigation
- Focus management

### Type Safety

- Define all prop types
- Use TypeScript strictly
- Avoid `any` types
- Export interfaces

### Error Handling

- User-friendly messages
- Graceful degradation
- Loading states
- Retry mechanisms

## Testing

### Component Testing

```typescript
import { render } from '@testing-library/svelte';
import MyComponent from './MyComponent.svelte';

test('renders correctly', () => {
  const { getByText } = render(MyComponent, {
    props: { value: 'test' },
  });
  expect(getByText('test')).toBeInTheDocument();
});
```

### Store Testing

```typescript
import { get } from 'svelte/store';
import { myStore } from './myStore';

test('updates value', () => {
  myStore.setValue('new');
  expect(get(myStore).value).toBe('new');
});
```

## Building for Production

The library is optimised during build:

- Tree shaking removes unused code
- CSS is purged and minified
- Components are compiled to vanilla JS
- Assets are fingerprinted

## Future Enhancements

- Component library documentation
- Storybook integration
- Visual regression testing
- Performance monitoring
- Accessibility audit tools
