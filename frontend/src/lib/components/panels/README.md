# Panel System

The panel system is the core UI architecture for SMILES2SPEC, providing a modular, reusable component system for organizing different types of content across pages.

## Architecture Overview

```
panels/
├── Panel.svelte           # Base panel component
├── PanelRenderer.svelte   # Panel orchestration system
├── index.js               # Main exports
├── about/                 # About page panels
├── chat/                  # Chat-specific panels
├── components/            # Shared panel components
├── info/                  # Information and documentation panels
├── layout/                # Panel layout systems
├── simulation/            # Spectrum simulation panels
└── utils/                 # Panel utilities and navigation
```

## Core Components

### Panel.svelte

Base panel component that provides consistent styling and behavior:

```svelte
<script>
  import Panel from '$lib/components/panels/Panel.svelte';
</script>

<Panel title="My Panel" className="custom-panel">
  <p>Panel content goes here</p>
</Panel>
```

**Props:**

- `title` (string) - Panel header title
- `className` (string) - Additional CSS classes
- `padding` (boolean) - Whether to include default padding

### PanelRenderer.svelte

Dynamic panel rendering system that manages panel lifecycle:

```svelte
<script>
  import { PanelRenderer } from '$lib/components/panels';

  const panels = [
    { component: 'SpectrumPlot', props: { data: spectrumData } },
    { component: 'StructurePanel', props: { smiles: 'CCO' } },
  ];
</script>

<PanelRenderer {panels} category="simulation" />
```

**Features:**

- Dynamic component loading
- Prop management
- Error boundaries
- Performance optimization

## Panel Categories

### Simulation Panels (`/simulation`)

Core functionality for spectrum prediction and visualization:

- **SpectrumPlot** - Interactive mass spectrum visualization
- **StructurePanel** - Molecular structure display
- **PeakTable** - Tabular peak data
- **ConsoleOutput** - Processing logs and status
- **ExportCenter** - Data export functionality
- **IonFragmentation** - Fragmentation pattern analysis

### Information Panels (`/info`)

Educational and technical content:

- **SystemOverviewPanel** - High-level system explanation
- **WorkflowPanel** - Step-by-step process guide
- **MLEnginePanel** - Machine learning details
- **AccuracyValidationPanel** - Model performance metrics
- **TechnicalRequirementsPanel** - System requirements
- **AdvancedFeaturesPanel** - Advanced functionality overview
- **FragmentationSciencePanel** - Scientific background
- **VisualizationFeaturesPanel** - Visualization capabilities

### About Panels (`/about`)

Project and developer information:

- **DeveloperHeroPanel** - Developer introduction
- **DeveloperJourneyPanel** - Development story
- **ProjectMetricsPanel** - Project statistics
- **TechStackBanner** - Technology showcase
- **FeaturesShowcasePanel** - Feature highlights
- **ContactConnectRevampedPanel** - Contact information
- **AcknowledgmentsStreamlinedPanel** - Credits and thanks
- **OpenSourceCommitmentPanel** - Open source philosophy

### Chat Panels (`/chat`)

Chat-specific functionality:

- **CompoundStatsPanel** - Molecular statistics display

### Layout Panels (`/layout`)

Panel organization systems:

- **PanelCarousel** - Horizontal scrolling panel navigation
- **PanelGrid** - Grid-based panel layout
- **PanelOverlay** - Modal and overlay panels

## Shared Components (`/components`)

Reusable panel building blocks:

- **InnerContainer** - Standard content wrapper
- **MiniPanel** - Compact panel variant
- **PlaceholderState** - Loading and empty states

### Info Panel Components (`/info/components`)

Specialized components for information panels:

- **BentoGrid** - Grid layout system
- **MoleculeViewer** - 3D molecular visualization
- **ProcessFlow** - Workflow diagram component

## Panel System Features

### Dynamic Loading

Panels can be loaded dynamically based on configuration:

```javascript
// In panelStore.ts
const panels = {
  simulation: [
    { id: 'spectrum', component: 'SpectrumPlot', title: 'Mass Spectrum' },
    { id: 'structure', component: 'StructurePanel', title: 'Molecular Structure' },
  ],
};
```

### Responsive Design

Panels automatically adapt to different screen sizes:

- Desktop: Grid layout with multiple panels visible
- Tablet: Carousel navigation
- Mobile: Single panel with navigation controls

### Navigation System

Panel navigation is handled through the carousel store:

```javascript
import { setCarouselMode, nextCarouselPanel } from '$lib/stores/carouselStore';

// Enter carousel mode
setCarouselMode('simulation', true);

// Navigate panels
nextCarouselPanel('simulation');
```

### State Management

Panel state is managed through:

- **Panel Store**: Panel definitions and props
- **Carousel Store**: Navigation and focus
- **App State**: Data and content

## Usage Patterns

### Static Panel Usage

For fixed panels with static content:

```svelte
<script>
  import { SystemOverviewPanel } from '$lib/components/panels/info';
</script>

<SystemOverviewPanel />
```

### Dynamic Panel Rendering

For dynamic panel systems:

```svelte
<script>
  import { PanelRenderer } from '$lib/components/panels';
  import { simulationPanels } from '$lib/stores/panelStore';
</script>

<PanelRenderer panels={$simulationPanels} category="simulation" />
```

### Custom Panel Creation

To create a new panel:

1. Create component in appropriate category directory
2. Add to category index.js file
3. Register in panelStore.ts
4. Add navigation logic if needed

Example:

```svelte
<!-- MyCustomPanel.svelte -->
<script>
  import Panel from '../Panel.svelte';

  export let title = 'Custom Panel';
  export let data = null;
</script>

<Panel {title}>
  <div class="custom-content">
    <!-- Panel content -->
  </div>
</Panel>

<style>
  .custom-content {
    /* Panel-specific styles */
  }
</style>
```

## Performance Considerations

### Lazy Loading

Panels can be lazy-loaded to improve initial page load:

```javascript
// Dynamic import in PanelRenderer
const component = await import(`./simulation/${panelConfig.component}.svelte`);
```

### Virtualization

For large panel lists, implement virtual scrolling:

```svelte
<!-- Only render visible panels -->
{#each visiblePanels as panel}
  <svelte:component this={panel.component} {...panel.props} />
{/each}
```

### Memory Management

Panels clean up resources when destroyed:

```svelte
<script>
  import { onDestroy } from 'svelte';

  onDestroy(() => {
    // Cleanup plot instances, event listeners, etc.
  });
</script>
```

## Styling System

### CSS Variables

Panels use consistent CSS tokens:

```css
.panel {
  background: var(--surface-primary);
  border: 1px solid var(--surface-stroke);
  border-radius: var(--radius-md);
  padding: var(--spacing-lg);
}
```

### Theme Support

Panels automatically adapt to theme changes through CSS variables.

### Animation System

Panel transitions are handled in `/utils/animations.css`:

- Slide transitions for carousel navigation
- Fade effects for panel switching
- Loading state animations

## Best Practices

### Panel Design

1. **Single Responsibility**: Each panel should focus on one specific function
2. **Consistent Interface**: Use standard props and events
3. **Responsive Design**: Ensure panels work on all screen sizes
4. **Accessibility**: Include proper ARIA labels and keyboard navigation

### Performance

1. **Lazy Load**: Only load panels when needed
2. **Cleanup**: Properly dispose of resources
3. **Memoization**: Cache expensive computations
4. **Virtual Scrolling**: For large lists of panels

### State Management

1. **Minimize Props**: Use stores for complex state
2. **Reactive Updates**: Use reactive statements for derived data
3. **Event Communication**: Use events for panel-to-panel communication

## Testing

### Component Testing

Test individual panels in isolation:

```javascript
import { render } from '@testing-library/svelte';
import SpectrumPlot from './SpectrumPlot.svelte';

test('renders spectrum plot', () => {
  const { getByTestId } = render(SpectrumPlot, {
    props: { data: mockSpectrumData },
  });

  expect(getByTestId('spectrum-plot')).toBeInTheDocument();
});
```

### Integration Testing

Test panel system interactions:

```javascript
import { render, fireEvent } from '@testing-library/svelte';
import PanelRenderer from './PanelRenderer.svelte';

test('panel navigation works', async () => {
  const { getByTestId } = render(PanelRenderer, {
    props: { panels: mockPanels },
  });

  await fireEvent.click(getByTestId('next-panel'));
  // Assert panel changed
});
```

## Future Enhancements

- **Drag & Drop**: Panel reordering
- **Split Panes**: Resizable panel layouts
- **Panel Templates**: Pre-configured panel combinations
- **Progressive Enhancement**: Better offline support
- **Panel Persistence**: Save panel configurations
