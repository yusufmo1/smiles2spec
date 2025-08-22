# SVG Icon System (Bundle Optimised)

This directory contains a comprehensive SVG icon system for the Smiles2Spec application with optimised bundle loading. The system provides consistent icon usage with proper type safety, customisation options, and performance optimisation.

## Bundle Optimisation Features

- **Dynamic Loading**: Icons are loaded only when needed
- **Caching**: Loaded icons are cached for subsequent use
- **Bundle Splitting**: Prevents all 66+ icons from being bundled unnecessarily
- **Size Reduction**: Potential 50-100KB bundle size reduction

## Icon Categories

- **UI**: General UI controls (close, search, etc.)
- **Actions**: User actions (upload, download, generate, etc.)
- **Scientific**: Domain-specific icons (molecule, spectrum, etc.)
- **Navigation**: Direction indicators (arrows, chevrons, etc.)
- **Status**: Status indicators (check, error, loading, etc.)

## Usage (Bundle Optimized)

### Recommended: DynamicIcon Component

```svelte
<script>
  import { DynamicIcon } from '$lib/components/icons';
</script>

<!-- Basic usage with default size and color -->
<DynamicIcon name="CloseIcon" />

<!-- With custom size -->
<DynamicIcon name="UploadIcon" size={24} />

<!-- With custom color -->
<DynamicIcon name="GenerateIcon" size={20} color="#7879ff" />

<!-- With custom class -->
<DynamicIcon name="CloseIcon" className="my-custom-icon" />

<!-- In a button -->
<button>
  <DynamicIcon name="UploadIcon" size={16} />
  Upload File
</button>
```

### Alternative: Direct Import (for frequently used icons)

```svelte
<script>
  // Only import icons used across many components
  import CloseIcon from '$lib/components/icons/CloseIcon.svelte';
</script>

<CloseIcon size={24} />
```

## Props

All icon components accept the following props:

| Prop      | Type             | Default        | Description                    |
| --------- | ---------------- | -------------- | ------------------------------ |
| size      | number \| string | 20             | Size of icon (in px if number) |
| color     | string           | 'currentColor' | Color of icon                  |
| className | string           | ''             | Additional CSS classes         |

## Adding New Icons

To add a new icon, create a new Svelte component following this template:

```svelte
<!-- NewIcon.svelte -->
<script lang="ts">
  export let size: number | string = 20;
  export let color: string = 'currentColor';
  export let className: string = '';
</script>

<svg
  width={size}
  height={size}
  viewBox="0 0 24 24"
  fill="none"
  stroke={color}
  stroke-width="2"
  stroke-linecap="round"
  stroke-linejoin="round"
  class={className}
>
  <!-- SVG path data here -->
</svg>
```

Then add it to the `index.ts` file:

```typescript
export { default as NewIcon } from './NewIcon.svelte';
```

## Benefits of this Approach

1. **Type Safety**: Full TypeScript support for all icon components
2. **IDE Autocomplete**: Better developer experience with component discovery
3. **Performance**: SVG is embedded directly in the component (no extra network requests)
4. **Consistency**: All icons share the same props and behavior
5. **Customization**: Easy to override size, color, and add custom classes
