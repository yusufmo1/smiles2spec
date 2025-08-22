<script lang="ts">
  import type { ComponentType } from 'svelte';
  import { IconSizes, IconColors, type IconSize, type IconColor } from './iconConfig';

  export let icon: ComponentType;
  export let size: IconSize | number = 'md';
  export let color: IconColor | string = 'primary';
  export let className: string = '';

  $: iconSize = typeof size === 'number' ? size : IconSizes[size];
  $: iconColor = color in IconColors ? IconColors[color as IconColor] : color;
</script>

<svelte:component this={icon} size={iconSize} color={iconColor} {className} />
