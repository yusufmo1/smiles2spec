<script lang="ts">
  import type { ComponentType } from 'svelte';
  import { onMount } from 'svelte';
  import { ICON_MAP, type IconName } from './index';
  import { IconSizes, IconColors, type IconSize, type IconColor } from './iconConfig';

  export let name: IconName;
  export let size: IconSize | number = 'md';
  export let color: IconColor | string = 'primary';
  export let className: string = '';

  let IconComponent: ComponentType | null = null;
  let loading = true;
  let error = false;

  $: iconSize = typeof size === 'number' ? size : IconSizes[size];
  $: iconColor = color in IconColors ? IconColors[color as IconColor] : color;

  async function loadIconFromMap() {
    loading = true;
    error = false;
    IconComponent = null;

    try {
      // Use the ICON_MAP instead of dynamic string import
      const iconLoader = ICON_MAP[name];
      if (!iconLoader) {
        throw new Error(`Icon "${name}" not found in ICON_MAP`);
      }

      IconComponent = await iconLoader();
      loading = false;
    } catch (err) {
      console.error(`Failed to load icon: ${name}`, err);
      error = true;
      loading = false;
    }
  }

  onMount(() => {
    loadIconFromMap();
  });

  // Reactive loading when icon name changes
  $: if (name) {
    loadIconFromMap();
  }
</script>

{#if loading}
  <div
    class="icon-placeholder {className}"
    style="width: {iconSize}px; height: {iconSize}px; background-color: {iconColor}; opacity: 0.3;"
  ></div>
{:else if error}
  <div
    class="icon-error {className}"
    style="width: {iconSize}px; height: {iconSize}px; color: {iconColor};"
    title="Failed to load icon: {name}"
  >
    ?
  </div>
{:else if IconComponent}
  <svelte:component this={IconComponent} size={iconSize} color={iconColor} {className} />
{/if}

<style>
  .icon-placeholder {
    border-radius: 2px;
    animation: pulse 1.5s ease-in-out infinite;
  }

  .icon-error {
    display: flex;
    align-items: center;
    justify-content: center;
    border: 1px solid currentColor;
    border-radius: 2px;
    font-size: 10px;
    font-weight: bold;
  }

  @keyframes pulse {
    0%,
    100% {
      opacity: 0.3;
    }
    50% {
      opacity: 0.6;
    }
  }
</style>
