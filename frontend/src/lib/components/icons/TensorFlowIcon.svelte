<script lang="ts">
  import { IconSizes, type IconSize, IconColors, type IconColor } from './iconConfig';

  export let size: IconSize | number = 'md';
  export let color: IconColor = 'primary';

  $: sizeValue = typeof size === 'number' ? size : IconSizes[size];
  $: fill = IconColors[color] || color;
</script>

<svg width={sizeValue} height={sizeValue} viewBox="0 0 32 32" xmlns="http://www.w3.org/2000/svg">
  <path
    {fill}
    d="M29.399 14.959l-0.017-6.636-12.806-7.319v29.991l5.117-2.972v-8.438l3.864 2.234-0.023-5.771-3.841-2.194v-3.35l7.708 4.454zM2.619 8.322l12.806-7.318v29.991l-5.117-2.972v-17.519l-7.708 4.454 0.019-6.636z"
  />
</svg>
