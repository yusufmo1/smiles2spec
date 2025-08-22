<script lang="ts">
  import { IconSizes, type IconSize, IconColors, type IconColor } from './iconConfig';

  export let size: IconSize | number = 'md';
  export let color: IconColor = 'primary';

  $: sizeValue = typeof size === 'number' ? size : IconSizes[size];
  $: fill = IconColors[color] || color;
</script>

<svg width={sizeValue} height={sizeValue} viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
  <!-- Central Carbon atom -->
  <circle cx="12" cy="12" r="2.5" {fill} />

  <!-- Four Hydrogen atoms positioned tetrahedrally -->
  <circle cx="12" cy="5" r="1.5" {fill} />
  <circle cx="19" cy="12" r="1.5" {fill} />
  <circle cx="12" cy="19" r="1.5" {fill} />
  <circle cx="5" cy="12" r="1.5" {fill} />

  <!-- Bonds from Carbon to each Hydrogen -->
  <line x1="12" y1="9.5" x2="12" y2="6.5" stroke={fill} stroke-width="2" />
  <line x1="14.5" y1="12" x2="17.5" y2="12" stroke={fill} stroke-width="2" />
  <line x1="12" y1="14.5" x2="12" y2="17.5" stroke={fill} stroke-width="2" />
  <line x1="9.5" y1="12" x2="6.5" y2="12" stroke={fill} stroke-width="2" />
</svg>
