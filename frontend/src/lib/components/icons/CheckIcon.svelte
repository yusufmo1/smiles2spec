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
  <polyline points="20 6 9 17 4 12"></polyline>
</svg>
