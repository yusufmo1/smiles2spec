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
  <rect x="2" y="4" width="20" height="16" rx="2" />
  <path d="M22 7l-10 7L2 7" />
</svg>
