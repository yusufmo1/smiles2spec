<script lang="ts">
  export let size: number | string = 24;
  export let color: string = 'currentColor';
  export let className: string = '';
</script>

<svg
  width={size}
  height={size}
  viewBox="0 0 24 24"
  fill="none"
  class={className}
  xmlns="http://www.w3.org/2000/svg"
>
  <path
    d="M3 3v18h18"
    stroke={color}
    stroke-width="2"
    stroke-linecap="round"
    stroke-linejoin="round"
  />
  <path
    d="M18 17l-5-5-4 4-5-5"
    stroke={color}
    stroke-width="2"
    stroke-linecap="round"
    stroke-linejoin="round"
  />
</svg>
