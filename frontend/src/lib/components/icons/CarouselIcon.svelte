<script lang="ts">
  export let size: number = 20;
  export let color: string = 'currentColor';
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
>
  <!-- Carousel/slides representation -->
  <rect x="2" y="6" width="4" height="12" rx="1" opacity="0.4" />
  <rect x="7" y="4" width="10" height="16" rx="2" />
  <rect x="18" y="6" width="4" height="12" rx="1" opacity="0.4" />
  <!-- Navigation dots -->
  <circle cx="10" cy="22" r="1" />
  <circle cx="12" cy="22" r="1" fill={color} />
  <circle cx="14" cy="22" r="1" />
</svg>
