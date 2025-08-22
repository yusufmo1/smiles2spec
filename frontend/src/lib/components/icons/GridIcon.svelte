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
  <!-- 3x3 grid representation -->
  <rect x="3" y="3" width="5" height="5" rx="1" />
  <rect x="10" y="3" width="5" height="5" rx="1" />
  <rect x="17" y="3" width="5" height="5" rx="1" />
  <rect x="3" y="10" width="5" height="5" rx="1" />
  <rect x="10" y="10" width="5" height="5" rx="1" />
  <rect x="17" y="10" width="5" height="5" rx="1" />
  <rect x="3" y="17" width="5" height="5" rx="1" />
  <rect x="10" y="17" width="5" height="5" rx="1" />
  <rect x="17" y="17" width="5" height="5" rx="1" />
</svg>
