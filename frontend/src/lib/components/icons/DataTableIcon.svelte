<script lang="ts">
  export let size: number | string = 20;
  export let color: string = 'currentColor';
  export let className: string = '';
</script>

<svg
  width={size}
  height={size}
  viewBox="0 0 32 32"
  fill="none"
  stroke={color}
  stroke-width="0"
  class={className}
>
  <rect x="8" y="18" width="4" height="2" fill={color} />
  <rect x="14" y="18" width="4" height="2" fill={color} />
  <rect x="8" y="14" width="4" height="2" fill={color} />
  <rect x="14" y="22" width="4" height="2" fill={color} />
  <rect x="20" y="14" width="4" height="2" fill={color} />
  <rect x="20" y="22" width="4" height="2" fill={color} />
  <path
    d="M27,3H5A2.0025,2.0025,0,0,0,3,5V27a2.0025,2.0025,0,0,0,2,2H27a2.0025,2.0025,0,0,0,2-2V5A2.0025,2.0025,0,0,0,27,3Zm0,2,0,4H5V5ZM5,27V11H27l0,16Z"
    fill={color}
  />
</svg>
