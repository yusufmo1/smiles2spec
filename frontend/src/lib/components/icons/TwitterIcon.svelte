<script lang="ts">
  export let size: number | string = 20;
  export let color: string = 'currentColor';
  export let className: string = '';
</script>

<svg
  width={size}
  height={size}
  viewBox="0 0 50 50"
  fill="none"
  xmlns="http://www.w3.org/2000/svg"
  class={className}
>
  <path
    d="M39,45H11c-3.314,0-6-2.686-6-6V11c0-3.314,2.686-6,6-6h28c3.314,0,6,2.686,6,6v28C45,42.314,42.314,45,39,45z"
    fill="none"
    stroke={color}
    stroke-width="2"
    stroke-miterlimit="10"
  />
  <path
    d="M37.914,37h-7.937L13.086,13h7.937L37.914,37z M31.022,35h3.063L19.978,15h-3.063L31.022,35z"
    fill={color}
  />
  <polygon points="15.5,37 23.973,27.08 22.844,25.476 13,37" fill={color} />
  <polygon points="26.115,21.646 27.244,23.251 36,13 33.5,13" fill={color} />
</svg>
