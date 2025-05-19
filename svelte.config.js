import preprocess from 'svelte-preprocess';

export default {
  preprocess: preprocess({
    typescript: true,
    // optional: scss, postcss, etc.
  })
}; 