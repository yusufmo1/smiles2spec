import { sveltekit } from '@sveltejs/kit/vite';

/** @type {import('vite').UserConfig} */
const config = {
  plugins: [sveltekit()],
  optimizeDeps: {
    include: ['plotly.js-dist-min'],
  },
  build: {
    rollupOptions: {
      output: {
        manualChunks: {
          plotly: ['plotly.js-dist-min'],
          panels: ['$lib/components/panels/index.js'],
        },
      },
    },
    chunkSizeWarningLimit: 1000,
  },
  resolve: {
    alias: {
      $icons: '/src/lib/components/icons',
    },
  },
  ssr: {
    noExternal: ['plotly.js-dist-min'],
  },
  server: {
    proxy: {
      '/api': {
        target: 'http://localhost:5050',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, ''), // Remove /api prefix when forwarding to backend
        // Note: This proxy only works in development
        // Production uses direct subdomain calls
      },
    },
  },
};

export default config;
