<script>
  import { panelDefinitions } from '../stores/panelData.js';
  import { carouselIndex } from '../stores.js';
  import Panel from './Panel.svelte';
  import { register } from 'swiper/element/bundle';
  import { onMount, onDestroy, tick } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  register();

  let swiperEl;
  let initialized = false;

  function handleSlideChange() {
    if (swiperEl?.swiper) {
      const newIndex = swiperEl.swiper.activeIndex;
      carouselIndex.set(newIndex);
      setTimeout(resizePlotsInActiveSlide, 300);
    }
  }

  function resizePlotsInActiveSlide() {
    const activeSlide = swiperEl?.querySelector('.swiper-slide-active');
    if (activeSlide) {
      const plots = activeSlide.querySelectorAll('.js-plotly-plot');
      plots.forEach(plot => Plotly.Plots.resize(plot));
    }
  }

  onMount(async () => {
    if (!swiperEl) return;

    Object.assign(swiperEl, {
      slidesPerView: 'auto',
      centeredSlides: true,
      effect: 'coverflow',
      grabCursor: true,
      initialSlide: $carouselIndex,
      watchSlidesProgress: true,
      loop: false,
      coverflowEffect: {
        rotate: 0,
        stretch: 120,
        depth: 250,
        modifier: 1,
        slideShadows: false
      }
    });

    await tick();
    swiperEl.initialize();

    if (swiperEl.swiper) {
      swiperEl.swiper.slideTo($carouselIndex, 0);
      initialized = true;
      setTimeout(resizePlotsInActiveSlide, 100);
    }

    swiperEl.addEventListener('slidechange', handleSlideChange);
    window.addEventListener('resize', resizePlotsInActiveSlide);
  });

  $: if (initialized && swiperEl?.swiper && $carouselIndex !== swiperEl.swiper.activeIndex) {
    swiperEl.swiper.slideTo($carouselIndex);
  }

  onDestroy(() => {
    window.removeEventListener('resize', resizePlotsInActiveSlide);
    swiperEl?.removeEventListener('slidechange', handleSlideChange);
    swiperEl?.swiper?.destroy(true, true);
  });
</script>

<div class="carousel-container">
  <swiper-container bind:this={swiperEl} init="false">
    {#each $panelDefinitions as panel (panel.id)}
      <swiper-slide>
        <div class="slide-content">
          <Panel {...panel} clickable={false} isCarousel={true} />
        </div>
      </swiper-slide>
    {/each}
  </swiper-container>
</div>

<style>
  .carousel-container {
    width: 100%;
    height: 100%;
    position: relative;
  }

  :global(swiper-container) {
    width: 100% !important;
    height: 100% !important;
  }

  :global(.swiper-slide) {
    height: auto !important;
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
  }

  .slide-content {
    width: 600px;
    max-width: 90vw;
    height: 500px;
    max-height: 80vh;
  }

  :global(.slide-content *) {
    pointer-events: auto;
  }

  @media (max-width: 768px) {
    .slide-content {
      width: 90vw;
      height: 70vh;
    }
  }
</style>
