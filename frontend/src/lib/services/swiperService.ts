import { browser } from '$app/environment';

/**
 * Swiper Singleton Service
 * 
 * Pre-initializes Swiper.js once and reuses it across all carousel components.
 * This eliminates the heavy bundle import and registration on every carousel activation.
 */
class SwiperService {
  private isRegistered = false;
  private registrationPromise: Promise<void> | null = null;

  /**
   * Initialize Swiper registration once globally
   * This handles the heavy lifting of importing and registering Swiper
   */
  async ensureRegistered(): Promise<void> {
    if (!browser) return;
    
    if (this.isRegistered) return;
    
    // Return existing promise if already loading
    if (this.registrationPromise) {
      return this.registrationPromise;
    }

    // Start registration process
    this.registrationPromise = this.doRegistration();
    return this.registrationPromise;
  }

  private async doRegistration(): Promise<void> {
    try {
      console.log('🎠 Initializing Swiper singleton...');
      const { register } = await import('swiper/element/bundle');
      register();
      this.isRegistered = true;
      console.log('✅ Swiper singleton ready');
    } catch (error) {
      console.error('❌ Failed to register Swiper:', error);
      this.registrationPromise = null; // Reset to allow retry
      throw error;
    }
  }

  /**
   * Configure a swiper element with optimised settings
   */
  configureSwiper(swiperEl: any, options: Record<string, any> = {}): void {
    if (!swiperEl) return;

    const defaultConfig = {
      slidesPerView: 1.25,
      centeredSlides: true,
      effect: 'coverflow',
      spaceBetween: 40,
      grabCursor: true,
      initialSlide: 0,
      // Optimised interaction settings
      allowTouchMove: true,
      simulateTouch: true,
      touchStartPreventDefault: false,
      preventClicks: false,
      preventClicksPropagation: false,
      slideToClickedSlide: false,
      touchReleaseOnEdges: true,
      resistance: true,
      resistanceRatio: 0.5,
      keyboard: {
        enabled: true,
        onlyInViewport: true,
        pageUpDown: false,
      },
      mousewheel: false,
      // Performance optimisations
      updateOnWindowResize: false, // We'll handle resize manually
      watchSlidesProgress: false,
      watchSlidesVisibility: false,
    };

    Object.assign(swiperEl, { ...defaultConfig, ...options });
  }

  /**
   * Check if Swiper is ready for use
   */
  isReady(): boolean {
    return this.isRegistered;
  }
}

// Export singleton instance
export const swiperService = new SwiperService();

// Pre-warm Swiper in browser environments for better UX
if (browser) {
  // Use requestIdleCallback if available, otherwise setTimeout
  const preWarm = () => swiperService.ensureRegistered();
  
  if ('requestIdleCallback' in window) {
    requestIdleCallback(preWarm);
  } else {
    setTimeout(preWarm, 100);
  }
}