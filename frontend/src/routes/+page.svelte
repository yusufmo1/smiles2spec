<script lang="ts">
  import HeroSection from '$lib/components/landing/HeroSection.svelte';
  import FeaturesGrid from '$lib/components/landing/FeaturesGrid.svelte';
  import TechnicalCapabilitiesSection from '$lib/components/landing/TechnicalCapabilitiesSection.svelte';
  import CallToActionSection from '$lib/components/landing/CallToActionSection.svelte';
  import QuickLinksSection from '$lib/components/landing/QuickLinksSection.svelte';
</script>

<svelte:head>
  <title>SMILES2SPEC - AI-Powered Mass Spectrum Prediction</title>
  <meta
    name="description"
    content="Transform SMILES molecular notation into accurate high-resolution mass spectra using advanced machine learning. Open-source tool for computational chemistry research."
  />
  <meta
    name="keywords"
    content="mass spectrometry, SMILES, machine learning, spectrum prediction, computational chemistry, open source"
  />
  <meta property="og:title" content="SMILES2SPEC - AI-Powered Mass Spectrum Prediction" />
  <meta
    property="og:description"
    content="Free, open-source tool for predicting mass spectra from SMILES notation using machine learning."
  />
  <meta property="og:type" content="website" />
</svelte:head>

<div class="landing-page">
  <HeroSection />
  <FeaturesGrid />
  <TechnicalCapabilitiesSection />
  <CallToActionSection />
  <QuickLinksSection />
</div>

<style>
  .landing-page {
    width: 100%;
    min-height: 100vh;
    display: flex;
    flex-direction: column;
    gap: 1rem;
  }

  /* Ensure smooth scrolling and proper spacing */
  :global(html) {
    scroll-behavior: smooth;
  }

  /* Add some breathing room */
  @media (max-width: 768px) {
    .landing-page {
      gap: 1rem;
    }
  }
</style>
