# Landing Page Components

The landing page components provide a modular system for building an engaging and informative homepage for SMILES2SPEC, with clear navigation paths to the main application features.

## Architecture Overview

```
landing/
├── HeroSection.svelte          # Main hero/banner section
├── FeaturesGrid.svelte         # Features showcase grid
├── QuickLinksSection.svelte    # Navigation links to main tools
└── CallToActionSection.svelte  # Final CTA and getting started
```

## Component Overview

### HeroSection.svelte

The primary hero section that introduces SMILES2SPEC and its core value proposition:

```svelte
<script>
  import { HeroSection } from '$lib/components/landing';
</script>

<HeroSection
  title="SMILES2SPEC"
  subtitle="AI-Powered Mass Spectrum Prediction"
  description="Predict high-resolution electron ionization mass spectra directly from SMILES strings using machine learning."
  primaryAction="Try Prediction"
  secondaryAction="Learn More"
  on:primaryClick={navigateToSimulation}
  on:secondaryClick={navigateToHowItWorks}
/>
```

**Features:**

- Responsive hero layout
- Animated title and subtitle
- Clear value proposition
- Dual action buttons
- Background graphics/animations
- Mobile-optimised design

**Props:**

- `title` (string) - Main heading text
- `subtitle` (string) - Secondary heading
- `description` (string) - Value proposition text
- `primaryAction` (string) - Primary button text
- `secondaryAction` (string) - Secondary button text
- `showAnimation` (boolean) - Enable hero animations
- `backgroundVariant` (string) - Background style variant

**Events:**

- `primaryClick` - Primary button clicked
- `secondaryClick` - Secondary button clicked
- `heroVisible` - Hero section entered viewport

### FeaturesGrid.svelte

Showcases key features and capabilities in a grid layout:

```svelte
<FeaturesGrid
  title="Key Features"
  features={featuresData}
  columns={3}
  showIcons={true}
  on:featureClick={handleFeatureClick}
/>
```

**Features:**

- Responsive grid layout (1-4 columns)
- Icon integration with each feature
- Hover effects and animations
- Click-through to relevant sections
- Mobile-first responsive design
- Customizable grid spacing

**Feature Data Structure:**

```typescript
interface Feature {
  id: string;
  title: string;
  description: string;
  icon: string; // Icon component name
  link?: string; // Optional navigation link
  highlighted?: boolean; // Highlight important features
  comingSoon?: boolean; // Mark upcoming features
}

const featuresData: Feature[] = [
  {
    id: 'ml-prediction',
    title: 'ML-Powered Prediction',
    description:
      'Advanced machine learning models trained on 280,000+ spectra for accurate predictions.',
    icon: 'CpuIcon',
    link: '/spectral-simulation',
  },
  {
    id: 'real-time-viz',
    title: 'Real-time Visualization',
    description:
      'Interactive spectrum plots with peak identification and molecular structure display.',
    icon: 'ChartIcon',
    highlighted: true,
  },
  {
    id: 'ai-chat',
    title: 'AI-Powered Analysis',
    description: 'Chat with AI about your spectra for interpretation and insights.',
    icon: 'ChatWithSpectrumIcon',
    link: '/chat-with-spectrum',
  },
  // ... more features
];
```

**Grid Layouts:**

- **Desktop (>1024px)**: 3-4 columns
- **Tablet (768-1024px)**: 2 columns
- **Mobile (<768px)**: 1 column

### QuickLinksSection.svelte

Provides direct navigation to main application features:

```svelte
<QuickLinksSection
  title="Get Started"
  links={quickLinks}
  layout="cards"
  showDescriptions={true}
  on:linkClick={handleNavigation}
/>
```

**Features:**

- Multiple layout options (cards, list, grid)
- Icon-based navigation
- Brief descriptions for each link
- Visual hierarchy for primary actions
- Responsive design
- Loading states for navigation

**Link Data Structure:**

```typescript
interface QuickLink {
  id: string;
  title: string;
  description: string;
  icon: string;
  href: string;
  isPrimary?: boolean; // Primary action styling
  isExternal?: boolean; // External link indicator
  badge?: string; // Optional badge text
}

const quickLinks: QuickLink[] = [
  {
    id: 'predict',
    title: 'Predict Spectrum',
    description: 'Enter a SMILES string and get instant mass spectrum prediction',
    icon: 'SpectrumIcon',
    href: '/spectral-simulation',
    isPrimary: true,
  },
  {
    id: 'chat',
    title: 'Chat with AI',
    description: 'Ask questions about mass spectrometry and get expert insights',
    icon: 'ChatWithSpectrumIcon',
    href: '/chat-with-spectrum',
  },
  {
    id: 'learn',
    title: 'How It Works',
    description: 'Learn about the science and technology behind the predictions',
    icon: 'InfoIcon',
    href: '/how-it-works',
  },
  {
    id: 'about',
    title: 'About Project',
    description: 'Meet the developer and learn about the project',
    icon: 'UserIcon',
    href: '/about',
  },
];
```

**Layout Options:**

- `cards` - Card-based layout with icons and descriptions
- `grid` - Simple grid layout
- `list` - Vertical list layout
- `hero-cards` - Large card layout for primary actions

### CallToActionSection.svelte

Final section encouraging user engagement and conversion:

```svelte
<CallToActionSection
  title="Ready to Explore?"
  description="Start predicting mass spectra from SMILES strings today"
  primaryAction="Start Predicting"
  secondaryAction="View Examples"
  showBenefits={true}
  benefits={ctaBenefits}
  on:primaryClick={navigateToSimulation}
  on:secondaryClick={showExamples}
/>
```

**Features:**

- Compelling call-to-action messaging
- Benefits list highlighting
- Multiple action buttons
- Social proof elements (optional)
- Newsletter signup (optional)
- Contact information

**Benefits Structure:**

```typescript
interface Benefit {
  id: string;
  text: string;
  icon?: string;
  highlighted?: boolean;
}

const ctaBenefits: Benefit[] = [
  {
    id: 'instant',
    text: 'Instant predictions from any SMILES string',
    icon: 'FlashIcon',
  },
  {
    id: 'accurate',
    text: 'High accuracy with R² > 0.85 performance',
    icon: 'TargetIcon',
    highlighted: true,
  },
  {
    id: 'free',
    text: 'Completely free and open-source',
    icon: 'HeartIcon',
  },
];
```

## Responsive Design

### Breakpoint Strategy

```css
/* Mobile First Approach */
.section {
  padding: 1rem;
}

/* Tablet */
@media (min-width: 768px) {
  .section {
    padding: 2rem;
  }
}

/* Desktop */
@media (min-width: 1024px) {
  .section {
    padding: 3rem;
  }
}

/* Large Desktop */
@media (min-width: 1440px) {
  .section {
    padding: 4rem;
  }
}
```

### Component Adaptation

Each component adapts its layout based on screen size:

#### HeroSection Responsive Behavior

- **Mobile**: Single column, reduced padding, smaller typography
- **Tablet**: Maintain single column with larger spacing
- **Desktop**: Potential two-column layout with illustration
- **Large**: Full-width hero with maximum visual impact

#### FeaturesGrid Responsive Behavior

- **Mobile**: 1 column, stack features vertically
- **Tablet**: 2 columns, maintain readability
- **Desktop**: 3-4 columns, optimal information density
- **Large**: 4 columns with maximum spacing

## Animation and Interactions

### Entrance Animations

Components support viewport-based animations:

```typescript
// Intersection Observer for scroll animations
import { onMount } from 'svelte';

let sectionRef: HTMLElement;
let isVisible = false;

onMount(() => {
  const observer = new IntersectionObserver(
    (entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          isVisible = true;
          observer.unobserve(entry.target);
        }
      });
    },
    {
      threshold: 0.1,
      rootMargin: '50px',
    }
  );

  if (sectionRef) {
    observer.observe(sectionRef);
  }

  return () => observer.disconnect();
});
```

### Hover Effects

Interactive elements include subtle hover effects:

```css
.feature-card {
  transition:
    transform 0.2s ease,
    box-shadow 0.2s ease;
}

.feature-card:hover {
  transform: translateY(-4px);
  box-shadow: 0 8px 25px rgba(120, 121, 255, 0.15);
}

.quick-link:hover {
  background: var(--surface-secondary);
  border-color: var(--accent-primary);
}
```

### Loading States

Components handle loading states gracefully:

```svelte
<script>
  export let isLoading = false;
</script>

<section class="features-grid" class:loading={isLoading}>
  {#if isLoading}
    <!-- Skeleton loading state -->
    {#each Array(6) as _, i}
      <div class="feature-skeleton">
        <div class="skeleton-icon"></div>
        <div class="skeleton-title"></div>
        <div class="skeleton-description"></div>
      </div>
    {/each}
  {:else}
    <!-- Actual features -->
    {#each features as feature}
      <FeatureCard {feature} on:click={handleFeatureClick} />
    {/each}
  {/if}
</section>
```

## Content Management

### Dynamic Content Support

Components support dynamic content loading:

```typescript
// Content configuration
interface LandingConfig {
  hero: {
    title: string;
    subtitle: string;
    description: string;
    actions: ActionButton[];
  };
  features: Feature[];
  quickLinks: QuickLink[];
  cta: {
    title: string;
    description: string;
    benefits: Benefit[];
  };
}

// Load from CMS or configuration
async function loadLandingContent(): Promise<LandingConfig> {
  // Could load from headless CMS, API, or local config
  return await fetch('/api/landing-content').then((r) => r.json());
}
```

### Internationalization Support

Components are designed for i18n compatibility:

```svelte
<script>
  import { _ } from 'svelte-i18n';
</script>

<h1>{$_('landing.hero.title')}</h1><p>{$_('landing.hero.description')}</p>
```

## SEO and Performance

### Meta Data Integration

Landing page components support SEO optimisation:

```svelte
<!-- In +page.svelte -->
<svelte:head>
  <title>{pageTitle} | SMILES2SPEC</title>
  <meta name="description" content={pageDescription} />
  <meta property="og:title" content={pageTitle} />
  <meta property="og:description" content={pageDescription} />
  <meta property="og:image" content={ogImage} />
  <meta name="twitter:card" content="summary_large_image" />
</svelte:head>
```

### Performance Optimization

```typescript
// Lazy loading for non-critical components
import { onMount } from 'svelte';

let FeaturesGrid: any;
let CallToActionSection: any;

onMount(async () => {
  // Load below-the-fold components after initial render
  const [featuresModule, ctaModule] = await Promise.all([
    import('./FeaturesGrid.svelte'),
    import('./CallToActionSection.svelte'),
  ]);

  FeaturesGrid = featuresModule.default;
  CallToActionSection = ctaModule.default;
});
```

### Image Optimization

```svelte
<script>
  import { onMount } from 'svelte';

  let imageLoaded = false;

  function handleImageLoad() {
    imageLoaded = true;
  }
</script>

<div class="hero-image">
  <img
    src="/images/hero-spectrum.webp"
    alt="Mass spectrum visualization"
    loading="lazy"
    on:load={handleImageLoad}
    class:loaded={imageLoaded}
  />
  {#if !imageLoaded}
    <div class="image-skeleton"></div>
  {/if}
</div>
```

## Accessibility

### Semantic HTML

Components use proper semantic markup:

```svelte
<section aria-labelledby="features-heading">
  <h2 id="features-heading">Key Features</h2>

  <div role="grid" aria-label="Feature grid">
    {#each features as feature, index}
      <article role="gridcell" tabindex="0">
        <h3>{feature.title}</h3>
        <p>{feature.description}</p>
      </article>
    {/each}
  </div>
</section>
```

### Keyboard Navigation

Full keyboard accessibility:

```typescript
function handleKeyDown(event: KeyboardEvent, action: () => void) {
  if (event.key === 'Enter' || event.key === ' ') {
    event.preventDefault();
    action();
  }
}
```

### Screen Reader Support

```svelte
<button aria-label="Navigate to spectrum prediction tool" aria-describedby="predict-description">
  Start Predicting
</button>
<div id="predict-description" class="sr-only">
  Opens the spectral simulation page where you can enter SMILES strings
</div>
```

## Usage Examples

### Complete Landing Page

```svelte
<script>
  import {
    HeroSection,
    FeaturesGrid,
    QuickLinksSection,
    CallToActionSection,
  } from '$lib/components/landing';

  import { goto } from '$app/navigation';

  const features = [
    // feature data
  ];

  const quickLinks = [
    // link data
  ];

  function navigateToSimulation() {
    goto('/spectral-simulation');
  }

  function navigateToHowItWorks() {
    goto('/how-it-works');
  }
</script>

<main class="landing-page">
  <HeroSection
    title="SMILES2SPEC"
    subtitle="AI-Powered Mass Spectrum Prediction"
    description="Predict high-resolution electron ionization mass spectra directly from SMILES strings using machine learning."
    primaryAction="Try Prediction"
    secondaryAction="Learn More"
    on:primaryClick={navigateToSimulation}
    on:secondaryClick={navigateToHowItWorks}
  />

  <FeaturesGrid title="Why Choose SMILES2SPEC?" {features} columns={3} showIcons={true} />

  <QuickLinksSection
    title="Get Started"
    links={quickLinks}
    layout="cards"
    showDescriptions={true}
    on:linkClick={({ detail }) => goto(detail.href)}
  />

  <CallToActionSection
    title="Ready to Predict Spectra?"
    description="Join researchers worldwide using SMILES2SPEC for mass spectrum prediction"
    primaryAction="Start Now"
    secondaryAction="View Documentation"
    on:primaryClick={navigateToSimulation}
    on:secondaryClick={() => goto('/how-it-works')}
  />
</main>
```

### Custom Configuration

```svelte
<script>
  import { HeroSection } from '$lib/components/landing';

  const customHeroConfig = {
    title: 'Welcome to Advanced Spectroscopy',
    subtitle: 'Cutting-edge molecular analysis',
    description:
      'Our platform provides state-of-the-art tools for mass spectrometry prediction and analysis.',
    primaryAction: 'Explore Tools',
    secondaryAction: 'Learn More',
    showAnimation: true,
    backgroundVariant: 'gradient',
  };
</script>

<HeroSection
  {...customHeroConfig}
  on:primaryClick={handlePrimaryAction}
  on:secondaryClick={handleSecondaryAction}
/>
```

## Testing

### Component Testing

```javascript
import { render, fireEvent } from '@testing-library/svelte';
import HeroSection from './HeroSection.svelte';

describe('HeroSection', () => {
  test('renders with required props', () => {
    const { getByText } = render(HeroSection, {
      props: {
        title: 'Test Title',
        subtitle: 'Test Subtitle',
        description: 'Test Description',
      },
    });

    expect(getByText('Test Title')).toBeInTheDocument();
    expect(getByText('Test Subtitle')).toBeInTheDocument();
    expect(getByText('Test Description')).toBeInTheDocument();
  });

  test('triggers navigation events', async () => {
    const primaryClick = vi.fn();
    const { getByTestId, component } = render(HeroSection, {
      props: {
        title: 'Test',
        primaryAction: 'Click Me',
      },
    });

    component.$on('primaryClick', primaryClick);

    await fireEvent.click(getByTestId('primary-action'));
    expect(primaryClick).toHaveBeenCalled();
  });
});
```

### Accessibility Testing

```javascript
import { axe, toHaveNoViolations } from 'jest-axe';

expect.extend(toHaveNoViolations);

test('landing components are accessible', async () => {
  const { container } = render(HeroSection, {
    props: { title: 'Test', subtitle: 'Test' },
  });

  const results = await axe(container);
  expect(results).toHaveNoViolations();
});
```

## Best Practices

### Content Strategy

1. **Clear value proposition** in hero section
2. **Benefit-focused** feature descriptions
3. **Action-oriented** call-to-action text
4. **Scannable content** with proper hierarchy

### Performance

1. **Lazy load** below-the-fold components
2. **Optimize images** with WebP format
3. **Minimize layout shift** with skeleton states
4. **Progressive enhancement** for animations

### User Experience

1. **Clear navigation paths** to main features
2. **Consistent visual hierarchy** across sections
3. **Mobile-first design** approach
4. **Loading states** for all interactions

### Accessibility

1. **Semantic HTML** structure
2. **Keyboard navigation** support
3. **Screen reader** compatibility
4. **High contrast** color schemes

## Future Enhancements

- **A/B Testing**: Component variants for conversion optimisation
- **Analytics Integration**: Track user interactions and conversions
- **Progressive Web App**: Enhanced mobile experience
- **Content Management**: Dynamic content updates
- **Personalization**: User-specific content recommendations
- **Multi-language**: International market support
