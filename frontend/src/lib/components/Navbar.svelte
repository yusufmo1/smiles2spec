<script lang="ts">
  import { page } from '$app/stores';
  import MenuIcon from './icons/MenuIcon.svelte';
  import CloseIcon from './icons/CloseIcon.svelte';
  import { onMount } from 'svelte';

  let isMobileMenuOpen = false;
  let isDesktop = false;

  // Determine active page based on current location using SvelteKit's page store
  $: activePage = getActivePage($page.url.pathname);

  function getActivePage(pathname: string) {
    if (pathname === '/') return 'home';
    if (pathname === '/spectral-simulation') return 'spectral-simulation';
    if (pathname === '/how-it-works') return 'how-it-works';
    if (pathname === '/about') return 'about';
    if (pathname === '/chat-with-spectrum') return 'chat-with-spectrum';
    return 'home';
  }

  function toggleMobileMenu() {
    isMobileMenuOpen = !isMobileMenuOpen;
    // Prevent body scroll when menu is open
    document.body.style.overflow = isMobileMenuOpen ? 'hidden' : '';
  }

  function closeMobileMenu() {
    isMobileMenuOpen = false;
    document.body.style.overflow = '';
  }

  function handleResize() {
    isDesktop = window.innerWidth >= 768;
    if (isDesktop && isMobileMenuOpen) {
      closeMobileMenu();
    }
  }

  onMount(() => {
    handleResize();
    window.addEventListener('resize', handleResize);
    return () => {
      window.removeEventListener('resize', handleResize);
      document.body.style.overflow = '';
    };
  });
</script>

<nav class="navbar glass-card">
  <div class="navbar-container">
    <div class="logo">
      <span class="text-gradient prosto-one-regular">SMILES2SPEC</span>
    </div>

    <!-- Desktop Navigation -->
    <ul class="nav-links desktop-nav">
      <li>
        <a
          href="/"
          class={activePage === 'home' ? 'active' : ''}
          data-sveltekit-preload-data="hover"
        >
          HOME
        </a>
      </li>
      <li>
        <a
          href="/spectral-simulation"
          class={activePage === 'spectral-simulation' ? 'active' : ''}
          data-sveltekit-preload-data="hover"
        >
          SPECTRAL SIMULATION
        </a>
      </li>
      <li>
        <a
          href="/chat-with-spectrum"
          class={activePage === 'chat-with-spectrum' ? 'active' : ''}
          data-sveltekit-preload-data="hover"
        >
          CHAT WITH SPECTRUM
        </a>
      </li>
      <li>
        <a
          href="/how-it-works"
          class={activePage === 'how-it-works' ? 'active' : ''}
          data-sveltekit-preload-data="hover"
        >
          HOW IT WORKS
        </a>
      </li>
      <li>
        <a
          href="/about"
          class={activePage === 'about' ? 'active' : ''}
          data-sveltekit-preload-data="hover"
        >
          ABOUT
        </a>
      </li>
    </ul>

    <!-- Mobile Menu Button -->
    <button class="mobile-menu-btn" on:click={toggleMobileMenu} aria-label="Toggle menu">
      <MenuIcon size={24} color="var(--text-secondary)" />
    </button>
  </div>
</nav>

<!-- Mobile Menu Backdrop -->
{#if isMobileMenuOpen}
  <div
    class="mobile-menu-backdrop"
    on:click={closeMobileMenu}
    on:keydown={(e) => e.key === 'Escape' && closeMobileMenu()}
    tabindex="-1"
    role="button"
  ></div>
{/if}

<!-- Mobile Slide-out Menu -->
<div class="mobile-menu" class:open={isMobileMenuOpen}>
  <div class="mobile-menu-header">
    <div class="logo">
      <span class="text-gradient prosto-one-regular">SMILES2SPEC</span>
    </div>
    <button class="mobile-menu-close" on:click={closeMobileMenu} aria-label="Close menu">
      <CloseIcon size={24} color="var(--text-secondary)" />
    </button>
  </div>

  <ul class="mobile-nav-links">
    <li>
      <a
        href="/"
        class={activePage === 'home' ? 'active' : ''}
        data-sveltekit-preload-data="hover"
        on:click={closeMobileMenu}
      >
        HOME
      </a>
    </li>
    <li>
      <a
        href="/spectral-simulation"
        class={activePage === 'spectral-simulation' ? 'active' : ''}
        data-sveltekit-preload-data="hover"
        on:click={closeMobileMenu}
      >
        SPECTRAL SIMULATION
      </a>
    </li>
    <li>
      <a
        href="/chat-with-spectrum"
        class={activePage === 'chat-with-spectrum' ? 'active' : ''}
        data-sveltekit-preload-data="hover"
        on:click={closeMobileMenu}
      >
        CHAT WITH SPECTRUM
      </a>
    </li>
    <li>
      <a
        href="/how-it-works"
        class={activePage === 'how-it-works' ? 'active' : ''}
        data-sveltekit-preload-data="hover"
        on:click={closeMobileMenu}
      >
        HOW IT WORKS
      </a>
    </li>
    <li>
      <a
        href="/about"
        class={activePage === 'about' ? 'active' : ''}
        data-sveltekit-preload-data="hover"
        on:click={closeMobileMenu}
      >
        ABOUT
      </a>
    </li>
  </ul>
</div>

<style>
  .navbar {
    width: 100%;
    max-width: 90%;
    margin: 1.5rem auto 1rem;
    padding: 0.75rem 1.75rem;
    border-radius: var(--enforce-pill);
    display: flex;
    justify-content: center;
    transition: all var(--transition-smooth);
    backdrop-filter: blur(30px) saturate(160%);
    -webkit-backdrop-filter: blur(30px) saturate(160%);
    background: rgba(255, 255, 255, 0.35);
    box-shadow:
      0 25px 55px rgba(0, 0, 0, 0.08),
      0 3px 12px rgba(0, 0, 0, 0.04);
    z-index: 10;
  }

  .navbar-container {
    width: 100%;
    display: flex;
    justify-content: space-between;
    align-items: center;
  }

  .logo {
    font-size: 1.25rem;
    font-weight: 600;
  }

  /* Desktop Navigation */
  .desktop-nav {
    display: flex;
    list-style: none;
    margin: 0;
    padding: 0;
    gap: 1.5rem;
  }

  .nav-links a {
    background: transparent;
    border: none;
    color: var(--text-secondary);
    font-family:
      'SF Pro Display',
      system-ui,
      -apple-system,
      BlinkMacSystemFont,
      'Segoe UI',
      Roboto,
      sans-serif;
    font-size: 1rem;
    font-weight: 500;
    cursor: pointer;
    padding: 0.5rem 1rem;
    border-radius: var(--enforce-pill);
    transition: all var(--transition-smooth);
    text-decoration: none;
    display: block;
  }

  .nav-links a:hover {
    color: var(--accent);
    background: rgba(120, 121, 255, 0.08);
  }

  .nav-links a.active {
    color: white;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    box-shadow: 0 4px 12px rgba(120, 121, 255, 0.2);
  }

  /* Mobile Menu Button */
  .mobile-menu-btn {
    display: none;
    background: transparent;
    border: none;
    cursor: pointer;
    padding: 0.5rem;
    border-radius: 8px;
    transition: background-color var(--transition-smooth);
  }

  .mobile-menu-btn:hover {
    background: rgba(120, 121, 255, 0.08);
  }

  /* Mobile Menu Backdrop */
  .mobile-menu-backdrop {
    position: fixed;
    top: 0;
    left: 0;
    width: 100vw;
    height: 100vh;
    background: rgba(0, 0, 0, 0.5);
    z-index: 40;
    backdrop-filter: blur(2px);
    -webkit-backdrop-filter: blur(2px);
  }

  /* Mobile Slide-out Menu */
  .mobile-menu {
    position: fixed;
    top: 1rem;
    right: -320px;
    width: 300px;
    height: calc(100vh - 2rem);
    background: rgba(255, 255, 255, 0.35);
    backdrop-filter: blur(30px) saturate(160%);
    -webkit-backdrop-filter: blur(30px) saturate(160%);
    border-radius: var(--enforce-pill);
    box-shadow:
      0 25px 55px rgba(0, 0, 0, 0.08),
      0 3px 12px rgba(0, 0, 0, 0.04);
    transition: right 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    z-index: 50;
    overflow-y: auto;
  }

  .mobile-menu.open {
    right: 1rem;
  }

  .mobile-menu-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 1.5rem 1.75rem 1rem;
    border-bottom: 1px solid rgba(255, 255, 255, 0.2);
  }

  .mobile-menu-close {
    background: transparent;
    border: none;
    cursor: pointer;
    padding: 0.5rem;
    border-radius: 8px;
    transition: background-color var(--transition-smooth);
  }

  .mobile-menu-close:hover {
    background: rgba(120, 121, 255, 0.08);
  }

  .mobile-nav-links {
    list-style: none;
    margin: 0;
    padding: 1rem 0;
  }

  .mobile-nav-links li {
    margin: 0;
  }

  .mobile-nav-links a {
    display: block;
    padding: 1rem 1.5rem;
    color: var(--text-secondary);
    font-family:
      'SF Pro Display',
      system-ui,
      -apple-system,
      BlinkMacSystemFont,
      'Segoe UI',
      Roboto,
      sans-serif;
    font-size: 1rem;
    font-weight: 500;
    text-decoration: none;
    transition: all var(--transition-smooth);
    border-left: 3px solid transparent;
  }

  .mobile-nav-links a:hover {
    background: rgba(120, 121, 255, 0.08);
    color: var(--accent);
    border-left-color: var(--accent);
  }

  .mobile-nav-links a.active {
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    color: white;
    border-left-color: var(--accent-secondary);
  }

  /* Responsive Breakpoints */
  @media (max-width: 767px) {
    .desktop-nav {
      display: none;
    }

    .mobile-menu-btn {
      display: block;
    }

    .navbar {
      max-width: 95%;
      padding: 0.75rem 1.25rem;
      margin: 1rem auto 0.75rem;
    }

    .logo {
      font-size: 1.1rem;
    }
  }

  @media (min-width: 768px) {
    .mobile-menu-btn {
      display: none;
    }

    .mobile-menu {
      display: none;
    }

    .mobile-menu-backdrop {
      display: none;
    }
  }

  /* Small mobile adjustments */
  @media (max-width: 480px) {
    .mobile-menu {
      width: 280px;
      right: -300px;
    }

    .mobile-menu.open {
      right: 0.5rem;
    }

    .mobile-nav-links a {
      font-size: 0.95rem;
      padding: 0.9rem 1.25rem;
    }
  }

  /* Extra-large screen adjustments */
  @media (min-width: 1920px) {
    .navbar {
      max-width: 85%;
    }
  }
</style>
