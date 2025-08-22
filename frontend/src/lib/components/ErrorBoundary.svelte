<script lang="ts">
  import { onMount, createEventDispatcher } from 'svelte';
  import { DynamicIcon } from '$lib/components/icons';

  const dispatch = createEventDispatcher();

  export const fallback = null;
  let error: any = null;
  let hasError = false;
  let errorDetails: any = '';
  let showDetails = false;

  function handleError(event: any) {
    console.error('Application Error:', event.error || event);

    // Add error reporting service integration here if needed
    // reportError(error, errorDetails);

    error = event.error || event;
    hasError = true;

    // Capture error details
    errorDetails = {
      message: error?.message || 'Unknown error',
      stack: error?.stack || 'No stack trace available',
      url: window.location.href,
      userAgent: navigator.userAgent,
      timestamp: new Date().toISOString(),
      componentStack: error?.componentStack || null,
      errorBoundary: 'ErrorBoundary.svelte',
    };

    // Dispatch error event
    dispatch('error', { error, details: errorDetails });

    // Prevent default error handling
    if (event.preventDefault) event.preventDefault();
  }

  function retry() {
    hasError = false;
    error = null;
    errorDetails = '';
    window.location.reload();
  }

  function goHome() {
    hasError = false;
    error = null;
    errorDetails = '';
    window.location.href = '/';
  }

  function toggleDetails() {
    showDetails = !showDetails;
  }

  function copyErrorDetails() {
    const errorInfo = JSON.stringify(errorDetails, null, 2);
    navigator.clipboard.writeText(errorInfo).then(() => {
      alert('Error details copied to clipboard');
    });
  }

  onMount(() => {
    // Handle unhandled promise rejections
    const handleUnhandledRejection = (event: any) => {
      handleError(event);
    };

    // Handle JavaScript errors
    const handleWindowError = (event: any) => {
      handleError(event);
    };

    window.addEventListener('unhandledrejection', handleUnhandledRejection);
    window.addEventListener('error', handleWindowError);

    return () => {
      window.removeEventListener('unhandledrejection', handleUnhandledRejection);
      window.removeEventListener('error', handleWindowError);
    };
  });
</script>

{#if hasError}
  <div class="error-boundary" role="alert">
    <div class="error-content">
      <div class="error-icon">
        <DynamicIcon name="ErrorIcon" size={48} color="#ff453a" />
      </div>
      <h2>Oops! Something went wrong</h2>
      <p class="error-message">
        {error?.message || 'An unexpected error occurred while navigating the application.'}
      </p>

      <div class="error-actions">
        <button class="primary-button" on:click={retry}>
          <DynamicIcon name="LoadingIcon" size={16} />
          Retry
        </button>
        <button class="secondary-button" on:click={goHome}>
          <DynamicIcon name="HomeIcon" size={16} />
          Go Home
        </button>
        <button class="detail-button" on:click={toggleDetails}>
          {#if showDetails}
            <DynamicIcon name="BookIcon" size={16} />
            Hide Details
          {:else}
            <DynamicIcon name="SearchIcon" size={16} />
            Show Details
          {/if}
        </button>
      </div>

      {#if showDetails}
        <div class="error-details">
          <h3>Error Details</h3>
          <div class="error-info">
            <p><strong>URL:</strong> {errorDetails.url}</p>
            <p><strong>Time:</strong> {errorDetails.timestamp}</p>
            <p><strong>Message:</strong> {errorDetails.message}</p>
          </div>

          {#if errorDetails.stack}
            <div class="stack-trace">
              <h4>Stack Trace:</h4>
              <pre>{errorDetails.stack}</pre>
            </div>
          {/if}

          <button class="copy-button" on:click={copyErrorDetails}>
            <DynamicIcon name="CopyIcon" size={16} />
            Copy Error Details
          </button>
        </div>
      {/if}
    </div>
  </div>
{:else}
  <slot />
{/if}

<style>
  .error-boundary {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: rgba(255, 255, 255, 0.95);
    backdrop-filter: blur(10px);
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 10000;
    padding: 2rem;
  }

  .error-content {
    background: white;
    padding: 3rem;
    border-radius: 24px;
    box-shadow: 0 25px 55px rgba(0, 0, 0, 0.15);
    text-align: center;
    max-width: 600px;
    width: 100%;
    border: 1px solid rgba(255, 75, 85, 0.2);
  }

  .error-icon {
    font-size: 4rem;
    margin-bottom: 1rem;
  }

  .error-content h2 {
    color: #333;
    margin: 0 0 1rem;
    font-size: 1.5rem;
  }

  .error-message {
    color: #666;
    margin-bottom: 2rem;
    line-height: 1.5;
  }

  .error-actions {
    display: flex;
    gap: 1rem;
    justify-content: center;
    flex-wrap: wrap;
    margin-bottom: 2rem;
  }

  .primary-button {
    background: #7879ff;
    color: white;
    border: none;
    padding: 0.75rem 1.5rem;
    border-radius: 24px;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s;
  }

  .primary-button:hover {
    background: #6969ee;
    transform: translateY(-2px);
  }

  .secondary-button,
  .detail-button,
  .copy-button {
    background: transparent;
    color: #666;
    border: 1px solid #666;
    padding: 0.75rem 1.5rem;
    border-radius: 24px;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s;
  }

  .secondary-button:hover,
  .detail-button:hover,
  .copy-button:hover {
    background: #666;
    color: white;
  }

  .error-details {
    background: rgba(0, 0, 0, 0.05);
    padding: 1.5rem;
    border-radius: 12px;
    text-align: left;
    margin-top: 1rem;
  }

  .error-details h3,
  .error-details h4 {
    margin: 0 0 1rem;
    color: #333;
  }

  .error-info p {
    margin: 0.5rem 0;
    font-size: 0.9rem;
    word-break: break-all;
  }

  .stack-trace {
    margin-top: 1rem;
  }

  .stack-trace pre {
    background: rgba(0, 0, 0, 0.1);
    padding: 1rem;
    border-radius: 8px;
    font-size: 0.8rem;
    overflow-x: auto;
    white-space: pre-wrap;
    word-break: break-all;
  }

  .copy-button {
    margin-top: 1rem;
    font-size: 0.9rem;
  }

  @media (max-width: 640px) {
    .error-boundary {
      padding: 1rem;
    }

    .error-content {
      padding: 2rem 1.5rem;
    }

    .error-actions {
      flex-direction: column;
      align-items: center;
    }

    .primary-button,
    .secondary-button,
    .detail-button {
      width: 100%;
      max-width: 200px;
    }
  }
</style>
