<script>
	import { onMount, tick } from 'svelte';
	import Header from './components/Header.svelte';
	import SmilesInput from './components/SmilesInput.svelte';
	import SpectrumPlot from './components/SpectrumPlot.svelte';
	import IonFragmentation from './components/IonFragmentation.svelte';
	import PeakTable from './components/PeakTable.svelte';
	import ConsoleOutput from './components/ConsoleOutput.svelte';
	import Panel from './components/Panel.svelte';
	import Background from './components/Background.svelte';
	import StructurePanel from './components/StructurePanel.svelte';
	import BulkNavigator from './components/BulkNavigator.svelte';
	import PanelOverlay from './components/PanelOverlay.svelte';
	import ExportCenter from './components/ExportCenter.svelte';
	import Navbar from './components/Navbar.svelte';
	import HowItWorks from './components/HowItWorks.svelte';
	import ChatWithSpectrum from './components/ChatWithSpectrum.svelte';
	import { predictSpectrum } from './services/api.js';
	import { isCarouselMode, focusedPanel, carouselIndex } from './stores.js';
	import { get } from 'svelte/store';
	
	// Add responsive design variables
	let windowWidth;
	
	// Add page navigation state
	let currentPage = 'home';
	
	let spectrumData = null;
	let peakData = [];
	let smilesInputComponent;
	let error = null;
	let consoleText = "";
	let structurePNG = null;
	let bulkList = [];          // all SMILES from txt
	let idx = 0;                // current pointer
	let currentName = "";       // from prediction response
	let currentSmiles = "";
	let hasFirstPrediction = false;  // tracking if we have at least one successful prediction
	
	// No default banner - the console starts empty
	
	async function handlePredict(event) {
		const { smiles } = event.detail;
		currentSmiles = smiles;
		error = null;
		const startTime = performance.now();
		const timestamp = new Date().toISOString().split('T')[1].split('.')[0];
		consoleText += `\n[${timestamp}] > PROCESSING SMILES: ${smiles}\n`;
		
		try {
			// Log API request start
			const apiStartTime = performance.now();
			consoleText += `[${timestamp}] SENDING API REQUEST...\n`;
			
			// Make API request
			const result = await predictSpectrum(smiles);
			
			// Log API response time
			const apiEndTime = performance.now();
			const apiDuration = Math.round(apiEndTime - apiStartTime);
			consoleText += `[${timestamp}] API RESPONSE RECEIVED (${apiDuration}ms)\n`;
			
			// Store the data
			spectrumData = result.spectrum;
			peakData = result.peaks;
			structurePNG = result.structure_png;
			currentName = result.chemical_name;
			hasFirstPrediction = true; // Mark that we have at least one successful prediction
			
			// Calculate and report peak statistics
			const peakCount = peakData.length;
			const significantPeaks = peakData.filter(p => p.intensity > 0.05).length;
			const maxIntensity = Math.max(...peakData.map(p => p.intensity));
			const baselineNoise = Math.min(...peakData.filter(p => p.intensity > 0).map(p => p.intensity));
			
			// Detailed chemical information
			consoleText += `[${timestamp}] MOLECULAR FORMULA: C${Math.floor(result.molecular_weight/12)}H${Math.floor(result.molecular_weight/2)}\n`;
			consoleText += `[${timestamp}] MOLECULAR WEIGHT: ${result.molecular_weight.toFixed(4)} g/mol\n`;
			consoleText += `[${timestamp}] EXACT MASS: ${result.exact_mass.toFixed(4)} amu\n`;
			
			// Peak analysis
			consoleText += `[${timestamp}] SPECTRAL ANALYSIS RESULTS:\n`;
			consoleText += `[${timestamp}] - TOTAL PEAKS: ${peakCount}\n`;
			consoleText += `[${timestamp}] - SIGNIFICANT PEAKS (>5%): ${significantPeaks}\n`;
			consoleText += `[${timestamp}] - BASE PEAK INTENSITY: ${maxIntensity.toFixed(4)}\n`;
			consoleText += `[${timestamp}] - SIGNAL-TO-NOISE: ${(maxIntensity/baselineNoise).toFixed(2)}\n`;
			
			// Top peaks listing
			const topPeaks = [...peakData]
				.sort((a, b) => b.intensity - a.intensity)
				.slice(0, 5);
			
			consoleText += `[${timestamp}] TOP 5 FRAGMENT PEAKS:\n`;
			topPeaks.forEach((peak, idx) => {
				consoleText += `[${timestamp}] ${idx+1}. M/Z ${peak.mz.toFixed(1)}: ${peak.intensity.toFixed(4)} (${Math.round(peak.intensity*100)}%)\n`;
			});
			
			// Final timing information
			const endTime = performance.now();
			const totalDuration = Math.round(endTime - startTime);
			consoleText += `[${timestamp}] ANALYSIS COMPLETED IN ${totalDuration}ms\n`;
		} catch (err) {
			const endTime = performance.now();
			const errorDuration = Math.round(endTime - startTime);
			const timestamp = new Date().toISOString().split('T')[1].split('.')[0];
			
			error = `ERROR: ${err.message}`;
			consoleText += `[${timestamp}] ERROR OCCURRED AFTER ${errorDuration}ms: ${err.message}\n`;
			consoleText += `[${timestamp}] ANALYSIS FAILED\n`;
			consoleText += `[${timestamp}] > RETRY WITH VALID SMILES\n`;
			console.error(err);
		} finally {
			smilesInputComponent.setLoading(false);
		}
	}
	
	// handle bulk list coming up from SmilesInput
	function handleBulk(event) {
		bulkList = event.detail.list;
		idx = 0;
		smilesInputComponent.setLoading(true);
		currentSmiles = bulkList[0];  // Set the current SMILES for display
		handlePredict({ detail: { smiles: bulkList[0] } });
	}
	
	// navigator arrows
	function nav(offset) {
		if (!bulkList.length) return;
		idx = (idx + offset + bulkList.length) % bulkList.length;
		const s = bulkList[idx];
		currentSmiles = s;  // Update the current SMILES
		handlePredict({ detail: { smiles: s } });
	}
	
	// Handle carousel toggle
	async function handleToggle(val) {
		isCarouselMode.set(val);
		
		if (val) {
			// wait for next tick to ensure DOM is updated
			await tick();
			
			// If no panel is focused yet, select the first one
			if (!$focusedPanel) {
				const firstPanel = document.querySelector('.panel');
				if (firstPanel) {
					focusedPanel.set(firstPanel);
					
					// Also set the carousel index to 0 for the first panel
					carouselIndex.set(0);
				}
			}
		} else {
			// close overlay
			focusedPanel.set(null);
		}
	}
	
	// Handle page navigation from the navbar
	function handlePageChange(newPage) {
		currentPage = newPage;
	}
</script>

<svelte:window bind:innerWidth={windowWidth}/>

<svelte:head>
  <link href="https://fonts.googleapis.com/css2?family=SF+Pro+Display:wght@400;500;600&family=SF+Pro+Text:wght@400;500;600&family=SF+Mono:wght@400;600&display=swap" rel="stylesheet">
  <link href="https://fonts.googleapis.com/css2?family=Prosto+One&display=swap" rel="stylesheet">
  <meta name="theme-color" content="#f4f5f7">
</svelte:head>

<Background />

<!-- New Navbar moved here -->
<Navbar 
  activePage={currentPage}
  on:click={(e) => handlePageChange(e.detail)}
/>

<!-- NEW master pill -->
<div class="app-shell glass-card">
  {#if currentPage === 'home'}
    <Header 
      className="prosto-one-regular" 
      showToggle={hasFirstPrediction}
      toggleValue={$isCarouselMode}
      onToggle={handleToggle}
    />
    
    <div class="body-wrapper">
      <!-- SMILES Input on its own row -->
      <div class="row">
        <div class="col-full">
          <SmilesInput 
            on:predict={handlePredict}
            on:bulk={handleBulk}
            bind:this={smilesInputComponent}
          />
          
          <!-- NEW navigator pill right below -->
          <BulkNavigator name={currentName}
                         index={bulkList.length ? idx+1 : 0}
                         total={bulkList.length}
                         on:prev={() => nav(-1)}
                         on:next={() => nav(1)}/>
        </div>
      </div>
      
      <!-- Error message if any -->
      {#if error}
        <div class="row">
          <div class="col-full">
            <div class="error-message">
              <span class="error-icon">⚠️</span>
              {error}
            </div>
          </div>
        </div>
      {/if}
      
      <!-- Main visualization section with two equal columns -->
      <div class="row">
        <div class="col-half">
          <Panel title="MASS SPECTRUM">
            <SpectrumPlot spectrumData={spectrumData} />
          </Panel>
        </div>
        <div class="col-half">
          <Panel title="TOP FRAGMENT IONS">
            <IonFragmentation data={spectrumData} />
          </Panel>
        </div>
      </div>
      
      <!-- Chemical structure and fragment ions row -->
      <div class="row">
        <div class="col-half">
          <Panel title="CHEMICAL STRUCTURE">
            <StructurePanel png={structurePNG} />
          </Panel>
        </div>
        <div class="col-half">
          <Panel title="PEAK DATA TABLE">
            <PeakTable peaks={peakData} smiles={currentSmiles} />
          </Panel>
        </div>
      </div>
      
      <!-- Bottom row with Console Output and Chat -->
      <div class="row">
        <div class="col-half">
          <Panel title="ANALYSIS CONSOLE">
            <ConsoleOutput output={consoleText} />
          </Panel>
        </div>
        <div class="col-half">
          <Panel title="CHAT WITH SPECTRUM">
            <ChatWithSpectrum 
              hasSmilesPrediction={hasFirstPrediction}
              currentSmiles={currentSmiles} 
            />
          </Panel>
        </div>
      </div>
      
      <!-- Export Center row -->
      <div class="row">
        <div class="col-half">
          <Panel title="EXPORT CENTER">
            <ExportCenter 
              spectrumData={spectrumData}
              peakData={peakData}
              structurePNG={structurePNG}
              smiles={currentSmiles}
              chemicalName={currentName}
              smilesList={bulkList}
            />
          </Panel>
        </div>
        <div class="col-half">
          <!-- Intentionally left blank for UI consistency -->
        </div>
      </div>
    </div>
  {:else if currentPage === 'how-it-works'}
    <Header 
      className="prosto-one-regular" 
      showToggle={false}
    />
    
    <HowItWorks />
  {/if}
</div>

<!-- Moved PanelOverlay outside the app-shell -->
<PanelOverlay />

<style>
	/* Removed main styles */
	
	/* Outer pill styling */
	.app-shell {
		width: 100%;
		max-width: 90%;           /* Changed from 1440px to 90% to use most of the screen width */
		margin: 0 auto 3rem;      /* Adjusted margin-top to 0 to account for the navbar */
		padding: 3.5rem 3rem 4rem; /* room for header and panels */
		border-radius: var(--enforce-pill);
		background: rgba(255,255,255,0.35);   /* much more transparent white */
		box-shadow:
			0 25px 55px rgba(0,0,0,0.08),
			0 3px 12px rgba(0,0,0,0.04);
		backdrop-filter: blur(30px) saturate(160%);
		-webkit-backdrop-filter: blur(30px) saturate(160%);
		display: flex;
		flex-direction: column;
		gap: 2.5rem;               /* space between header and body */
		overflow: visible;         /* keep drop-shadows of inner panels */
	}

	/* Container that holds rows/cols */
	.body-wrapper {
		display: flex;
		flex-direction: column;
		gap: 2rem;                 /* vertical rhythm same as before */
	}
	
	/* Grid system */
	.row {
		display: flex;
		flex-wrap: wrap;
		margin: 0;                /* removed side margins as we're inside a card */
		margin-bottom: 0;         /* gap handles spacing now */
	}
	
	.col-half {
		width: 50%;
		padding: 0 1.5rem;        /* Increased padding from 1rem to 1.5rem for better spacing */
		margin-bottom: 0;
	}
	
	.col-full {
		width: 100%;
		padding: 0 1.5rem;        /* Increased padding from 1rem to 1.5rem for consistency */
	}
	
	/* Error message */
	.error-message {
		background: rgba(255, 75, 85, 0.1);
		color: #ff453a;
		padding: 1rem 1.5rem;
		margin-bottom: 1rem;
		font-size: 0.9rem;
		font-weight: 500;
		display: flex;
		align-items: center;
		gap: 0.75rem;
		border-radius: var(--enforce-pill);
		border: 1px solid rgba(255, 75, 85, 0.25);
		box-shadow: 0 4px 12px rgba(255, 75, 85, 0.1);
		overflow: hidden !important;
	}
	
	.error-icon {
		font-size: 1.2rem;
	}
	
	/* Responsive adjustments */
	@media (max-width: 1024px) {
		.col-half {
			width: 100%;
			margin-bottom: 2rem;
		}
		
		.row .col-half:last-child {
			margin-bottom: 0;
		}
		
		.row {
			margin-bottom: 0;
		}
		
		.col-half, .col-full {
			padding: 0 1rem;      /* Keep 1rem padding on smaller screens */
		}
	}
	
	/* Additional responsive tweaks for small screens */
	@media (max-width: 640px) {
		.app-shell {
			margin: 0 auto 1rem;  /* Adjusted top margin for navbar */
			max-width: 95%;       /* Use even more width on small screens */
			padding: 2rem 1.25rem 2.5rem;
		}
	}

	/* Extra-large screen adjustments */
	@media (min-width: 1920px) {
		.app-shell {
			max-width: 85%;      /* Slightly smaller percentage on very large screens */
		}
		
		.col-half {
			padding: 0 2rem;     /* More padding between panels on very large screens */
		}
		
		.col-full {
			padding: 0 2rem;
		}
	}
	
	/* Ensure overlay can break out of any container */
	:global(#app), :global(body), :global(html) {
		overflow-x: hidden; /* Prevent horizontal scrolling */
		position: relative;
	}
	
	:global(body), :global(html) {
		min-height: 100vh;
		width: 100%;
		margin: 0;
		padding: 0;
	}
</style>