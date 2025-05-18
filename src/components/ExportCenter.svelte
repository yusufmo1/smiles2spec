<script>
  import { exportMsp, exportMspBatch } from '../services/api.js';
  import Plotly from 'plotly.js-dist-min';
  import JSZip from 'jszip';
  import { predictSpectrum } from '../services/api.js';
  
  // Export center component for all export functionality
  export let spectrumData = null;
  export let peakData = [];
  export let structurePNG = null;
  export let smiles = "";
  export let chemicalName = "";
  export let smilesList = [];
  
  // Track available export options
  $: hasSpectrum = spectrumData !== null;
  $: hasPeaks = Array.isArray(peakData) && peakData.length > 0;
  $: hasStructure = structurePNG !== null;
  $: hasChemicalInfo = !!smiles || !!chemicalName;
  $: hasBatchData = Array.isArray(smilesList) && smilesList.length > 1;
  
  async function exportMspData() {
    try {
      let blob;
      
      if (hasBatchData) {
        // Export all entries as batch
        blob = await exportMspBatch(smilesList);
      } else if (smiles) {
        // Export single entry
        blob = await exportMsp(smiles);
      } else {
        throw new Error("No chemical data available to export");
      }
      
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = hasBatchData ? `batch_export.msp` : `${smiles || chemicalName || "export"}.msp`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      alert(`Could not export MSP: ${e.message}`);
    }
  }

  function downloadSpectrumPNG() {
    try {
      // Find the spectrum plot element
      const plotElement = document.querySelector('.plot-container');
      if (!plotElement) {
        throw new Error("Spectrum plot not found");
      }

      // Access the Plotly.js instance
      const fileName = `${smiles || chemicalName || "spectrum"}_mass_spectrum.png`;
      
      // Create custom layout with title including SMILES and white background
      const layout = {
        title: {
          text: `Mass Spectrum - ${smiles}`, 
          font: {
            size: 16,
            color: 'black'
          }
        },
        paper_bgcolor: 'white',
        plot_bgcolor: 'white',
        autosize: true
      };
      
      // Use Plotly's toImage function to generate a PNG with modified layout
      Plotly.toImage(plotElement, {
        format: 'png',
        width: 1200,
        height: 800,
        scale: 2, // Higher resolution
        layout: layout
      }).then(dataUrl => {
        // Create a link element and trigger download
        const a = document.createElement("a");
        a.href = dataUrl;
        a.download = fileName;
        a.click();
      }).catch(err => {
        throw new Error(`Failed to generate PNG: ${err.message}`);
      });
    } catch (e) {
      alert(`Could not download PNG: ${e.message}`);
    }
  }
  
  async function downloadBatchSpectrumPNGs() {
    try {
      if (!hasBatchData || !smilesList || smilesList.length === 0) {
        throw new Error("No batch data available to export");
      }
      
      // Create a zip file
      const zip = new JSZip();
      let completed = 0;
      const total = smilesList.length;
      
      // Use progress indicator
      const progressIndicator = document.createElement("div");
      progressIndicator.style.position = "fixed";
      progressIndicator.style.top = "10px";
      progressIndicator.style.left = "50%";
      progressIndicator.style.transform = "translateX(-50%)";
      progressIndicator.style.padding = "10px 20px";
      progressIndicator.style.background = "rgba(0,0,0,0.7)";
      progressIndicator.style.color = "white";
      progressIndicator.style.borderRadius = "4px";
      progressIndicator.style.zIndex = "9999";
      progressIndicator.textContent = `Preparing PNGs: 0/${total}`;
      document.body.appendChild(progressIndicator);
      
      // Process each SMILES string to generate a spectrum PNG
      for (const [index, currentSmiles] of smilesList.entries()) {
        try {
          // Predict spectrum for this SMILES
          const result = await predictSpectrum(currentSmiles);
          
          if (result && result.spectrum) {
            // Create a temporary plot element
            const tempDiv = document.createElement("div");
            document.body.appendChild(tempDiv);
            tempDiv.style.width = "1200px";
            tempDiv.style.height = "800px";
            tempDiv.style.position = "absolute";
            tempDiv.style.visibility = "hidden";
            
            // Generate plot data
            const filteredData = [];
            for (let i = 0; i < result.spectrum.x.length; i++) {
              if (result.spectrum.y[i] > 0.01) {
                filteredData.push({
                  x: result.spectrum.x[i],
                  y: result.spectrum.y[i]
                });
              }
            }
            
            const plotData = [{
              x: filteredData.map(d => d.x),
              y: filteredData.map(d => d.y),
              type: 'bar',
              marker: {
                color: 'var(--accent)',
                line: {
                  width: 1,
                  color: 'rgba(0,0,0,0.05)'
                }
              }
            }];
            
            // Create custom layout with title including SMILES and white background
            const layout = {
              title: {
                text: `Mass Spectrum - ${currentSmiles}`,
                font: {
                  size: 16,
                  color: 'black'
                }
              },
              xaxis: {
                title: {
                  text: 'M/Z',
                  font: {
                    size: 12,
                    color: 'black'
                  }
                },
                range: [0, 250]
              },
              yaxis: {
                title: {
                  text: 'Relative Intensity',
                  font: {
                    size: 12,
                    color: 'black'
                  }
                }
              },
              paper_bgcolor: 'white',
              plot_bgcolor: 'white',
              autosize: true,
              margin: {
                l: 60,
                r: 40,
                b: 50,
                t: 50,
                pad: 10
              }
            };
            
            // Create the plot
            await Plotly.newPlot(tempDiv, plotData, layout, {displayModeBar: false});
            
            // Convert plot to image
            const dataUrl = await Plotly.toImage(tempDiv, {
              format: 'png',
              width: 1200,
              height: 800,
              scale: 2
            });
            
            // Convert data URL to blob
            const response = await fetch(dataUrl);
            const blob = await response.blob();
            
            // Add to zip file
            const chemName = result.chemical_name || `compound_${index+1}`;
            const safeSmiles = currentSmiles.replace(/[\/\\:*?"<>|]/g, '_').substring(0, 30);
            zip.file(`${safeSmiles}_${chemName}_mass_spectrum.png`, blob);
            
            // Clean up temp element
            document.body.removeChild(tempDiv);
            
            // Update progress
            completed++;
            progressIndicator.textContent = `Preparing PNGs: ${completed}/${total}`;
          }
        } catch (innerError) {
          console.error(`Error processing SMILES ${currentSmiles}:`, innerError);
          // Continue with next SMILES even if this one fails
        }
      }
      
      // Generate zip file
      const zipBlob = await zip.generateAsync({type: 'blob'});
      
      // Download zip file
      const url = URL.createObjectURL(zipBlob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "mass_spectra_batch.zip";
      a.click();
      URL.revokeObjectURL(url);
      
      // Remove progress indicator
      document.body.removeChild(progressIndicator);
    } catch (e) {
      console.error("Error in batch download:", e);
      alert(`Could not download batch PNGs: ${e.message}`);
    }
  }
</script>

<div class="export-center">
  <div class="export-status">
    {#if hasChemicalInfo}
      <span class="status-item available">✓ Chemical data</span>
    {:else}
      <span class="status-item unavailable">✗ No chemical data</span>
    {/if}
    
    {#if hasSpectrum}
      <span class="status-item available">✓ Spectrum data</span>
    {:else}
      <span class="status-item unavailable">✗ No spectrum data</span>
    {/if}
    
    {#if hasStructure}
      <span class="status-item available">✓ Structure image</span>
    {:else}
      <span class="status-item unavailable">✗ No structure image</span>
    {/if}
  </div>
  
  <div class="export-options">
    <button class="export-btn" 
            on:click={exportMspData} 
            disabled={!smiles && (!smilesList || smilesList.length === 0)} 
            title="Export as MSP format">
      <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
        <path d="M4 18h16v2H4v-2zM4 14h16v2H4v-2zM8 4v6h3v6h2v-6h3V4H8z"/>
      </svg>
      {hasBatchData ? 'Export All MSP' : 'Export MSP'}
    </button>

    <button class="export-btn" 
            on:click={downloadSpectrumPNG} 
            disabled={!hasSpectrum} 
            title="Download mass spectrum as PNG">
      <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
        <path d="M5 20h14v-2H5v2zm7-13l5 5h-3v4h-4v-4H7l5-5z"/>
      </svg>
      Download Spectrum PNG
    </button>
    
    {#if hasBatchData}
      <button class="export-btn" 
              on:click={downloadBatchSpectrumPNGs} 
              title="Download all mass spectra as PNG files in a ZIP">
        <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
          <path d="M20 6h-8l-2-2H4c-1.1 0-1.99.9-1.99 2L2 18c0 1.1.9 2 2 2h16c1.1 0 2-.9 2-2V8c0-1.1-.9-2-2-2zm0 12H4V8h16v10zM8 13.01l1.41 1.41L11 12.84V17h2v-4.16l1.59 1.59L16 13.01 12.01 9 8 13.01z"/>
        </svg>
        Export All Spectra as PNGs
      </button>
    {/if}
  </div>
</div>

<style>
  .export-center {
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
    gap: 1.5rem;
  }
  
  .export-status {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    justify-content: center;
    font-size: 0.75rem;
  }
  
  .status-item {
    padding: 0.35rem 0.75rem;
    border-radius: var(--enforce-pill);
  }
  
  .available {
    background: rgba(52, 199, 89, 0.15);
    color: #34c759;
  }
  
  .unavailable {
    background: rgba(142, 142, 147, 0.15);
    color: #8e8e93;
  }
  
  .export-options {
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
    justify-content: center;
  }
  
  .export-btn {
    background: var(--accent-soft);
    border: none;
    color: var(--accent);
    font-size: 0.8rem;
    font-weight: 500;
    padding: 0.5rem 1rem;
    border-radius: var(--enforce-pill);
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    transition: all 0.2s;
  }

  .export-btn:hover:not(:disabled) {
    background: var(--accent);
    color: white;
  }

  .export-btn:disabled {
    opacity: 0.4;
    cursor: not-allowed;
  }
</style>