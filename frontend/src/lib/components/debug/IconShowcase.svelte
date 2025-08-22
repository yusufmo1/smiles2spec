<script lang="ts">
  import {
    DynamicIcon,
    IconSizes,
    IconColors,
    type IconSize,
    type IconColor,
    type IconName,
  } from '$lib/components/icons';
  import { getIconByType, getIconAriaLabel, type IconType } from '$lib/utils/iconMapping';

  // Options for displaying icons
  let selectedSize: IconSize = 'md';
  let selectedColor: IconColor = 'primary';
  let showSemanticMapping = false;

  // Group icons by category
  const iconCategories: Record<string, IconName[]> = {
    Actions: [
      'CloseIcon',
      'ArrowLeftIcon',
      'ArrowRightIcon',
      'SearchIcon',
      'UploadIcon',
      'DownloadIcon',
      'SendIcon',
      'CopyIcon',
      'ExportIcon',
      'RefreshIcon',
    ],
    Layout: ['GridIcon', 'CarouselIcon'],
    Science: [
      'MoleculeIcon',
      'SpectrumIcon',
      'CpuIcon',
      'ChartIcon',
      'DatabaseIcon',
      'FragmentationIcon',
      'DataTableIcon',
      'FlaskIcon',
      'MicroscopeIcon',
      'AtomIcon',
      'BondIcon',
      'TestTubeIcon',
      'BeakerIcon',
    ],
    Technologies: [
      'PythonIcon',
      'NPMIcon',
      'RDKitIcon',
      'ScikitLearnIcon',
      'DockerIcon',
      'PlotlyJSIcon',
      'MatplotlibIcon',
      'PyTorchIcon',
      'TensorFlowIcon',
      'SvelteKitIcon',
      'TypeScriptIcon',
      'CSSIcon',
    ],
    Status: ['CheckIcon', 'ErrorIcon', 'InfoIcon', 'LoadingIcon'],
    Features: [
      'WandIcon',
      'ChatWithSpectrumIcon',
      'LightbulbIcon',
      'TargetIcon',
      'FlashIcon',
      'ScaleIcon',
      'ValidationIcon',
      'SettingsIcon',
    ],
    'Social & UI': [
      'UserIcon',
      'GraduationCapIcon',
      'BookIcon',
      'HeartIcon',
      'StarIcon',
      'HomeIcon',
      'PhoneIcon',
      'MailIcon',
      'TwitterIcon',
      'LinkedinIcon',
      'GithubIcon',
    ],
  };

  // Semantic mappings for documentation
  const semanticMappings: Array<{ type: IconType; name: string; description: string }> = [
    { type: 'error', name: 'ErrorIcon', description: 'Used for errors and failures' },
    { type: 'success', name: 'CheckIcon', description: 'Used for success states' },
    { type: 'info', name: 'InfoIcon', description: 'Used for informational notes' },
    { type: 'loading', name: 'LoadingIcon', description: 'Used for loading states' },
    { type: 'cpu', name: 'CpuIcon', description: 'Used for processing and computation' },
    { type: 'chart', name: 'SpectrumIcon', description: 'Used for data visualization' },
    { type: 'generate', name: 'WandIcon', description: 'Used for generation actions' },
    { type: 'flask', name: 'FlaskIcon', description: 'Used for experiments' },
    { type: 'microscope', name: 'MicroscopeIcon', description: 'Used for detailed analysis' },
  ];
</script>

<div class="icon-showcase">
  <div class="controls">
    <div class="control-group">
      <label for="icon-size">Size:</label>
      <select id="icon-size" bind:value={selectedSize}>
        {#each Object.keys(IconSizes) as size}
          <option value={size}>{size} ({IconSizes[size as IconSize]}px)</option>
        {/each}
      </select>
    </div>

    <div class="control-group">
      <label for="icon-color">Color:</label>
      <select id="icon-color" bind:value={selectedColor}>
        {#each Object.keys(IconColors) as color}
          <option value={color}>{color}</option>
        {/each}
      </select>
    </div>

    <div class="control-group">
      <label>
        <input type="checkbox" bind:checked={showSemanticMapping} />
        Show semantic mapping
      </label>
    </div>
  </div>

  {#if showSemanticMapping}
    <div class="semantic-mapping">
      <h2>Semantic Icon Mapping</h2>
      <p>For consistent icon usage, use the semantic mapping utility:</p>

      <pre>import &#123; getIconByType &#125; from '$lib/utils/iconMapping';</pre>

      <div class="mapping-examples">
        {#each semanticMappings as mapping}
          <div class="mapping-example">
            <div class="mapping-icon">
              <DynamicIcon
                name={getIconByType(mapping.type)}
                size={IconSizes[selectedSize]}
                color={IconColors[selectedColor]}
              />
            </div>
            <div class="mapping-details">
              <div class="mapping-type">{mapping.type}</div>
              <div class="mapping-name">{mapping.name}</div>
              <div class="mapping-desc">{mapping.description}</div>
            </div>
          </div>
        {/each}
      </div>

      <pre>// Example usage:
const SuccessIcon = getIconByType('success');
const ariaLabel = getIconAriaLabel('success'); // Returns "Success"</pre>
    </div>
  {:else}
    {#each Object.entries(iconCategories) as [category, iconNames]}
      <div class="icon-category">
        <h2>{category}</h2>
        <div class="icon-grid">
          {#each iconNames as iconName}
            <div class="icon-card">
              <div class="icon-preview">
                <DynamicIcon
                  name={iconName}
                  size={IconSizes[selectedSize]}
                  color={IconColors[selectedColor]}
                />
              </div>
              <div class="icon-name">{iconName}</div>
              <div class="icon-usage">
                <code
                  >{`<DynamicIcon name="${iconName}" size={IconSizes.${selectedSize}} color={IconColors.${selectedColor}} />`}</code
                >
              </div>
            </div>
          {/each}
        </div>
      </div>
    {/each}
  {/if}
</div>

<style>
  .icon-showcase {
    padding: 1rem;
    color: var(--text-primary);
  }

  .controls {
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
    margin-bottom: 2rem;
    padding: 1rem;
    background: rgba(0, 0, 0, 0.03);
    border-radius: 8px;
  }

  .control-group {
    display: flex;
    align-items: center;
    gap: 0.5rem;
  }

  select,
  input {
    padding: 0.5rem;
    border-radius: 4px;
    border: 1px solid var(--surface-stroke);
    background: white;
  }

  .icon-category {
    margin-bottom: 3rem;
  }

  h2 {
    margin: 0 0 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--surface-stroke);
    font-size: 1.2rem;
  }

  .icon-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
    gap: 1rem;
  }

  .icon-card {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 1rem;
    border: 1px solid var(--surface-stroke);
    border-radius: 8px;
    background: white;
    transition: all 0.2s;
  }

  .icon-card:hover {
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    transform: translateY(-2px);
  }

  .icon-preview {
    height: 40px;
    width: 40px;
    display: flex;
    align-items: center;
    justify-content: center;
    margin-bottom: 0.75rem;
  }

  .icon-name {
    font-size: 0.875rem;
    font-weight: 500;
    margin-bottom: 0.5rem;
  }

  .icon-usage {
    font-size: 0.75rem;
    padding: 0.5rem;
    background: rgba(0, 0, 0, 0.03);
    border-radius: 4px;
    width: 100%;
    overflow-x: auto;
    white-space: nowrap;
    color: var(--text-secondary);
  }

  code {
    font-family: monospace;
  }

  .semantic-mapping {
    background: white;
    border: 1px solid var(--surface-stroke);
    padding: 1.5rem;
    border-radius: 8px;
    margin-bottom: 2rem;
  }

  pre {
    background: rgba(0, 0, 0, 0.03);
    padding: 1rem;
    border-radius: 4px;
    overflow-x: auto;
    font-family: monospace;
    margin: 1rem 0;
  }

  .mapping-examples {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
    gap: 1rem;
    margin: 1.5rem 0;
  }

  .mapping-example {
    display: flex;
    gap: 1rem;
    padding: 0.75rem;
    border: 1px solid var(--surface-stroke);
    border-radius: 6px;
  }

  .mapping-icon {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 32px;
    flex-shrink: 0;
  }

  .mapping-details {
    flex-grow: 1;
  }

  .mapping-type {
    font-weight: 700;
    font-size: 0.875rem;
    margin-bottom: 0.25rem;
  }

  .mapping-name {
    font-size: 0.75rem;
    color: var(--text-secondary);
    margin-bottom: 0.25rem;
  }

  .mapping-desc {
    font-size: 0.75rem;
    color: var(--text-tertiary);
  }
</style>
