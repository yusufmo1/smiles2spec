<script>
  export let png = null;
</script>

<div class="structure-wrapper">
  {#if png}
    <img src={`data:image/png;base64,${png}`} alt="Molecular structure" />
  {:else}
    <div class="placeholder">
      <span class="placeholder-icon">🧪</span>
      <p>Awaiting molecule</p>
    </div>
  {/if}
</div>

<style>
.structure-wrapper{
  display: flex;
  align-items: center;
  justify-content: center;
  width: 100%;
  height: 100%;
  overflow: hidden;
}
img {
  max-width: 100%;
  max-height: 100%;
  object-fit: contain;
}
.placeholder{
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100%;
  width: 100%;
  background: rgba(0, 0, 0, 0.05);
  border-radius: 8px;
  color: rgba(0, 0, 0, 0.4);
}
.placeholder-icon {
  font-size: 3rem;
  margin-bottom: 1rem;
}
p {
  font-size: 0.9rem;
  font-weight: 500;
}
</style> 