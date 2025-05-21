<script>
</script>

<div class="bg"></div>

<style>
  .bg {
    position: fixed;
    inset: 0;
    z-index: -1;
    background: #f7f7f8; 
    overflow: hidden;
    background: 
      radial-gradient(circle at 30% 30%, rgba(114, 135, 253, 0.8), transparent 70%),
      radial-gradient(circle at 70% 70%, rgba(146, 113, 255, 0.8), transparent 70%);
    filter: blur(60px);
    animation: shift-bg 20s ease-in-out infinite alternate;
  }
  
  .bg::before {
    content: '';
    position: absolute;
    width: 180%;
    height: 180%;
    top: -40%;
    left: -40%;
    background-blend-mode: overlay;
    animation: breathe 15s ease-in-out infinite alternate;
  }
  
  @keyframes breathe {
    from { transform: scale(1) rotate(0deg); opacity: 0.8; }
    to   { transform: scale(1.2) rotate(5deg); opacity: 1; }
  }
  
  @keyframes shift-bg {
    from { background-position: 0% 0%; }
    to   { background-position: 3% 3%; }
  }
</style>