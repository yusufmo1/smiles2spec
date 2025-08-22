"""
Gunicorn Configuration for SMILES2SPEC Backend

Optimized for ML workload with RDKit and scikit-learn models.
Uses sync workers to ensure compatibility with ML libraries.
"""

import os
import multiprocessing

# Server socket
bind = f"0.0.0.0:{os.getenv('PORT', '5050')}"
backlog = 2048

# Worker processes
# Limited to 2 workers for ML model memory efficiency
workers = int(os.getenv('GUNICORN_WORKERS', min(2, multiprocessing.cpu_count())))
worker_class = "sync"  # Sync workers for ML model compatibility
worker_connections = 1000
max_requests = int(os.getenv('GUNICORN_MAX_REQUESTS', '1000'))
max_requests_jitter = int(os.getenv('GUNICORN_MAX_REQUESTS_JITTER', '100'))

# Timeout settings - extended for ML inference
timeout = int(os.getenv('GUNICORN_TIMEOUT', '300'))  # 5 minutes for ML operations
keepalive = int(os.getenv('GUNICORN_KEEPALIVE', '5'))
graceful_timeout = int(os.getenv('GUNICORN_GRACEFUL_TIMEOUT', '30'))

# Application
preload_app = True  # Critical for ML model loading efficiency
sendfile = False

# Logging
loglevel = os.getenv('GUNICORN_LOG_LEVEL', 'info').lower()
accesslog = '-'  # Log to stdout
errorlog = '-'   # Log to stderr
access_log_format = '%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s" %(D)s'

# Process naming
proc_name = 'smiles2spec-backend'

# Security
limit_request_line = 4094
limit_request_fields = 100
limit_request_field_size = 8190

# Performance tuning
worker_tmp_dir = '/dev/shm'  # Use shared memory for better performance

def on_starting(server):
    """Called just before the master process is initialized."""
    server.log.info("Starting SMILES2SPEC Backend with Gunicorn")

def when_ready(server):
    """Called just after the server is started."""
    server.log.info(f"Server is ready. Listening on {bind}")

def on_reload(server):
    """Called to recycle workers during a reload via SIGHUP."""
    server.log.info("Reloading workers")

def worker_int(worker):
    """Called just after a worker exited on SIGINT or SIGQUIT."""
    worker.log.info(f"Worker {worker.pid} received INT or QUIT signal")

def pre_fork(server, worker):
    """Called just before a worker is forked."""
    server.log.info(f"Worker {worker.pid} spawned")

def post_fork(server, worker):
    """Called just after a worker has been forked."""
    server.log.info(f"Worker {worker.pid} spawned successfully")

def worker_abort(worker):
    """Called when a worker receives the SIGABRT signal."""
    worker.log.info(f"Worker {worker.pid} received SIGABRT signal")