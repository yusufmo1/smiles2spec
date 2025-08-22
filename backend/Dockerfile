FROM python:3.10-slim

# System dependencies including X11 libs needed by RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    libxrender1 \
    libxext6 \
    libfontconfig1 \
    libcairo2 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf2.0-0 \
    libglib2.0-0 \
    libx11-6 \
    libx11-xcb1 \
    libxcb1 \
    libxcomposite1 \
    libxcursor1 \
    libxdamage1 \
    libxfixes3 \
    libxi6 \
    libxrandr2 \
    libxss1 \
    libxtst6 \
    libnss3 \
    && rm -rf /var/lib/apt/lists/*

# Working directory
WORKDIR /app

# Upgrade pip first to avoid version conflicts
RUN pip install --upgrade pip==25.1.1

# Copy requirements first for better caching
COPY backend/requirements.txt .

# Install all dependencies from requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend code
COPY backend /app/backend

# Environment variables
ENV FLASK_ENV=production \
    PYTHONUNBUFFERED=1 \
    PORT=5050

EXPOSE 5050

# Add health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:5050/health || exit 1

# Use Gunicorn for production WSGI server with module import
CMD ["gunicorn", "--config", "backend/gunicorn.conf.py", "backend.wsgi:app"] 