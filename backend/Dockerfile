FROM python:3.10-slim

# System dependencies including X11 libs needed by RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
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

# Copy requirements first for better caching
COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt && \
    pip install flask flask-cors

# Copy backend code
COPY backend /app/backend

# Environment variables
ENV FLASK_ENV=production \
    PYTHONUNBUFFERED=1 \
    PORT=5050

EXPOSE 5050
CMD ["python", "-m", "backend.app"] 