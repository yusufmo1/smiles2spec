FROM python:3.10-slim

# System dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    libxrender1 \
    libxext6 \
    libfontconfig1 \
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