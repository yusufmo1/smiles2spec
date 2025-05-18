#!/bin/bash

# Simple deployment script for the SMILES to Spectrum application

# Make script exit on first error
set -e

# Display banner
echo "================================================"
echo "SMILES to Spectrum Docker Deployment"
echo "================================================"

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Docker is not installed. Please install Docker first."
    exit 1
fi

# Check if Docker Compose is installed
if ! command -v docker-compose &> /dev/null; then
    echo "Docker Compose is not installed. Please install Docker Compose first."
    exit 1
fi

# Change to the project root directory (assumes script is in docker/)
cd "$(dirname "$0")/.."

# Display info
echo "Building and deploying containers..."
echo "This may take a few minutes..."

# Build and start the containers
docker-compose build
docker-compose up -d

# Check if containers are running
echo "Checking container status..."
docker-compose ps

# Display access information
echo ""
echo "================================================"
echo "Deployment completed successfully!"
echo "Access the application at: http://localhost"
echo "API health check: http://localhost:5000/api/health"
echo "================================================"
echo ""
echo "To view logs: docker-compose logs -f"
echo "To stop: docker-compose down"
echo "================================================" 