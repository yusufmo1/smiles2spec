# Docker Deployment for SMILES to Spectrum Application

This directory contains the Docker configuration for deploying both the frontend and backend of the SMILES to Spectrum application.

## Structure

- `Dockerfile.frontend` - Docker configuration for the Svelte frontend
- `Dockerfile.backend` - Docker configuration for the Flask backend
- `docker-compose.yml` - Orchestration file for both services
- `nginx.conf` - Nginx configuration for serving the frontend and proxying API requests
- `deploy.sh` - Deployment script for easy setup

**Note:** For convenience, symlinks to `Dockerfile.frontend`, `Dockerfile.backend`, and `docker-compose.yml` are created in the project root directory, allowing you to run Docker commands from there.

## Requirements

- Docker and Docker Compose installed on your mini PC
- Sufficient disk space for the application and its dependencies

## Deployment Steps

### Option 1: Using the Deployment Script (Recommended)

The easiest way to deploy is using the provided script:

```bash
cd /path/to/smiles2spec
./docker/deploy.sh
```

The script will:
1. Check if Docker and Docker Compose are installed
2. Build and start the containers
3. Check container status and display access information

### Option 2: Manual Deployment

1. Make sure Docker and Docker Compose are installed on your system:

```bash
docker --version
docker-compose --version
```

2. Clone this repository to your mini PC (if not already done).

3. Navigate to the project root directory:

```bash
cd /path/to/smiles2spec
```

4. Build and start the containers:

```bash
docker-compose up -d
```

5. To check if the containers are running properly:

```bash
docker-compose ps
```

6. Access the application:
   - Web interface: http://localhost or http://[your-mini-pc-ip]
   - API directly: http://localhost:5000/api/health or http://[your-mini-pc-ip]:5000/api/health

7. To stop the services:

```bash
docker-compose down
```

## Troubleshooting

### Logs

To view logs from the containers:

```bash
# View logs from both services
docker-compose logs

# View logs from a specific service
docker-compose logs frontend
docker-compose logs backend

# Follow logs (continuous display)
docker-compose logs -f
```

### Common Issues

1. **Port conflicts**: If ports 80 or 5000 are already in use on your mini PC, modify the port mappings in `docker-compose.yml`.

2. **Memory issues**: The backend requires adequate memory for loading ML models. If your mini PC has limited resources, you might need to adjust Docker's resource limits.

3. **Model files**: Ensure that the model files are available in the `backend/models` directory.

## Persistent Data

Models and other data are stored in volumes that persist between container restarts. The configuration mounts the local `backend/models` directory into the container. 