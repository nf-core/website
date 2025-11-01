FROM mcr.microsoft.com/devcontainers/javascript-node:4-20-bookworm@sha256:b32c7d3f03c813ff653ca5859aa74946a7c98663a2aa2c0affe1265e17e5a2e7

# Set working directory
WORKDIR /workspace

# Copy root package files
COPY package.json package-lock.json* ./

# Copy only package.json files from each workspace
COPY sites/configs/package.json ./sites/configs/package.json
COPY sites/docs/package.json ./sites/docs/package.json
COPY sites/main-site/package.json ./sites/main-site/package.json
COPY sites/modules-subworkflows/package.json ./sites/modules-subworkflows/package.json
COPY sites/pipeline-results/package.json ./sites/pipeline-results/package.json
COPY sites/pipelines/package.json ./sites/pipelines/package.json

# Install dependencies
RUN npm install --workspaces
