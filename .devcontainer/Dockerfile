FROM mcr.microsoft.com/devcontainers/javascript-node:1-20-bookworm@sha256:428d4ed6264774619bc353fc7dc98e41b4c01b603a638db9a549e66ed86edd15

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
