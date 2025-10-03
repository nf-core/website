FROM mcr.microsoft.com/devcontainers/javascript-node:1-20-bookworm

# Copy package files
COPY package.json package-lock.json* ./
COPY sites/main-site/package.json ./sites/main-site/
COPY sites/docs/package.json ./sites/docs/
COPY sites/configs/package.json ./sites/configs/
COPY sites/pipelines/package.json ./sites/pipelines/
COPY sites/modules-subworkflows/package.json ./sites/modules-subworkflows/
COPY sites/pipeline-results/package.json ./sites/pipeline-results/

# Install dependencies
RUN npm install --workspaces
