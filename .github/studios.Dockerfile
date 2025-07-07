# Test Dockerfile for Studios-compatible Nextflow Training Container
# This uses the existing pre-built training container and adds Studios components

# Add a default Connect client version
ARG CONNECT_CLIENT_VERSION="0.8"

# Get the Seqera connect client
FROM public.cr.seqera.io/platform/connect-client:${CONNECT_CLIENT_VERSION} AS connect

# Use the existing pre-built training container as base
FROM ghcr.io/nextflow-io/training:latest

# Set default training tag version to match container version
ARG IMAGE_TAG=latest
ENV TRAINING_TAG=${IMAGE_TAG}

# Install code-server and wget
RUN apt-get update && apt-get install -y wget && \
    curl -fsSL https://code-server.dev/install.sh | sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the connect binary from the connect client image
COPY --from=connect /usr/bin/connect-client /usr/bin/connect-client

# Install connect dependencies
RUN /usr/bin/connect-client --install && \
    /usr/bin/connect-client --setToolIdentifier vscode

# Set up workspace environment
ENV HOME=/workspace \
    EDITOR=code \
    VISUAL=code \
    GIT_EDITOR="code --wait" \
    SHELL=bash

# Create and configure workspace
WORKDIR /workspace
RUN mkdir -p /workspace/.nextflow && \
    mkdir -p /workspace/.vscode-server && \
    mkdir -p /workspace/training

# Create setup script
RUN echo '#!/bin/bash\n\
rm -rf /workspace/training/*\n\
wget -qO- https://github.com/nextflow-io/training/archive/refs/tags/${TRAINING_TAG}.tar.gz | tar xz --strip-components=1 -C /workspace/training\n\
nextflow help > /dev/null\n\
exec "$@"' > /usr/local/bin/setup.sh && \
    chmod +x /usr/local/bin/setup.sh

# Pre-heat nextflow
RUN nextflow help > /dev/null

# Override the entrypoint to use connect-client
ENTRYPOINT ["/usr/bin/connect-client", "--entrypoint", "/usr/local/bin/setup.sh"]

# Start code-server on the port specified by Studios, in the workspace directory
CMD ["/usr/bin/bash", "-c", "code-server --host 0.0.0.0 --port ${CONNECT_TOOL_PORT} --auth none --user-data-dir /workspace/.vscode-server /workspace/training"]
