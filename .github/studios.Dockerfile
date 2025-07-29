
ARG CONNECT_CLIENT_VERSION="0.8-rc"

FROM cr.seqera.io/public/data-studio-vscode:1.101.2-${CONNECT_CLIENT_VERSION}

ARG IMAGE_TAG=2.2.1
ENV TRAINING_TAG=${IMAGE_TAG}

# Create and configure workspace
WORKDIR /workspaces
RUN mkdir -p /workspaces/.nextflow && \
    mkdir -p /workspaces/.vscode-server && \
    mkdir -p /workspaces/training

# Configure VS Code settings for dark mode and auto-open terminal
RUN mkdir -p /workspaces/training/.vscode
COPY <<EOF /workspaces/training/.vscode/settings.json
{
  "workbench.colorTheme": "Default Dark Modern",
  "terminal.integrated.defaultProfile.linux": "bash",
  "terminal.integrated.cwd": "${workspaceFolder}",
  "workbench.startupEditor": "none",
  "workbench.panel.defaultLocation": "bottom",
  "workbench.panel.opensMaximized": "never",
  "terminal.integrated.showOnStartup": "always",
  "cSpell.diagnosticLevel": "Hint",
  "nextflow.java.home": "/usr",
  "nextflow.debug": false
}
EOF

# Pre-install nf-core extension pack
RUN /home/.vscode/bin/openvscode-server --install-extension nf-core.nf-core-extensionpack --force

# Patch init script to use training directory as default folder
RUN sed -i 's|--default-folder '"'"'/workspace'"'"'|--default-folder '"'"'/workspaces/training'"'"'|g' /init

# Create setup script that downloads training materials then runs init
RUN echo -e '#!/bin/bash\n\
echo "Downloading training materials..."\n\
wget -qO- https://github.com/nextflow-io/training/archive/refs/tags/${TRAINING_TAG}.tar.gz | tar xz --strip-components=1 -C /workspaces/training\n\
echo "Initializing Nextflow..."\n\
# Set conservative JVM flags for memory management only\n\
export NXF_OPTS="-XX:+UseContainerSupport -XX:MaxRAMPercentage=50.0"\n\
# Pre-warm Nextflow\n\
nextflow help > /dev/null 2>&1\n\
sleep 2\n\
nextflow info > /dev/null 2>&1\n\
sleep 2\n\
echo "Starting VS Code..."\n\
exec /init' > /usr/local/bin/setup.sh && \
    chmod +x /usr/local/bin/setup.sh

# Override the entrypoint to use connect-client
ENTRYPOINT ["/usr/bin/connect-client", "--entrypoint", "/usr/local/bin/setup.sh"]

# Default arguments (empty since setup script handles everything)
CMD []
