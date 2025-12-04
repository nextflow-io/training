
ARG CONNECT_CLIENT_VERSION="0.8-rc"

FROM cr.seqera.io/public/data-studio-vscode:1.101.2-${CONNECT_CLIENT_VERSION}

# Environment variable for repository URL with default
ENV REPO_URL="https://github.com/nextflow-io/training.git"
ENV REPO_REF="2.3.0"

# Create and configure workspace
WORKDIR /workspaces
RUN mkdir -p /workspaces/.nextflow && \
    mkdir -p /workspaces/.vscode-server

# Pre-install nf-core extension pack
RUN /home/.vscode/bin/openvscode-server --install-extension nf-core.nf-core-extensionpack --force

# Create setup script that clones repository then runs init
RUN cat > /usr/local/bin/setup.sh <<'BASH' && \
    chmod +x /usr/local/bin/setup.sh
#!/bin/bash
set -euo pipefail

echo "Cloning repository: ${REPO_URL}"

# Derive workspace dir from repo name
REPO_NAME=$(basename "${REPO_URL}" .git)
WORKSPACE_DIR="/workspaces/${REPO_NAME}"

# Check if repository already exists
if [ -d "${WORKSPACE_DIR}/.git" ]; then
    echo "Repository already exists at ${WORKSPACE_DIR}, skipping setup"
else
    # Clone the specified ref (branch or tag)
    git clone --depth 1 --branch "${REPO_REF}" "${REPO_URL}" "${WORKSPACE_DIR}"

    # VS Code settings scoped to the repo
    mkdir -p "${WORKSPACE_DIR}/.vscode"
    cat > "${WORKSPACE_DIR}/.vscode/settings.json" << 'EOF'
{
    "workbench.colorTheme": "Default Dark Modern",
    "terminal.integrated.defaultProfile.linux": "bash",
    "terminal.integrated.cwd": "${workspaceFolder}",
    "workbench.startupEditor": "none",
    "workbench.panel.defaultLocation": "bottom",
    "workbench.panel.opensMaximized": "never",
    "terminal.integrated.hideOnStartup": "never",
    "cSpell.diagnosticLevel": "Hint",
    "nextflow.java.home": "/usr",
    "nextflow.debug": false
}
EOF

    # Default folder points at the cloned repo
    sed -i "s|--default-folder '/workspace'|--default-folder '${WORKSPACE_DIR}'|g" /init

    # Pre-warm Nextflow
    export NXF_OPTS="-XX:+UseContainerSupport -XX:MaxRAMPercentage=50.0"
    nextflow help > /dev/null 2>&1
    sleep 2
    nextflow info > /dev/null 2>&1
    sleep 2
fi

exec /init
BASH

# Override the entrypoint to use connect-client
ENTRYPOINT ["/usr/bin/connect-client", "--entrypoint", "/usr/local/bin/setup.sh"]

# Default arguments (empty since setup script handles everything)
CMD []
