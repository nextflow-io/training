#!/usr/bin/env bash

# Install the tooling used by the platform_automation side quest:
#   - terraform : provision compute environments and pipelines declaratively
#   - tw        : the Seqera Platform CLI (tower-cli)
#   - seqerakit : YAML-driven wrapper around tw for declarative Platform setup
#   - az        : the Azure CLI, used by Terraform's azurerm provider (az login)
#
# These run at container startup (via setup.sh -> onCreateCommand) rather than
# as devcontainer features. The production image (ghcr.io/nextflow-io/training)
# is only rebuilt on pushes to master, so feature additions don't reach an
# already-published image. Installing here makes the tools available in every
# container immediately.
#
# Idempotent (skips tools that are already present) and architecture-aware
# (amd64 + arm64) so re-runs and both CPU families are safe.
set -euo pipefail

TERRAFORM_VERSION="1.9.8"
TW_VERSION="0.32.0"

ARCH="$(uname -m)"

# Make sure the download/extract tools the installs below rely on exist.
ensure_pkg() {
    local cmd="$1" pkg="$2"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Installing ${pkg}..."
        apt-get update -qq && apt-get install -y -qq "${pkg}" >/dev/null
    fi
}

ensure_pkg curl curl
ensure_pkg unzip unzip

# --- Terraform -------------------------------------------------------------
if command -v terraform >/dev/null 2>&1; then
    echo "terraform already installed: $(terraform version | head -n1)"
else
    case "${ARCH}" in
        x86_64 | amd64) tf_arch="amd64" ;;
        aarch64 | arm64) tf_arch="arm64" ;;
        *) echo "Unsupported architecture for terraform: ${ARCH}" >&2; exit 1 ;;
    esac
    echo "Installing Terraform ${TERRAFORM_VERSION} (${tf_arch})..."
    tmp="$(mktemp -d)"
    curl -fSL "https://releases.hashicorp.com/terraform/${TERRAFORM_VERSION}/terraform_${TERRAFORM_VERSION}_linux_${tf_arch}.zip" -o "${tmp}/terraform.zip"
    unzip -q "${tmp}/terraform.zip" -d "${tmp}"
    install -m 0755 "${tmp}/terraform" /usr/local/bin/terraform
    rm -rf "${tmp}"
fi

# --- Seqera CLI (tw) -------------------------------------------------------
# tower-cli only publishes a native Linux binary for x86_64; on arm64 install
# the architecture-independent JAR (Java is provided by the java feature).
if command -v tw >/dev/null 2>&1; then
    echo "tw already installed: $(tw --version 2>/dev/null | head -n1)"
else
    base="https://github.com/seqeralabs/tower-cli/releases/download/v${TW_VERSION}"
    case "${ARCH}" in
        x86_64 | amd64)
            echo "Installing Seqera CLI (tw) ${TW_VERSION} native binary for x86_64..."
            curl -fSL "${base}/tw-linux-x86_64" -o /usr/local/bin/tw
            chmod +x /usr/local/bin/tw
            ;;
        aarch64 | arm64)
            echo "Installing Seqera CLI (tw) ${TW_VERSION} JAR for arm64..."
            mkdir -p /usr/local/lib/tw
            curl -fSL "${base}/tw-jar.jar" -o /usr/local/lib/tw/tw.jar
            cat > /usr/local/bin/tw <<'EOF'
#!/usr/bin/env bash
exec java -jar /usr/local/lib/tw/tw.jar "$@"
EOF
            chmod +x /usr/local/bin/tw
            ;;
        *) echo "Unsupported architecture for tw: ${ARCH}" >&2; exit 1 ;;
    esac
fi

# --- seqerakit -------------------------------------------------------------
# Pure-Python (arch-independent) wrapper that drives tw from YAML. Installed
# with the conda pip (python.defaultInterpreterPath = /opt/conda/bin/python),
# so no PEP 668 override is needed. Depends on tw being on PATH (installed above).
if command -v seqerakit >/dev/null 2>&1; then
    echo "seqerakit already installed: $(seqerakit --version 2>/dev/null | head -n1)"
else
    echo "Installing seqerakit..."
    pip install --quiet seqerakit
fi

# --- Azure CLI (az) --------------------------------------------------------
# Microsoft's deb installer sets up the apt repo and pulls the right package
# for the host architecture (amd64 + arm64 are both published).
if command -v az >/dev/null 2>&1; then
    echo "az already installed: $(az version --output tsv --query '"azure-cli"' 2>/dev/null | head -n1)"
else
    echo "Installing Azure CLI (az)..."
    curl -sL https://aka.ms/InstallAzureCLIDeb | bash
fi
