#!/usr/bin/env bash

# Install the Seqera Platform CLI (tw).
# Arch-aware so the multi-architecture image build works on both amd64 and arm64.
# Only downloads the binary (no execution), so it is safe under QEMU emulation.
set -eu

VERSION="${VERSION:-0.32.0}"

case "$(uname -m)" in
    x86_64 | amd64) ARCH="x86_64" ;;
    aarch64 | arm64) ARCH="arm64" ;;
    *)
        echo "Unsupported architecture: $(uname -m)" >&2
        exit 1
        ;;
esac

if [ "${VERSION}" = "latest" ]; then
    URL="https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-${ARCH}"
else
    URL="https://github.com/seqeralabs/tower-cli/releases/download/v${VERSION}/tw-linux-${ARCH}"
fi

echo "Installing Seqera CLI (tw) ${VERSION} for ${ARCH}..."
curl -fSL "${URL}" -o /usr/local/bin/tw
chmod +x /usr/local/bin/tw
