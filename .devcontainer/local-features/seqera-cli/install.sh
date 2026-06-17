#!/usr/bin/env bash

# Install the Seqera Platform CLI (tw).
# Arch-aware so the multi-architecture image build works on both amd64 and arm64.
# tower-cli only publishes a native Linux binary for x86_64; on arm64 we install
# the architecture-independent JAR (Java is provided by the java feature).
# Only downloads files (no execution), so it is safe under QEMU emulation.
set -eu

VERSION="${VERSION:-0.32.0}"

if [ "${VERSION}" = "latest" ]; then
    BASE="https://github.com/seqeralabs/tower-cli/releases/latest/download"
else
    BASE="https://github.com/seqeralabs/tower-cli/releases/download/v${VERSION}"
fi

case "$(uname -m)" in
    x86_64 | amd64)
        echo "Installing Seqera CLI (tw) ${VERSION} native binary for x86_64..."
        curl -fSL "${BASE}/tw-linux-x86_64" -o /usr/local/bin/tw
        chmod +x /usr/local/bin/tw
        ;;
    aarch64 | arm64)
        # No native Linux arm64 binary is published; use the JAR via Java.
        echo "Installing Seqera CLI (tw) ${VERSION} JAR for arm64..."
        mkdir -p /usr/local/lib/tw
        curl -fSL "${BASE}/tw-jar.jar" -o /usr/local/lib/tw/tw.jar
        cat > /usr/local/bin/tw <<'EOF'
#!/usr/bin/env bash
exec java -jar /usr/local/lib/tw/tw.jar "$@"
EOF
        chmod +x /usr/local/bin/tw
        ;;
    *)
        echo "Unsupported architecture: $(uname -m)" >&2
        exit 1
        ;;
esac
