#!/usr/bin/env bash

# Install Seqera AI CLI (BETA, BUSL-1.1).
# Tracks the @dev npm channel so the training image picks up fixes quickly.
# Requires node + npm to be available before this feature runs.

set -euo pipefail

if ! command -v npm >/dev/null 2>&1; then
    echo "seqera-ai feature: npm not found - skipping install." >&2
    exit 0
fi

npm install -g seqera@dev
