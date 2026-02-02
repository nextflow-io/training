#!/bin/bash
# Build multilingual MkDocs site for deployment with mike
# Usage: ./build-multilang.sh <version> [output_dir]

set -e

VERSION="${1:-latest}"
OUTPUT_DIR="${2:-/tmp/multilang-build}"
DOCS_DIR="$(cd "$(dirname "$0")/../docs" && pwd)"

# Dynamically discover languages from docs/ directories containing mkdocs.yml
LANGUAGES=($(find "$DOCS_DIR" -maxdepth 2 -name "mkdocs.yml" -exec dirname {} \; | xargs -n1 basename | sort))

echo "Building multilingual site for version: $VERSION"
echo "Output directory: $OUTPUT_DIR"
echo "Languages: ${LANGUAGES[*]}"

# Clean output directory
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Build each language
for lang in "${LANGUAGES[@]}"; do
    echo ""
    echo "=== Building $lang ==="

    CONFIG_FILE="$DOCS_DIR/$lang/mkdocs.yml"

    if [ ! -f "$CONFIG_FILE" ]; then
        echo "Warning: Config not found for $lang at $CONFIG_FILE, skipping..."
        continue
    fi

    # Determine the output path based on language
    if [ "$lang" = "en" ]; then
        LANG_OUTPUT="$OUTPUT_DIR"
    else
        LANG_OUTPUT="$OUTPUT_DIR/$lang"
    fi

    # Build with overridden site_url
    # For EN: site_url is base (e.g., https://training.nextflow.io/latest/)
    # For other langs: site_url includes lang (e.g., https://training.nextflow.io/latest/pt/)
    if [ "$lang" = "en" ]; then
        SITE_URL="https://training.nextflow.io/$VERSION/"
    else
        SITE_URL="https://training.nextflow.io/$VERSION/$lang/"
    fi

    echo "Building to: $LANG_OUTPUT"
    echo "Site URL: $SITE_URL"

    # Run mkdocs build with proper site_url
    cd "$DOCS_DIR/$lang"
    mkdocs build -f "$CONFIG_FILE" -d "$LANG_OUTPUT" 2>&1 | tail -5

    echo "Built $lang successfully"
done

echo ""
echo "=== Build complete ==="
echo "Output: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"
