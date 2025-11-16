#!/usr/bin/env bash

set -euo pipefail

# Check if CSV file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

csv_file="$1"

# Check if file exists
if [ ! -f "$csv_file" ]; then
    echo "Error: File '$csv_file' not found"
    exit 1
fi

# Skip header if present and process each line
tail -n +2 "$csv_file" | while IFS=, read -r url attribution; do
    # Remove leading/trailing whitespace and quotes
    url=$(echo "$url" | xargs)
    attribution=$(echo "$attribution" | xargs | sed 's/^"//;s/"$//')

    # Extract photo ID from URL (the part after /photos/)
    photo_id=$(echo "$url" | grep -oP '(?<=photos/)[^/]+' || echo "unknown")

    # Default to jpg for Unsplash images
    ext="jpg"

    echo "Downloading: $photo_id"

    # Download the image with wget (quiet mode)
    if wget -q -O "${photo_id}.${ext}" "$url"; then
        echo "✓ Saved image: ${photo_id}.${ext}"
    else
        echo "✗ Failed to download: $url"
        continue
    fi

    # Save attribution to text file
    echo "$attribution" > "${photo_id}.txt"
    echo "✓ Saved attribution: ${photo_id}.txt"
    echo ""
done

echo "Download complete!"
