#!/usr/bin/env bash

set -euo pipefail

# Show help message
show_help() {
    cat << EOF
Usage: $0 [-n|--count NUM] [-p|--prefix PATH] [-h|--help]

Download random cat images from cataas.com along with their tags.

Options:
    -c, --count NUM     Number of cats to download (default: 1)
    -p, --prefix PATH   Directory to download files to (default: current directory)
    -h, --help          Show this help message and exit

Examples:
    $0                      # Download 1 cat to current directory
    $0 -n 5                 # Download 5 cats to current directory
    $0 -p cats -n 10        # Download 10 cats to ./cats directory
    $0 --prefix data/cats   # Download 1 cat to ./data/cats directory

Output:
    For each cat, creates two files:
    - <ID>.jpg          The cat image
    - <ID>.txt          The tags (one per line)

EOF
}

# Default values
num_downloads=1
prefix="."

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -c|--count)
            num_downloads="$2"
            shift 2
            ;;
        -p|--prefix)
            prefix="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate number
if ! [[ "$num_downloads" =~ ^[0-9]+$ ]] || [ "$num_downloads" -lt 1 ]; then
    echo "Error: Number must be a positive integer"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$prefix"
echo "Downloading $num_downloads cat(s) to $prefix/"
echo ""

# Download loop
for i in $(seq 1 "$num_downloads"); do
    echo "[$i/$num_downloads]"

    # Get the JSON metadata
    json=$(curl -s 'https://cataas.com/cat?json=true')

    # Extract fields using jq
    cat_id=$(echo "$json" | jq -r '.id')
    mimetype=$(echo "$json" | jq -r '.mimetype')
    url=$(echo "$json" | jq -r '.url')
    tags=$(echo "$json" | jq -r '.tags[]')  # Extract tags, one per line

    # Map mimetype to extension - only accept jpg, png, gif
    case "$mimetype" in
        image/jpeg|image/jpg)
            ext="jpg"
            ;;
        image/png)
            ext="png"
            ;;
        image/gif)
            ext="gif"
            ;;
        *)
            echo "✗ Skipping unsupported type: $mimetype"
            echo ""
            continue
            ;;
    esac

    # Build filenames with prefix
    filename="${prefix}/${cat_id}.${ext}"
    tagfile="${prefix}/${cat_id}.txt"

    # Download the image
    curl -s "$url" -o "$filename"
    echo "✓ Saved as $filename"

    # Save tags to text file
    echo "$tags" > "$tagfile"
    echo "✓ Saved tags to $tagfile"
    echo ""
done

echo "Download complete! Downloaded $num_downloads cat(s) to $prefix/"
