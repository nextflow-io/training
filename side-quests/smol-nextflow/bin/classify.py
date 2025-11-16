#!/usr/bin/env -S uv run
# /// script
# dependencies = [
#   "torch>=2.0.0",
#   "pillow>=10.0.0",
#   "open-clip-torch>=2.20.0",
# ]
# ///
"""Classify images using MetaCLIP with local weights."""

from pathlib import Path
from PIL import Image
import argparse
import json
import open_clip
import sys
import torch

def classify_images(image_dir: Path, labels: list[str], model_path: Path,
                   json_output: bool = False, architecture: str = None):
    """Classify images in a directory using MetaCLIP."""
    # Auto-detect architecture from model filename if not provided
    if architecture is None:
        filename = model_path.name.lower()
        if 'b32' in filename or '32' in filename:
            architecture = 'ViT-B-32-quickgelu'
        elif 'b16' in filename or '16' in filename:
            architecture = 'ViT-B-16-quickgelu'
        elif 'l14' in filename:
            architecture = 'ViT-L-14-quickgelu'
        elif 'h14' in filename:
            architecture = 'ViT-H-14-quickgelu'
        else:
            raise ValueError(f"Cannot infer architecture from {model_path.name}. "
                           f"Please specify --architecture")

    print(f"Using architecture: {architecture}", file=sys.stderr)

    # Load model and preprocessing
    model, _, preprocess = open_clip.create_model_and_transforms(
        architecture,
        pretrained=str(model_path),
        weights_only=False  # Trust MetaCLIP checkpoint from Facebook Research
    )
    tokenizer = open_clip.get_tokenizer(architecture)
    # Prepare text labels
    text = tokenizer(labels)

    # Process each image
    results = []
    for img_path in image_dir.glob('*'):
        if img_path.suffix.lower() in ['.jpg', '.jpeg', '.png', '.gif', '.bmp']:
            try:
                image = preprocess(Image.open(img_path)).unsqueeze(0)

                with torch.no_grad():
                    image_features = model.encode_image(image)
                    text_features = model.encode_text(text)

                    # Normalize and compute similarity
                    image_features /= image_features.norm(dim=-1, keepdim=True)
                    text_features /= text_features.norm(dim=-1, keepdim=True)
                    similarity = (100.0 * image_features @ text_features.T).softmax(dim=-1)

                # Get all confidences for all labels
                confidences = {label: float(conf) for label, conf in zip(labels, similarity[0])}

                # Get top prediction
                values, indices = similarity[0].topk(1)
                prediction = labels[indices[0]]
                confidence = values[0].item()

                result = {
                    'file': img_path.name,
                    'path': str(img_path),
                    'prediction': prediction,
                    'confidence': confidence,
                    'all_confidences': confidences
                }

                results.append(result)

                if json_output:
                    # Print JSONL (one JSON object per line)
                    print(json.dumps(result))
                else:
                    # Print human-readable format
                    print(f"{img_path.name}: {prediction} ({confidence:.2%})")

            except Exception as e:
                error_result = {
                    'file': img_path.name,
                    'path': str(img_path),
                    'error': str(e)
                }
                if json_output:
                    print(json.dumps(error_result))
                else:
                    print(f"Error processing {img_path.name}: {e}")

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify images using MetaCLIP")
    parser.add_argument(
        '--image-dir',
        type=Path,
        default=Path("data/pics"),
        help='Directory containing images to classify (default: data/pics)'
    )
    parser.add_argument(
        '--model-path',
        type=Path,
        default=Path("data/models/b32_400m.pt"),
        help='Path to MetaCLIP model weights (default: data/models/b32_400m.pt)'
    )
    parser.add_argument(
        '--labels',
        nargs='+',
        default=["cute cat", "ugly cat"],
        help='Labels for classification (default: ["cute cat", "ugly cat"])'
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output results as JSONL (one JSON object per line) to stdout'
    )
    parser.add_argument(
        '--architecture',
        type=str,
        choices=['ViT-B-32-quickgelu', 'ViT-B-16-quickgelu', 'ViT-L-14-quickgelu', 'ViT-H-14-quickgelu'],
        help='Model architecture (auto-detected from filename if not specified)'
    )

    args = parser.parse_args()

    results = classify_images(
        args.image_dir,
        args.labels,
        args.model_path,
        args.json,
        args.architecture
    )

