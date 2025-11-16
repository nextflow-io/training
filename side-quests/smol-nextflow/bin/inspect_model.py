#!/usr/bin/env -S uv run
# /// script
# dependencies = [
#   "torch>=2.0.0",
#   "numpy>=1.24.0",
# ]
# ///
"""Inspect MetaCLIP model checkpoint metadata."""

import argparse
from pathlib import Path
import torch

def inspect_checkpoint(checkpoint_path: Path):
    """Load and inspect a PyTorch checkpoint file."""

    print(f"Loading checkpoint: {checkpoint_path}")
    print(f"File size: {checkpoint_path.stat().st_size / (1024**2):.2f} MB\n")

    # Load the checkpoint
    checkpoint = torch.load(checkpoint_path, map_location='cpu', weights_only=False)

    # Check what type of object it is
    print(f"Checkpoint type: {type(checkpoint)}")
    print()

    # If it's a dict, explore its keys
    if isinstance(checkpoint, dict):
        print("Top-level keys:")
        for key in checkpoint.keys():
            print(f"  - {key}")
        print()

        # Look for common metadata keys
        metadata_keys = ['epoch', 'metadata', 'config', 'args', 'model_config',
                        'train_config', 'version', 'architecture']

        print("Metadata found:")
        for key in metadata_keys:
            if key in checkpoint:
                print(f"\n{key}:")
                print(f"  {checkpoint[key]}")

        # Check state_dict structure
        if 'state_dict' in checkpoint:
            state_dict = checkpoint['state_dict']
            print(f"\nstate_dict contains {len(state_dict)} tensors")
            print("First 10 tensor names:")
            for i, name in enumerate(list(state_dict.keys())[:10]):
                shape = state_dict[name].shape
                print(f"  {name}: {shape}")

        # If checkpoint IS the state_dict directly
        elif all(isinstance(v, torch.Tensor) for v in list(checkpoint.values())[:5]):
            print(f"\nCheckpoint appears to be a state_dict with {len(checkpoint)} tensors")
            print("First 10 tensor names:")
            for i, (name, tensor) in enumerate(list(checkpoint.items())[:10]):
                print(f"  {name}: {tensor.shape}")

    else:
        print(f"Checkpoint is not a dict, it's: {checkpoint}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect PyTorch checkpoint metadata")
    parser.add_argument(
        'checkpoint',
        type=Path,
        help='Path to .pt checkpoint file'
    )

    args = parser.parse_args()

    if not args.checkpoint.exists():
        print(f"Error: {args.checkpoint} does not exist")
        exit(1)

    inspect_checkpoint(args.checkpoint)
