#!/usr/bin/env bash

# Install python cli tools using uv

uv tool install pre-commit
uv tool install nf-core==3.5.2
uv tool install "mkdocs-quiz>=1.5.2"
# seqerakit: used by the platform_automation side quest
uv tool install "seqerakit==0.5.7"
