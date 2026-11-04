#!/usr/bin/env bash

# Initial commands to run, before everything else

# Set up directories
mkdir -p /workspaces/.nextflow
mkdir -p /workspaces/training/

# Copy over the welcome message
cp welcome-message.txt /usr/local/etc/vscode-dev-containers/first-run-notice.txt
