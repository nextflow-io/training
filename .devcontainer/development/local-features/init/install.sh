#!/usr/bin/env bash

# Initial commands to run, before everything else

# Set up directories
mkdir -p /workspaces/.nextflow
mkdir -p /workspaces/training/

# Fix for Java options
printf 'unset JAVA_TOOL_OPTIONS\n' >> $HOME/.bashrc

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
