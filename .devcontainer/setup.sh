#!/usr/bin/env bash

# Fix for Java options
printf 'unset JAVA_TOOL_OPTIONS\n' >> $HOME/.bashrc
unset JAVA_TOOL_OPTIONS

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

# Update Nextflow
nextflow self-update
nextflow -version

# Debug message showing where we're running
if [ -z \"$CODESPACES\" ]
then
    echo " 🖥️ Running on: Devcontainers Development"
else
    echo " 🌐 Running on: Codespaces Development"
fi
