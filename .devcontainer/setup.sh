#!/usr/bin/env bash

# Fix for Java options
printf 'unset JAVA_TOOL_OPTIONS\n' >> $HOME/.bashrc
unset JAVA_TOOL_OPTIONS

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

# Force Java environment variables to use conda installation
printf 'export JAVA_HOME=/opt/conda\nexport JAVA_CMD=/opt/conda/bin/java\n' >> $HOME/.bashrc
export JAVA_HOME=/opt/conda
export JAVA_CMD=/opt/conda/bin/java

# Update Nextflow
nextflow self-update
nextflow -version

cat /usr/local/etc/vscode-dev-containers/first-run-notice.txt
