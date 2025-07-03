#!/bin/zsh
brew install openjdk@17 
echo 'export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"' >> ~/.zshrc  
source ~/.zshrc 
java -version

# Install Nextflow
export CAPSULE_LOG=none
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
chmod +x /usr/local/bin/nextflow

# Check info
nextflow -version
nextflow info   