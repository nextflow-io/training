#!/usr/bin/env bash

# Install Seqera Platform "Tower Agent"
curl -fSL https://github.com/seqeralabs/tower-agent/releases/latest/download/tw-agent-linux-x86_64 > tw-agent
chmod +x tw-agent
mv tw-agent /usr/local/bin/tw-agent
