FROM debian:stretch-slim

MAINTAINER <your name>

RUN apt-get update && apt-get install -y curl cowsay

ENV PATH=$PATH:/usr/games/