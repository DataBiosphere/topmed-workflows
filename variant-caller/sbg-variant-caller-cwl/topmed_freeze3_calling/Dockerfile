FROM ubuntu:16.04

ARG COMMIT=5f70239e390552d3a75b4e9fc098c342131958a6

RUN apt-get update && apt-get install -y \
    git \
    python-pip \
    python \
    python3.5 \
    wget \
    curl \
    unzip \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    htop \
    nano

WORKDIR root

# wget the Git repository by commit since this will create a smaller image than 
# git clone because it does not include the .git dir
RUN wget -O top_med_freeze3_calling.zip https://github.com/statgen/topmed_freeze3_calling/archive/$COMMIT.zip
RUN unzip top_med_freeze3_calling.zip
RUN mkdir topmed_freeze3_calling
RUN mv topmed_freeze3_calling-$COMMIT/* topmed_freeze3_calling/
RUN rm -rf topmed_freeze3_calling-$COMMIT

WORKDIR topmed_freeze3_calling

# Build the tools included in the image
RUN make SHELL=/bin/bash

COPY Dockerfile /opt/