FROM ubuntu:14.04

LABEL maintainer="jshands@ucsc.edu"

RUN apt-get update && apt-get install -y \
   sudo \
   nano \
   cmake \
   git \
   wget \
   vim \
   python \
   build-essential \
   zlib1g-dev \
   libcurl4-gnutls-dev \
   libbz2-dev \
   libssl-dev \
   libncurses5-dev \
   liblzma-dev \
   unzip \
   autoconf

WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -O samtools.tar.bz2
RUN tar -xjvf samtools.tar.bz2
WORKDIR samtools-1.9
RUN ./configure --prefix /opt/samtools/samtools-1.9
RUN make
RUN make install
ENV PATH /opt/samtools/samtools-1.9/bin:$PATH

WORKDIR	/opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2 
RUN tar -xjvf htslib.tar.bz2 
WORKDIR	htslib-1.3.2
RUN ./configure --prefix /opt/htslib/htslib-1.3.2
RUN make
RUN make install
ENV PATH /opt/htslib/htslib-1.3.2/bin:$PATH

# wget the Git repository by commit since this will create a smaller image than 
# git clone because it does not include the .git dir
# Commit c3c2c4445cb164b63462cbdfcc3efa3b54c4aff2 corresponds to release 1.0.5
ARG COMMIT=c3c2c4445cb164b63462cbdfcc3efa3b54c4aff2
WORKDIR /opt/verifybamid
RUN wget -O VerifyBamID.zip https://github.com/Griffan/VerifyBamID/archive/$COMMIT.zip
RUN unzip VerifyBamID.zip
WORKDIR /opt/verifybamid/VerifyBamID
RUN mv /opt/verifybamid/VerifyBamID-$COMMIT/* /opt/verifybamid/VerifyBamID
WORKDIR /opt/verifybamid/VerifyBamID/build
RUN cmake .. && make && make test
ENV PATH /opt/verifybamid/VerifyBamID/bin:$PATH

COPY Dockerfile /opt/
