FROM ubuntu:14.04

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    libncurses-dev \
    unzip \
    zlib1g-dev

RUN curl -L https://github.com/lh3/bwa/archive/v0.7.15.zip > /root/bwa.zip
RUN unzip /root/bwa.zip -d /root
RUN rm /root/bwa.zip
RUN make -C /root/bwa-0.7.15
RUN ln -s /root/bwa-0.7.15/bwa /usr/bin/bwa
RUN curl -L https://github.com/samtools/htslib/archive/1.3.1.zip > /root/htslib.zip
RUN unzip /root/htslib.zip -d /root
RUN rm /root/htslib.zip
RUN make -C /root/htslib-1.3.1
RUN curl -L https://github.com/samtools/samtools/archive/1.3.1.zip > /root/samtools.zip
RUN unzip /root/samtools.zip -d /root
RUN rm /root/samtools.zip
RUN make HTSDIR=/root/htslib-1.3.1 -C /root/samtools-1.3.1
RUN ln -s /root/samtools-1.3.1/samtools /usr/bin/samtools
RUN curl -L https://github.com/GregoryFaust/samblaster/archive/v.0.1.24.zip > /root/samblaster.zip
RUN unzip /root/samblaster.zip -d /root
RUN make -C /root/samblaster-v.0.1.24
RUN ln -s /root/samblaster-v.0.1.24/samblaster /usr/bin/samblaster
RUN git clone https://github.com/statgen/libStatGen /root/libStatGen
RUN make -C /root/libStatGen
RUN git clone https://github.com/statgen/bamUtil /root/bamUtil
RUN git --git-dir=/root/bamUtil/.git --work-tree=/root/bamUtil checkout c8eae40d7824769bf63390fe16be82cc146b8d6f #NonPrimaryDedup
RUN make -C /root/bamUtil
RUN cp /root/bamUtil/bin/bam /usr/bin/bam-non-primary-dedup
RUN git --git-dir=/root/bamUtil/.git --work-tree=/root/bamUtil checkout b6e4a7de6b7ce08d488f539ada4f1717cd4d12e4 #ExternalMemorySortManager
RUN make -C /root/bamUtil
RUN cp /root/bamUtil/bin/bam /usr/bin/bam-ext-mem-sort-manager
