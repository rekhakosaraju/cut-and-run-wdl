FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    default-jre wget unzip python3 python3-pip git build-essential \
    samtools bedtools awscli && \
    rm -rf /var/lib/apt/lists/*

# STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    cp STAR-2.7.11b/bin/Linux_x86_64/STAR /usr/local/bin/

# Trim Galore + cutadapt + FastQC
RUN pip3 install cutadapt && \
    wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz && \
    tar -xzf 0.6.10.tar.gz && cp TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && chmod +x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/fastqc

# Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -O /opt/picard.jar

# SEACR
RUN git clone https://github.com/MicrosoftFreddyao/SEACR.git /opt/SEACR || \
    git clone https://github.com/FredHutch/SEACR.git /opt/SEACR

# UCSC tools (bedGraphToBigWig)
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -O /usr/local/bin/bedGraphToBigWig && \
    chmod +x /usr/local/bin/bedGraphToBigWig

WORKDIR /data