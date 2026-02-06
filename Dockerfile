# Use an official Python runtime as base image, slim-bookworm carries java 17
FROM python:3.10-slim-bookworm

# Set working directory
WORKDIR /app

# Install system dependencies required for Hail and other packages
RUN apt-get update && apt-get install -y \
    openjdk-17-jdk-headless \
    gcc g++ \
    libgomp1 \
    git \
    build-essential \
    autoconf automake libtool pkg-config \
    zlib1g-dev \
    bzip2 libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*


# Set JAVA_HOME environment variable (required for Hail)
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64

# Install Python dependencies
RUN pip install --no-cache-dir \
    hail \
    pandas \
    matplotlib \
    seaborn \
    pyyaml

# Install HTSlib (required by grafanc)
RUN git clone --recursive https://github.com/samtools/htslib.git /tmp/htslib && \
    cd /tmp/htslib && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

RUN echo "/usr/local/lib" > /etc/ld.so.conf.d/local.conf && ldconfig

# Install grafanc tool from binary
COPY grafanc /usr/local/bin/grafanc
RUN chmod +x /usr/local/bin/grafanc

# Copy the pipeline scripts into the container
COPY vcf-af-pipeline.py .
COPY M1_preprocessing.py .
COPY M2_unrelated_samples.py .
COPY M3_ancestry.py .
COPY M4_af_annotation.py .
COPY utils.py .

# Copy any additional files 
COPY GrafAnc_SNPs/ ./GrafAnc_SNPs/
# Make sure GrafAnc data file is exactly where GrafAnc expects it
RUN mkdir -p /app/data
COPY AncSnpPopAFs.txt /app/data/AncSnpPopAFs.txt
RUN mkdir -p /usr/local/cpp/data \
 && ln -sf /app/data/AncSnpPopAFs.txt /usr/local/cpp/data/AncSnpPopAFs.txt