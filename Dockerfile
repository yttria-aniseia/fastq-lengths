FROM scratch

LABEL base_image="scratch"
LABEL version="1"
LABEL software="fastq-lengths"
LABEL software.version="0.1.0"
LABEL about.summary="quickly characterize fastq sequence lengths"
LABEL about.home="https://github.com/yttria-aniseia/fastq-lengths"

ENV PATH=/
COPY bin/x86_64/fastq-lengths /fastq-lengths
CMD ["fastq-lengths" "in.fastq"]