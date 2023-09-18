FROM busybox

LABEL base_image="busybox"
LABEL version="2"
LABEL software="fastq-lengths"
LABEL software.version="0.2.8"
LABEL about.summary="quickly characterize fastq sequence lengths"
LABEL about.home="https://github.com/yttria-aniseia/fastq-lengths"

COPY bin/x86_64/fastq-lengths /usr/local/bin/fastq-lengths
CMD ["fastq-lengths" "in.fastq"]
