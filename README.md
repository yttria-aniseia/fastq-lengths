# Introduction
`fastq-lengths` is a simple utility for quickly determining a limited set of
summary statistics related to the length of sequences in a FASTQ file.

`fastq-lengths` generally assumes well-formed fastq and does little validation,
but notably does *not* assume sequence and quality entries to be on one line
like some tools -- if the FASTQ file conforms to
[spec](https://maq.sourceforge.net/fastq.shtml) it should work.

No external dependencies other than functions in the POSIX C standard.  
It is approximately 12x faster to find the median sequence length than the
naive equivalent using Python and SeqIO, with runtime mem usage in the KB (vs
tens of MB for Python).

```
usage:
fastq-lengths [subcommand [stopafter]] file.fastq
    subcommand: lengths | median | count
        lengths - report all sequence lengths and occurence counts (default)
        median  - report the median sequence length
        count   - report the total number of records in the file

    stopafter:  only read the first <stopafter> records in the file, e.g.
                when 1000 records is enough to characterize the contents.
```
```
fastq-lengths lengths 100 example.fq
32  1
150 99
```
```
fastq-lengths median example.fq
150
```
```
fastq-lengths count example.fq
45101612
```


# Build
```sh
git clone git@github.com:yttria-aniseia/fastq-lengths.git
cd fastq-lengths
make
```


```
    Command being timed: "./fastq-lengths median ERR1757416_2.fastq"
    User time (seconds): 0.45
    System time (seconds): 0.00
    Percent of CPU this job got: 100%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.45
    Maximum resident set size (kbytes): 1436
```