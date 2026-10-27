# add_illumina_readgroups.py

```
usage: add_illumina_readgroups.py [-h] -c CN -l LB -s SM [-p PROC]
                                  [-t THREADS]
                                  bam

Add readgroups based on Illumina IDs to given BAM

positional arguments:
  bam                   Input BAM

optional arguments:
  -h, --help            show this help message and exit
  -c, --cn CN           Sequencing Center Name
  -l, --lb LB           Library
  -s, --sm SM           Sample
  -p, --proc PROC       Worker processes
  -t, --threads THREADS
                        Compression threads per worker process
```

With `-p N`, the BAM is split into N record-aligned chunks that are processed
in parallel (each with `-t` compression threads), then concatenated.
