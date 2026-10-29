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

The output is `<input stem>.readgroupsadded.bam` in the working directory.

## How it works

- The input is split into record-aligned chunks of about 1 GiB, a multiple of
  `-p` of them. `-p` worker processes tag them, each with `-t` compression
  threads.
- Each chunk is written straight into its own reserved region of one sparse
  file (`<output>.partial`). The gaps between them are then removed in place
  with `fallocate(FALLOC_FL_COLLAPSE_RANGE)`, behind the final header with the
  new `@RG` lines. On ext4/XFS this takes seconds and copies no data, so peak
  disk use is about the output size. Elsewhere (e.g. NFS, macOS) the chunks are
  moved down in place instead, which is slower.
- **Checkpointing:** finished chunks are recorded in `<output>.checkpoint`. A
  rerun in the same directory on the same input, e.g. after preemption,
  resumes from there, losing at most about one chunk of work per worker.

## Recommended settings (wolF / Slurm on GCP)

- **`cpus-per-task`:** about 1 per 40 GB of input, e.g. 8 for a 300 GB BAM. Run
  with `-p "$SLURM_CPUS_PER_TASK" -t 1`. Inside a wolF task, `os.cpu_count()`
  shows the whole node, not the job's share.
- **`mem`:** 1 GB is plenty. The tool's own memory is about 3 MB per worker
  (under 100 MB at `-p 32`); the rest is reclaimable page cache.
- **`scratch_disk_size`:** the expected output size + 10%, about 1.1× the input
  size, set from the start. pd-standard throughput scales with disk size, so a
  disk left to grow from the 10 GB default is slow early on.
- **`checkpoint = True`** on the wolF task. Otherwise wolF clears the working
  directory (`rm -rf *`) before each attempt, and a preempted job starts over.

**Why:** the output disk's write throughput sets how many cores are worth
paying for. At output + 10%, pd-standard takes about 43–47 MB/s for a 300 GB
BAM, which is what about 13 cores produce. Below that, a job uses the same
core-hours at any size but holds its disk longer. Above it, extra cores sit
idle. Both the disk and its throughput grow with the input, so the best size
scales with input size.

### Benchmarks

300 GB synthetic BAM, read over gcsfuse, written to a 313 GiB pd-standard disk
(output + 10%), in the Docker image, on n1 spot VMs in us-central1. Costs are
for the job's share of a node: $0.01544/vCPU-hour and $0.00207/GB-hour of
memory (n1 spot, October 2026), with `mem` 1 GB, plus pd-standard at
$0.04/GB-month.

| `cpus-per-task` | wall     | cores used | job cost |
|-----------------|----------|------------|----------|
| 8               | ~3.25 h* | 7.5        | ~$0.47   |
| 16              | 1.83 h   | 13.7       | ~$0.49   |
| 32              | 1.78 h   | 14.2       | ~$0.92   |

\* Projected from the first 2 hours of a run that was later preempted.

- **Fastest for the money:** 16 cores costs about 4% more than 8 and finishes
  1.8× sooner. 32 cores is disk-bound: no faster, and twice the cost.
- **Preemption:** a 16-core run preempted halfway resumed with `Resuming: 144
  of 288 chunks already done`, and its output passed a full read-back check.
- **Checkpointing overhead:** about 2% of run time.
