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
- **Threads:** with `-t 1`, each worker is effectively single-threaded and uses
  about one core. With `-t N` > 1, each worker gets two htslib thread pools of
  N threads, one decompressing its input and one compressing its output.
- **Checkpointing:** finished chunks are recorded in `<output>.checkpoint`. A
  rerun in the same directory on the same input, e.g. after preemption,
  resumes from there, losing at most about one chunk of work per worker.

## Recommended settings (wolF / Slurm on GCP)

- **`cpus-per-task`:** about 1 per 22 GB of input, e.g. 14 for a 300 GB BAM.
  That's about where the output disk's write throughput runs out (below). Run
  with `-p "$SLURM_CPUS_PER_TASK" -t 1`. Inside a wolF task, `os.cpu_count()`
  shows the whole node, but the job is confined to its own cores: wolF moves
  the container into the Slurm job's cgroups (cpuset and memory), and Slurm
  constrains both.
- **Why `-t 1`:** more workers beat more threads, and once there's a worker
  per core, extra threads only compete for the job's own cores. On the 300 GB
  BAM, `-t 2` was slower than `-t 1` at both sizes tested: 203.8 min against
  about 195 min on 8 cores, and 121.4 min against 109.5 min on 16.
- **`mem`:** 4 GB. The tool's own memory is tiny, about 3 MB per worker (under
  100 MB at `-p 32`). But Slurm caps the job's memory cgroup at `mem`, and the
  page cache from writing the output is charged to it. The benchmarks ran
  without a cap, so tighter settings are untested. 4 GB adds about $0.03 to a
  300 GB job.
- **`scratch_disk_size`:** the expected output size + 10%, about 1.1× the input
  size, set from the start. pd-standard throughput scales with disk size, so a
  disk left to grow from the 10 GB default is slow early on.
- **`checkpoint = True`** on the wolF task. Otherwise wolF clears the working
  directory (`rm -rf *`) before each attempt, and a preempted job starts over.

**Why:** the output disk's write throughput sets how many cores are worth
paying for. At output + 10%, pd-standard takes about 43–47 MB/s for a 300 GB
BAM, which is what about 14 cores produce. Below that, a job uses the same
core-hours at any size but holds its disk longer. Above it, extra cores sit
idle. Both the disk and its throughput grow with the input, so the best size
scales with input size. 14 cores itself wasn't run; it's interpolated from the
8-, 16- and 32-core runs.

**Node type:** a job that fills its node pays for the whole node, including
memory it doesn't use. So run it on a node with little memory per core:
n1-highcpu, or n1-custom with the minimum 0.9 GB per vCPU. On n1-standard-8,
the same 8-core job costs about 30% more.

### Benchmarks

300 GB synthetic BAM, read over gcsfuse, written to a 313 GiB pd-standard disk
(output + 10%), in the Docker image, on n1 spot VMs in us-central1, with
`-t 1`. Costs are for the whole node at n1 spot prices (October 2026), plus
pd-standard at $0.04/GB-month.

| node                       | `-p` | wall      | cores used | job cost |
|----------------------------|------|-----------|------------|----------|
| n1-highcpu-8               | 8    | ~3.25 h*  | 7.5        | ~$0.51   |
| n1-custom-14 (12.75 GB)    | 14   | ~1.9 h**  | –          | ~$0.52   |
| n1-highcpu-16              | 16   | 1.83 h    | 13.7       | ~$0.54   |
| n1-highcpu-32              | 32   | 1.78 h    | 14.2       | ~$1.02   |
| n1-standard-8 (30 GB)      | 8    | ~3.25 h*  | 7.5        | ~$0.66   |

\* Projected from the first 2 hours of a run that was later preempted.
\** Not run; interpolated.

- **Fastest for the money:** 14–16 cores costs about 3–6% more than 8 and
  finishes about 1.7× sooner. 32 cores is disk-bound: no faster, and twice the
  cost.
- **On a shared node,** where Slurm fills the rest of the node with other
  jobs, a job costs only its share: about $0.01544/vCPU-hour plus
  $0.00207/GB-hour of `mem`, about $0.46–0.50 for 8–16 cores.
- **Preemption:** a 16-core run preempted halfway resumed with `Resuming: 144
  of 288 chunks already done`, and its output passed a full read-back check.
- **Checkpointing overhead:** about 2% of run time.
