# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this repository is

`add_illumina_readgroups.py` rewrites every read's readgroup in a BAM from its Illumina read name,
and puts matching `@RG` lines in the header. It runs as the `repair_readgroups` task of the GDC
realignment pipeline (`../alignment_pipeline_wolF`), on BAMs from external centers whose
readgroups are missing or wrong.

- **The readgroup rule:** `readgroup_id_from_name()` makes the ID `<flowcell>.<lane>` from the
  3rd and 4th `:`-separated fields of `instrument:run:flowcell:lane:tile:x:y`, or the 1st and
  2nd of a 5-field name without instrument and run (`flowcell:lane:tile:x:y`). Each `@RG` line
  gets `ID`, `CN`, `LB`, `SM` (from `-c`, `-l`, `-s`) and `PL:ILLUMINA`. The rest of the input
  header is kept.
- **Parallelism:** `-p` worker processes, with `-t` compression threads each, tag record-aligned
  chunks of about 1 GiB (a multiple of `-p` of them). Each chunk is written into its own reserved
  region of one sparse `<output>.partial` file. The gaps are then removed in place with
  `fallocate(FALLOC_FL_COLLAPSE_RANGE)` behind the new header; where that's unsupported
  (not ext4/XFS), the chunks are moved down in place instead.
- **Checkpointing:** finished chunks are recorded in `<output>.checkpoint`. A rerun in the same
  directory on the same input (name, size and mtime) resumes from there. The final assembly
  rewrites the partial file in place, so a run interrupted during it starts over.
- **Output:** `<input stem>.readgroupsadded.bam`, in the working directory.
- **README.md:** usage, how it works, and recommended wolF/Slurm job settings with per-job costs
  (300 GB input, read over gcsfuse, written to pd-standard sized at output + 10%).

Remote: `git@github.com:getzlab/add_illumina_readgroups_TOOL.git`, branch `master`.

## Build and deploy

- **Image:** `Dockerfile` installs pysam from source, with htslib's libcurl support, on
  `google/cloud-sdk:slim`. `deploy.sh` builds for `linux/amd64`, and pushes
  `gcr.io/broad-getzlab-workflows/add_illumina_readgroups:v<commit count>` (prefixed with
  `<branch>_` off `master`).
- **The pipeline pins `v17`.** `v13` added `-p`/`-t`, the collapse merge and checkpointing;
  `v17` reads readgroups from names without instrument and run. A new image needs the pin bumped in the pipeline's task (below).
- **Building locally** on an arm64 Mac runs under QEMU emulation: slow for compiles, and not a
  guide to performance. Benchmark on GCE.

## How the pipeline uses it

The pipeline is `../alignment_pipeline_wolF`, on branch `gdc_realigner-wolf-2.0`. Its plan is
`MODERN_WOLF_2.0_PLAN.md` ("Step 2" covers readgroup repair). The task and workflow are in
`wolF/tasks.py`.

**The task** (`class repair_readgroups`) follows the README's recommended settings (job price
is the priority):
- runs `add_illumina_readgroups.py -p "$SLURM_CPUS_PER_TASK" -t 1 ...` on `v17`, and fails
  if `SLURM_CPUS_PER_TASK` is unset;
- `mem` 4G, `use_scratch_disk = True`, `checkpoint = True`;
- per BAM, the workflow sets `cpus-per-task` (about 1 per 22 GB, capped at 16) and
  `scratch_disk_size` (1.1× the input) with `repair_readgroups_sizing()`, and sends the job to
  an n1-highcpu-8/16 partition with `highcpu_slurm_args()`.

If the README's recommendations change, change those helpers too.

**Where it sits in `reversion_workflow`**, for BAM input:

1. **Localize:** `LocalizeToBucket` makes the input BAM available through a gcsfuse bucket
   mount. The task reads it from there, so reads are sequential over gcsfuse.
2. **Repair:** `repair_readgroups` writes the repaired BAM to its scratch disk.
3. **Free the input:** `DeleteLocalizedFiles` (`DeleteLocalizedBAMs`) evicts the localized
   input once the repaired BAM exists.
4. **Read it back:** `bam_readgroup_to_contents` reads its header (`samtools view -H`).
   `biobambam_bamtofastq` reverts it to per-readgroup FASTQs (`outputperreadgroup=1,
   collate=1`).
5. **Free the repaired BAM:** `DeleteDisk` (`DeleteRegroupedBAM`) deletes its scratch disk once
   the FASTQs exist.

**So the repaired BAM is an intermediate.** `bamtofastq` reads it once, start to finish, and it's
then deleted. Nothing downstream needs it sorted, indexed, or kept. What downstream does need:
- **Every read's `RG` tag names an `@RG` line in the header.** The FASTQ file names come from
  the readgroup IDs. The pipeline fails with a clear error if a FASTQ's readgroup isn't in the
  header, and reads without an `RG` tag go to `default_*` FASTQs.
- **Each `@RG` line has the fields used to build bwa's `-R`:** `ID`, `SM`, `LB`, `PL`, `CN`.
- **A valid BAM** that biobambam2 2.0.95 (libmaus2 2.0.610) can read.

**When it runs**, set by `repair_readgroups` on `revert_and_or_align_workflow`:
- **`False` (the default):** never. The BAM's readgroups are used as they are.
- **`True`:** always. `center_name` is required, and `library_name` defaults to `sample_name`.
- **`"auto"`:** the `check_readgroups` task decides. It reads the header and the first 1M reads
  straight from the bucket, and repairs if:
  - the header has no `@RG` lines;
  - a read has no `RG` tag, or names a readgroup the header lacks;
  - one readgroup holds reads from more than one `flowcell.lane`, e.g. a single readgroup
    added over several lanes.

  A non-coordinate-sorted BAM is always repaired, because its first reads could all be from one
  lane. Non-Illumina read names are an error.
  - **Keep the two rules identical.** `check_readgroups` must stay much cheaper than this tool,
    and it must derive readgroups with the same rule. If `readgroup_id_from_name()` changes,
    change the check's awk too.
- **Non-Illumina read names:** this tool has no guard. With fewer than 4 `:`-separated fields,
  the ID comes out empty or meaningless. The `"auto"` check rejects such BAMs, but `True` mode
  doesn't.

## Running inside a wolF task

wolF dispatches tasks through canine onto an autoscaling Slurm cluster on GCE
(`slurm_gcp_docker`). Sources are in `../wolf_update_repos/` (`wolF/` on `combine-2.0`,
`canine/` on `parallel-localization-fuse`, `slurm_gcp_docker/` on `fuse-localize`), each with
its own CLAUDE.md. What matters for this tool:

- **Containers:** each task runs under podman. wolF moves the container's processes into the
  Slurm job's cgroups, so the job is confined to its cores and `mem` (page cache included). But
  `os.cpu_count()` and the memory the tool can see are still the node's. Size `-p`/`-t`
  from `SLURM_CPUS_PER_TASK`, and memory from `SLURM_MEM_PER_NODE` (MB), at run time. wolF passes
  both into the container. A size written into the script doesn't follow wolF's out-of-memory
  retry, which resubmits with 1.5× `mem`. wolF's `UPGRADE_PLAN.md` ("How-to: sizing tool memory
  and threads from the job's environment") has the convention.
- **Node types:** Slurm places the job on a node type with enough CPUs and memory, and a job
  can take only part of a node. A job that fills its node pays for all of it, including unused
  memory, so this tool belongs on n1-highcpu or n1-custom nodes with minimal memory. Node types
  come from the run's `nodetypes.json`. The pipeline's includes n1-highcpu-8/16, preemptible
  and not, for this task. On n1-standard-8, an 8-core job would cost about 30% more.
  `slurm_gcp_docker` also accepts n1-custom types, e.g. 14 vCPUs with minimal memory, at about a
  5% premium per vCPU and GB.
- **Preemption:** nodes are preemptible by default. A preempted job is requeued and its script
  reruns from the start. The tool can resume from its own checkpoint only if the task sets
  wolF's `checkpoint = True`. Without it, wolF prepends `rm -rf *` to the script
  (`wolf/task.py`, after hashing, so it doesn't affect job avoidance), which wipes the partial
  output each attempt.
  - With it on, the job finds its files again because it gets the same scratch disk back: wolF
    names the disk from the job's hashes, and canine reattaches an unfinished disk without
    reformatting it. Not if a task sets a non-deterministic `scratch_disk_name`.
- **Exit codes canine reserves:** canine's job entrypoint treats exit **5** as "requeue at once
  and resume", and **15** as "skip, mark done". The tool must never exit with either by accident.
  Any other nonzero exit fails the shard, and wolF then retries the task (3 times by default).
- **Input:** a path in a gcsfuse-mounted bucket. Reads go over gcsfuse, sequentially, so
  readahead and request size matter.
- **Disks:** with `use_scratch_disk = True`, the output goes to a GCE persistent disk created
  for the job, `pd-standard` by default (`scratch_disk_type="ssd"` is available).
  - **Size:** it starts at `scratch_disk_size` GB, 10 by default. A daemon (canine
    `localization/base.py`, `_disk_resize_daemon_lines`) grows it by 60% whenever free space
    falls below 30%.
  - **Throughput:** pd-standard read and write share a cap that scales with disk size. So a
    disk left to grow from 10 GB is slow early on. To get the README's throughput, set
    `scratch_disk_size` large enough from the start.
  - **Filesystem:** ext4 (same `mkfs` flags as the benchmarks), so the collapse merge applies.
  - **One disk per job:** jobs packed on one node each have their own disk, so they don't share
    disk throughput.
  - **Job avoidance:** a scratch disk labelled finished makes a rerun of the job skip it.
- **Memory:** the tool needs very little (under 125 MB at `-p 32`). The job's `mem` mostly
  decides node placement and page cache. The 4G recommendation is untested under a cap.

## Conventions

- **Protected data:** BAMs and sample metadata are protected. Look only at schema, counts and
  types, never values, and keep copies only in private buckets (public access prevention
  enforced).
- **GCP project:** test and benchmark resources for the realignment pipeline go in
  `getzlab-cmi-pegs-testing`, not the gcloud default. Exception: benchmarks of this tool on its
  own go in `broad-getzlab-workflows`.
- **Commits:**
  - Scan before every commit and push, gated with `&&`:
    `git diff --cached | python3 /Users/dheiman/Documents/workspace/credscan.py && git commit ...`.
    It prints the number of lines scanned; 0 means the scan failed.
  - End commit messages with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - Commit and push only when asked.
- **Local shell:** zsh, where `$VAR` isn't word-split and `grep` is aliased to ugrep. Use
  `/usr/bin/grep` and arrays in scripts.
- **Times:** report clock times in Eastern time.
