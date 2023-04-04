#! /usr/bin/env python3

from argparse import ArgumentParser
from collections import namedtuple
from pathlib import Path
from sys import version_info
from multiprocessing import Pool, Manager, Process
import pysam

if version_info < (3, 7):
    raise RuntimeError("This tool requires Python 3.7+")

def process_read(read, readgroups, ReadGroup, read_queue):
    readgroup_id = '.'.join(read.query_name.split(':')[2:4])
    read.set_tag('RG', readgroup_id, 'Z')
    if readgroup_id not in readgroups:
        readgroups[readgroup_id] = ReadGroup(readgroup_id)
    read_queue.put(read)

def write_reads(read_queue, bam_fh):
    while True:
        read = read_queue.get()
        if read is None:
            break
        bam_fh.write(read)
        
def run_serial(in_bam, temp_bam, ReadGroup, threads):
    readgroups = {}
    # Add readgroup tag to each read and create readgroups for header
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False, threads=threads) as in_bam_fh, \
         pysam.AlignmentFile(temp_bam, 'wb', template=in_bam_fh, threads=threads) as temp_bam_fh:
        header_dict = in_bam_fh.header.to_dict()
        for read in in_bam_fh.fetch(until_eof=True):
            readgroup_id = '.'.join(read.query_name.split(':')[2:4])
            read.set_tag('RG', readgroup_id, 'Z')
            temp_bam_fh.write(read)
            if readgroup_id not in readgroups:
                readgroups[readgroup_id] = ReadGroup(readgroup_id)
    header_dict['RG'] = [readgroup._asdict() for readgroup in readgroups.values()]
    return header_dict

if __name__ == '__main__':
    parser = ArgumentParser(description='Add readgroups based on Illumina ' +
                            'IDs to given BAM')
    parser.add_argument('-c', '--cn', required=True,
                        help='Sequencing Center Name')
    parser.add_argument('-l', '--lb', required=True, help='Library')
    parser.add_argument('-s', '--sm', required=True, help='Sample')
    parser.add_argument('-p', '--proc', type=int, default=3, help='Processes')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Compression Threads')
    parser.add_argument('bam', help='Input BAM')

    args = parser.parse_args()

    assert args.proc > 0, f"There must be at least 1 process, {args.proc} given."
    assert args.threads > 0, f"There must be at least 1 compression thread, {args.threads} given."

    ReadGroup = namedtuple('ReadGroup', 'ID CN LB PL SM',
                           defaults=[args.cn, args.lb, 'ILLUMINA', args.sm])

    in_bam = args.bam
    in_bam_stem = Path(in_bam).stem
    threads = args.threads

    # Intermediate BAM with readgroup tagged reads, but no readgroups in header
    temp_bam = in_bam_stem + '.temp.bam'

    # Header file with readgroups
    header_sam = in_bam_stem + '.header.sam'

    # Final BAM
    final_bam = in_bam_stem + '.readgroupsadded.bam'

    header_dict = {}

    if args.proc > 3:
        with Manager() as manager:
            readgroups = manager.dict()
            read_queue = manager.Queue()

            # Add readgroup tag to each read and create readgroups for header
            with pysam.AlignmentFile(in_bam, 'rb', check_sq=False, threads=threads) as in_bam_fh, \
                 pysam.AlignmentFile(temp_bam, 'wb', template=in_bam_fh, threads=threads) as temp_bam_fh:
                writer = Process(target=write_reads, args=(read_queue, temp_bam_fh))
                header_dict = in_bam_fh.header.to_dict()
                writer.start()
                read_procs = Pool(args.proc - 2)
                read_procs.starmap_async(process_read,
                                         ((read, readgroups,
                                           ReadGroup, read_queue) \
                                          for read in \
                                          in_bam_fh.fetch(until_eof=True)))
                read_procs.close()
                read_procs.join()

                read_queue.put(None)
                writer.join()

            # Add readgroups to header
            header_dict['RG'] = [readgroup._asdict() for readgroup in readgroups.values()]
    else:
        header_dict = run_serial(in_bam, temp_bam, ReadGroup, threads)

    # Write updated header to file
    header_sam_fh = pysam.AlignmentFile(header_sam, 'wh', header=header_dict)
    header_sam_fh.close()

    # Create new bam with readgroups in header and read tags
    Path(final_bam).touch()
    pysam.reheader(header_sam, temp_bam, save_stdout=final_bam)

    # Delete intermediate BAM and header file
    Path(temp_bam).unlink()
    Path(header_sam).unlink()
