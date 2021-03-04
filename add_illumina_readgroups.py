#! /usr/bin/env python3

import pysam
from argparse import ArgumentParser
from collections import namedtuple
from pathlib import Path
from sys import version_info

if version_info < (3, 7):
    raise RuntimeError("This tool requires Python 3.7+")

def main():
    parser = ArgumentParser(description='Add readgroups based on Illumina ' +
                            'IDs to given BAM')
    parser.add_argument('-c', '--cn', required=True,
                        help='Sequencing Center Name')
    parser.add_argument('-l', '--lb', required=True, help='Library')
    parser.add_argument('-s', '--sm', required=True, help='Sample')
    parser.add_argument('bam', help='Input BAM')

    args = parser.parse_args()

    ReadGroup = namedtuple('ReadGroup', 'ID CN LB PL SM',
                           defaults=[args.cn, args.lb, 'ILLUMINA', args.sm])

    in_bam = args.bam
    in_bam_stem = Path(in_bam).stem

    # Intermediate BAM with readgroup tagged reads, but no readgroups in header
    temp_bam = in_bam_stem + '.temp.bam'

    # Header file with readgroups
    header_sam = in_bam_stem + '.header.sam'

    # Final BAM
    final_bam = in_bam_stem + '.readgroupsadded.bam'

    readgroups = dict()
    header_dict = dict()

    # Add readgroup tag to each read and create readgroups for header
    with pysam.AlignmentFile(in_bam, 'rb') as in_bam_fh, \
         pysam.AlignmentFile(temp_bam, 'wb', template=in_bam_fh) as temp_bam_fh:
        header_dict = in_bam_fh.header.to_dict()
        for read in in_bam_fh.fetch(until_eof=True):
            readgroup_id = '.'.join(read.query_name.split(':')[2:4])
            read.set_tag('RG', readgroup_id, 'Z')
            temp_bam_fh.write(read)
            if readgroup_id not in readgroups:
                readgroups[readgroup_id] = ReadGroup(readgroup_id)

    # Add readgroups to header
    header_dict['RG'] = [readgroup._asdict() for readgroup in readgroups.values()]
    # Write updated header to file
    header_sam_fh = pysam.AlignmentFile(header_sam, 'wh', header=header_dict)
    header_sam_fh.close()

    # Create new bam with readgroups in header and read tags
    Path(final_bam).touch()
    pysam.reheader(header_sam, temp_bam, save_stdout=final_bam)

    # Delete intermediate BAM and header file
    Path(temp_bam).unlink()
    Path(header_sam).unlink()

if __name__ == '__main__':
    main()
