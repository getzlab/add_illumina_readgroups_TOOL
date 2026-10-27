#! /usr/bin/env python3

import re
import struct
import zlib
from argparse import ArgumentParser
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from sys import version_info

import pysam

if version_info < (3, 7):
    raise RuntimeError("This tool requires Python 3.7+")

BGZF_MAGIC = b'\x1f\x8b\x08\x04'
BGZF_HEADER_LEN = 18
MAX_BGZF_BLOCK = 65536
MAX_RECORD_SIZE = 1 << 24
QNAME_RE = re.compile(rb'[!-?A-~]+')
SEARCH_WINDOW = 4 * MAX_BGZF_BLOCK


def readgroup_id_from_name(query_name):
    return '.'.join(query_name.split(':')[2:4])


def read_bgzf_block(fh, coffset):
    """Return (block_size, decompressed_data) for the BGZF block at coffset,
    or None if there is no valid BGZF block there."""
    fh.seek(coffset)
    header = fh.read(BGZF_HEADER_LEN)
    if len(header) < BGZF_HEADER_LEN or header[:4] != BGZF_MAGIC:
        return None
    xlen = struct.unpack_from('<H', header, 10)[0]
    if xlen != 6 or header[12:14] != b'BC' or header[14:16] != b'\x02\x00':
        return None
    bsize = struct.unpack_from('<H', header, 16)[0] + 1
    rest = fh.read(bsize - BGZF_HEADER_LEN)
    if len(rest) != bsize - BGZF_HEADER_LEN:
        return None
    try:
        data = zlib.decompress(rest[:-8], -15)
    except zlib.error:
        return None
    if len(data) != struct.unpack_from('<I', rest, len(rest) - 4)[0]:
        return None
    return bsize, data


def is_record_chain(buf, pos, n_ref):
    """Check that a chain of valid BAM records starts at buf[pos]."""
    n_valid = 0
    while pos + 36 <= len(buf):
        (block_size, ref_id, ref_pos, l_read_name, _mapq, _bin, n_cigar,
         _flag, l_seq, next_ref_id, next_pos, _tlen) = \
            struct.unpack_from('<iiiBBHHHiiii', buf, pos)
        if not (32 < block_size < MAX_RECORD_SIZE
                and -1 <= ref_id < n_ref and -1 <= next_ref_id < n_ref
                and ref_pos >= -1 and next_pos >= -1
                and l_read_name > 1 and l_seq >= 0
                and 32 + l_read_name + 4 * n_cigar + (l_seq + 1) // 2 + l_seq
                <= block_size):
            return False
        name_end = pos + 36 + l_read_name
        if name_end > len(buf):
            break
        if buf[name_end - 1] != 0 or \
                not QNAME_RE.fullmatch(buf, pos + 36, name_end - 1):
            return False
        n_valid += 1
        pos += 4 + block_size
    return n_valid > 0


def find_record_start(fh, approx_offset, n_ref, file_size):
    """Find the virtual offset of the first BAM record starting in the first
    BGZF block at or after approx_offset. Returns None if there is none."""
    fh.seek(approx_offset)
    window = fh.read(SEARCH_WINDOW)
    search_from = 0
    while True:
        idx = window.find(BGZF_MAGIC, search_from)
        if idx < 0:
            return None
        coffset = approx_offset + idx
        block = read_bgzf_block(fh, coffset)
        if block is None:
            search_from = idx + 1
            continue
        bsize, data = block
        # Pull in a few following blocks so we can validate records that
        # span block boundaries
        buf = bytearray(data)
        next_coffset = coffset + bsize
        for _ in range(3):
            if next_coffset >= file_size:
                break
            nxt = read_bgzf_block(fh, next_coffset)
            if nxt is None:
                break
            buf.extend(nxt[1])
            next_coffset += nxt[0]
        if not data:
            # Empty block (e.g. EOF marker); try the next one
            search_from = idx + bsize
            if search_from >= len(window):
                return None
            continue
        for uoffset in range(len(data)):
            if is_record_chain(buf, uoffset, n_ref):
                return (coffset << 16) | uoffset
        # No record starts in this block; move on to the next block
        search_from = idx + bsize
        if search_from >= len(window):
            approx_offset = coffset + bsize
            fh.seek(approx_offset)
            window = fh.read(SEARCH_WINDOW)
            search_from = 0


def compute_chunks(in_bam, n_chunks):
    """Split the BAM into up to n_chunks [start, end) virtual offset ranges,
    each beginning on a record boundary. end is None for the last chunk."""
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as bam_fh:
        first_record = bam_fh.tell()
        n_ref = bam_fh.nreferences
    file_size = Path(in_bam).stat().st_size
    starts = [first_record]
    with open(in_bam, 'rb') as raw_fh:
        for i in range(1, n_chunks):
            voffset = find_record_start(raw_fh, file_size * i // n_chunks,
                                        n_ref, file_size)
            if voffset is not None and voffset > starts[-1]:
                starts.append(voffset)
    return list(zip(starts, starts[1:] + [None]))


def process_chunk(in_bam, out_bam, start, end, threads):
    """Tag reads in [start, end) and write them to out_bam. Returns the
    readgroup IDs encountered, in order of first appearance."""
    readgroups = {}
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False, threads=threads) as in_bam_fh, \
         pysam.AlignmentFile(out_bam, 'wb', template=in_bam_fh, threads=threads) as out_bam_fh:
        in_bam_fh.seek(start)
        while end is None or in_bam_fh.tell() < end:
            try:
                read = next(in_bam_fh)
            except StopIteration:
                break
            readgroup_id = readgroup_id_from_name(read.query_name)
            read.set_tag('RG', readgroup_id, 'Z')
            out_bam_fh.write(read)
            readgroups[readgroup_id] = None
    return list(readgroups)


if __name__ == '__main__':
    parser = ArgumentParser(description='Add readgroups based on Illumina ' +
                            'IDs to given BAM')
    parser.add_argument('-c', '--cn', required=True,
                        help='Sequencing Center Name')
    parser.add_argument('-l', '--lb', required=True, help='Library')
    parser.add_argument('-s', '--sm', required=True, help='Sample')
    parser.add_argument('-p', '--proc', type=int, default=1,
                        help='Worker processes')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Compression threads per worker process')
    parser.add_argument('bam', help='Input BAM')

    args = parser.parse_args()

    assert args.proc > 0, f"There must be at least 1 process, {args.proc} given."
    assert args.threads > 0, f"There must be at least 1 compression thread, {args.threads} given."

    ReadGroup = namedtuple('ReadGroup', 'ID CN LB PL SM',
                           defaults=[args.cn, args.lb, 'ILLUMINA', args.sm])

    in_bam = args.bam
    in_bam_stem = Path(in_bam).stem

    # Header file with readgroups
    header_sam = in_bam_stem + '.header.sam'

    # Final BAM
    final_bam = in_bam_stem + '.readgroupsadded.bam'

    chunks = compute_chunks(in_bam, args.proc)
    # Intermediate BAMs with readgroup tagged reads, but no readgroups in header
    temp_bams = [f'{in_bam_stem}.temp.{idx}.bam' for idx in range(len(chunks))]

    if len(chunks) == 1:
        chunk_readgroups = [process_chunk(in_bam, temp_bams[0], *chunks[0],
                                          args.threads)]
    else:
        with ProcessPoolExecutor(max_workers=len(chunks)) as executor:
            futures = [executor.submit(process_chunk, in_bam, temp_bam,
                                       start, end, args.threads)
                       for temp_bam, (start, end) in zip(temp_bams, chunks)]
            chunk_readgroups = [future.result() for future in futures]

    readgroups = {}
    for readgroup_ids in chunk_readgroups:
        for readgroup_id in readgroup_ids:
            readgroups.setdefault(readgroup_id, ReadGroup(readgroup_id))

    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as in_bam_fh:
        header_dict = in_bam_fh.header.to_dict()

    # Add readgroups to header
    header_dict['RG'] = [readgroup._asdict() for readgroup in readgroups.values()]

    # Write updated header to file
    header_sam_fh = pysam.AlignmentFile(header_sam, 'wh', header=header_dict)
    header_sam_fh.close()

    # Concatenate tagged chunks under the new header (copies compressed
    # blocks without recompressing)
    pysam.cat('-h', header_sam, '-o', final_bam, *temp_bams)

    # Delete intermediate BAMs and header file
    for temp_bam in temp_bams:
        Path(temp_bam).unlink()
    Path(header_sam).unlink()
