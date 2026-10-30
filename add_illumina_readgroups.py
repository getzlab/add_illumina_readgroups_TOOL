#! /usr/bin/env python3

import ctypes
import ctypes.util
import json
import mmap
import os
import re
import struct
import sys
import zlib
from argparse import ArgumentParser
from collections import deque, namedtuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
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
BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
COPY_BUFFER = 16 << 20
COPY_DEPTH = 8
DIRECT_ALIGN = 4096
FALLOC_FL_KEEP_SIZE = 0x01
FALLOC_FL_PUNCH_HOLE = 0x02
FALLOC_FL_COLLAPSE_RANGE = 0x08
# Sparse space reserved after each chunk's expected output, and how close a
# worker may get to the next chunk's region before aborting
CHUNK_MARGIN = 64 << 20
CHUNK_GUARD = 32 << 20
MAX_BLOCK_INPUT = 0xff00
# Input bytes per chunk. Finished chunks are checkpointed, so this bounds the
# work a restart (e.g. after preemption) loses to about one chunk per worker.
CHUNK_INPUT = 1 << 30
CHECKPOINT_VERSION = 1


def readgroup_id_from_name(query_name):
    # <flowcell>.<lane> from instrument:run:flowcell:lane:tile:x:y (maybe with
    # a UMI after y), or from flowcell:lane:tile:x:y when instrument and run
    # are missing
    fields = query_name.split(':')
    if len(fields) == 5:
        return fields[0] + '.' + fields[1]
    return '.'.join(fields[2:4])


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
    crc, isize = struct.unpack_from('<II', rest, len(rest) - 8)
    if len(data) != isize or zlib.crc32(data) != crc:
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


def process_chunk(in_bam, out_path, offset, limit, start, end, threads):
    """Tag reads in [start, end) and write them as a BAM at byte offset of
    out_path, aborting if the output would get within CHUNK_GUARD of limit.
    Returns the readgroup IDs encountered, in order of first appearance, and
    the byte offset where the output ends."""
    readgroups = {}
    fd = os.open(out_path, os.O_WRONLY)
    try:
        # Clear output left by an interrupted earlier attempt: htslib's
        # unaligned writes over existing data would read each page back first,
        # which is slow on network disks
        try:
            fallocate(fd, FALLOC_FL_PUNCH_HOLE | FALLOC_FL_KEEP_SIZE,
                      offset, limit - offset)
        except (OSError, AttributeError):
            pass
        os.lseek(fd, offset, os.SEEK_SET)
        # htslib writes through a duplicate of fd, which shares its offset
        with pysam.AlignmentFile(in_bam, 'rb', check_sq=False, threads=threads) as in_bam_fh, \
             pysam.AlignmentFile(os.dup(fd), 'wb', template=in_bam_fh, threads=threads) as out_bam_fh:
            in_bam_fh.seek(start)
            n_reads = 0
            while end is None or in_bam_fh.tell() < end:
                try:
                    read = next(in_bam_fh)
                except StopIteration:
                    break
                readgroup_id = readgroup_id_from_name(read.query_name)
                read.set_tag('RG', readgroup_id, 'Z')
                out_bam_fh.write(read)
                readgroups[readgroup_id] = None
                n_reads += 1
                if n_reads % 65536 == 0 and \
                        os.lseek(fd, 0, os.SEEK_CUR) > limit - CHUNK_GUARD:
                    raise RuntimeError('Chunk output exceeded its reserved space')
        out_end = os.lseek(fd, 0, os.SEEK_CUR)
        # Persist the output before the chunk is checkpointed as done
        os.fsync(fd)
    finally:
        os.close(fd)
    if out_end > limit:
        raise RuntimeError('Chunk output exceeded its reserved space')
    return list(readgroups), out_end


def chunk_count(in_bam_size, n_proc):
    """A multiple of n_proc, so the workers finish together, giving chunks of
    about CHUNK_INPUT bytes."""
    return n_proc * max(1, round(in_bam_size / (n_proc * CHUNK_INPUT)))


def plan_layout(in_bam, temp_bam, n_proc, threads):
    """Split the input into chunks and choose where in temp_bam (which is
    truncated) each chunk's output goes. Returns the layout as a dict."""
    in_bam_size = Path(in_bam).stat().st_size
    chunks = compute_chunks(in_bam, chunk_count(in_bam_size, n_proc))

    open(temp_bam, 'wb').close()
    align = os.statvfs(temp_bam).f_bsize
    if not (DIRECT_ALIGN <= align <= MAX_BGZF_BLOCK and align & (align - 1) == 0):
        align = DIRECT_ALIGN

    # Measure the chunk header size so record blocks start block-aligned
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as in_bam_fh, \
         pysam.AlignmentFile(temp_bam, 'wb', template=in_bam_fh, threads=threads):
        pass
    chunk_header_size = os.path.getsize(temp_bam) - len(BGZF_EOF)
    in_header_size = len(bam_header_bytes(temp_bam))
    open(temp_bam, 'wb').close()

    # Room for the final header, which adds readgroups to the input header
    reserved = 2 * in_header_size + CHUNK_MARGIN
    offsets = []
    for start, end in chunks:
        body_start = align_up(reserved + chunk_header_size, align)
        offsets.append(body_start - chunk_header_size)
        chunk_size = (in_bam_size if end is None else end >> 16) - (start >> 16)
        reserved = body_start + 2 * chunk_size + CHUNK_MARGIN
    limits = offsets[1:] + [reserved]
    return dict(chunks=chunks, offsets=offsets, limits=limits, align=align)


def fsync_dir(path):
    fd = os.open(os.path.dirname(os.path.abspath(path)), os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_checkpoint(checkpoint, layout, done):
    """Atomically write the layout and the finished chunks' results."""
    temp = checkpoint + '.tmp'
    with open(temp, 'w') as fh:
        fh.write(json.dumps(layout) + '\n')
        for index, (readgroup_ids, out_end) in sorted(done.items()):
            fh.write(json.dumps(dict(chunk=index, readgroups=readgroup_ids,
                                     out_end=out_end)) + '\n')
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temp, checkpoint)
    fsync_dir(checkpoint)


def append_checkpoint(checkpoint, entry):
    with open(checkpoint, 'a') as fh:
        fh.write(json.dumps(entry) + '\n')
        fh.flush()
        os.fsync(fh.fileno())


def load_checkpoint(checkpoint, temp_bam, identity):
    """Return (layout, {chunk index: (readgroup IDs, out_end)}) left by an
    earlier run on the same input, or None if there's nothing to resume."""
    try:
        with open(checkpoint) as fh:
            lines = fh.read().splitlines()
        layout = json.loads(lines[0])
    except (OSError, ValueError, IndexError):
        return None
    if layout.get('identity') != identity or not os.path.exists(temp_bam):
        return None
    done = {}
    for line in lines[1:]:
        try:
            entry = json.loads(line)
        except ValueError:
            break  # last line partly written
        if entry.get('assembling'):
            return None  # temp_bam was being rewritten in place
        done[entry['chunk']] = (entry['readgroups'], entry['out_end'])
    return layout, done


def align_up(value, align):
    return -(-value // align) * align


def header_end(fh, offset):
    """Return the byte offset where the BGZF blocks holding the BAM header
    starting at offset end. htslib flushes after writing the header, so the
    header ends on a block boundary. Returns None if it doesn't."""
    data = bytearray()
    coffset = offset
    while True:
        # Header layout: magic, l_text, text, n_ref, (l_name, name, l_ref)*
        pos = 8 + struct.unpack_from('<i', data, 4)[0] if len(data) >= 8 else None
        if pos is not None and len(data) >= pos + 4:
            n_ref = struct.unpack_from('<i', data, pos)[0]
            pos += 4
            for _ in range(n_ref):
                if len(data) < pos + 4:
                    pos = None
                    break
                pos += 8 + struct.unpack_from('<i', data, pos)[0]
            if pos is not None and len(data) >= pos:
                return coffset if len(data) == pos else None
        block = read_bgzf_block(fh, coffset)
        if block is None or not block[1]:
            raise ValueError('Could not parse BAM header')
        coffset += block[0]
        data += block[1]


def bam_header_bytes(header_bam):
    """Return the uncompressed header of a BAM that has no records."""
    data = bytearray()
    coffset = 0
    with open(header_bam, 'rb') as fh:
        while True:
            block = read_bgzf_block(fh, coffset)
            if block is None:
                raise ValueError(f'Could not parse {header_bam}')
            if not block[1]:
                return bytes(data)
            coffset += block[0]
            data += block[1]


def last_bgzf_block(fh, start, end):
    """Return (offset, data) of the BGZF block in [start, end) ending at end."""
    window_start = max(start, end - MAX_BGZF_BLOCK)
    fh.seek(window_start)
    window = fh.read(end - window_start)
    idx = window.find(BGZF_MAGIC)
    while idx >= 0:
        block = read_bgzf_block(fh, window_start + idx)
        if block is not None and window_start + idx + block[0] == end:
            return window_start + idx, block[1]
        idx = window.find(BGZF_MAGIC, idx + 1)
    raise ValueError('Could not find last BGZF block of chunk')


def tail_bgzf_blocks(fh, start, end, min_data=1024):
    """Return (offset, data) of the last BGZF blocks in [start, end), taking
    enough blocks to hold at least min_data bytes where possible, so that the
    data can be re-encoded with padding to any alignment."""
    offset, data = last_bgzf_block(fh, start, end)
    while len(data) < min_data and offset > start:
        offset, more = last_bgzf_block(fh, start, offset)
        data = more + data
    return offset, data


def bgzf_block(cdata, data):
    return (struct.pack('<4sIBBHBBHH', BGZF_MAGIC, 0, 0, 0xff, 6, ord('B'),
                        ord('C'), 2, len(cdata) + 25)
            + cdata + struct.pack('<II', zlib.crc32(data), len(data)))


def bgzf_compress(data):
    out = bytearray()
    for pos in range(0, len(data), MAX_BLOCK_INPUT):
        chunk = data[pos:pos + MAX_BLOCK_INPUT]
        compressor = zlib.compressobj(zlib.Z_DEFAULT_COMPRESSION, zlib.DEFLATED, -15)
        cdata = compressor.compress(chunk) + compressor.flush()
        if len(cdata) + 26 > MAX_BGZF_BLOCK:
            cdata = b'\x01' + struct.pack('<HH', len(chunk), len(chunk) ^ 0xffff) + chunk
        out += bgzf_block(cdata, chunk)
    return bytes(out)


def bgzf_stored_exact(data, length):
    """Encode non-empty data as non-empty BGZF blocks of uncompressed (stored)
    deflate data, padded with empty stored deflate blocks to exactly length
    bytes in total. Returns None if that isn't possible."""
    # length = len(data) + 26 * n_blocks + 5 * n_stored_blocks
    extra = length - len(data)
    n_blocks = extra % 5 or 5
    while n_blocks <= len(data) and 31 * n_blocks <= extra:
        n_empty = (extra - 31 * n_blocks) // 5
        if 31 + -(-len(data) // n_blocks) + 5 * -(-n_empty // n_blocks) <= MAX_BGZF_BLOCK:
            out = bytearray()
            pos = 0
            for i in range(n_blocks):
                size = len(data) // n_blocks + (i < len(data) % n_blocks)
                chunk = data[pos:pos + size]
                pos += size
                empties = n_empty // n_blocks + (i < n_empty % n_blocks)
                out += bgzf_block(b'\x00\x00\x00\xff\xff' * empties + b'\x01'
                                  + struct.pack('<HH', size, size ^ 0xffff) + chunk,
                                  chunk)
            return bytes(out)
        n_blocks += 5
    return None


def stored_to_boundary(offset, data, align):
    """Encode data as BGZF blocks that, written at offset, end on a multiple
    of align."""
    end = align_up(offset + len(data) + 31, align)
    for _ in range(64):
        blocks = bgzf_stored_exact(data, end - offset)
        if blocks is not None:
            return blocks
        end += align
    raise ValueError('Could not pad BGZF data to alignment')


def encode_header(header, align):
    """BGZF-compress an uncompressed BAM header so it ends on a multiple of
    align."""
    head = bgzf_compress(header[:-MAX_BLOCK_INPUT])
    return head + stored_to_boundary(len(head), header[-MAX_BLOCK_INPUT:], align)


def fallocate(fd, mode, offset, length):
    libc = ctypes.CDLL(ctypes.util.find_library('c'), use_errno=True)
    func = libc.fallocate  # AttributeError if unavailable
    func.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int64, ctypes.c_int64]
    if func(fd, mode, offset, length) != 0:
        err = ctypes.get_errno()
        raise OSError(err, os.strerror(err))


def collapse_range(fd, offset, length):
    """Remove [offset, offset + length) from the file without copying, where
    the filesystem supports it (e.g. ext4, XFS)."""
    fallocate(fd, FALLOC_FL_COLLAPSE_RANGE, offset, length)


def open_for_copy(src, start):
    """Open src with O_DIRECT if the platform and filesystem support it.
    Returns (fd, direct)."""
    if hasattr(os, 'O_DIRECT'):
        fd = None
        buf = mmap.mmap(-1, DIRECT_ALIGN)
        try:
            fd = os.open(src, os.O_RDONLY | os.O_DIRECT)
            with memoryview(buf) as view:
                os.preadv(fd, [view], start - start % DIRECT_ALIGN)
            return fd, True
        except OSError:
            if fd is not None:
                os.close(fd)
        finally:
            buf.close()
    return os.open(src, os.O_RDONLY), False


def read_full(fd, buf, offset):
    """Read into buf from offset until buf is full or EOF; return bytes read."""
    total = 0
    with memoryview(buf) as view:
        while total < len(buf):
            with view[total:] as rest:
                n = os.preadv(fd, [rest], offset + total)
            if n == 0:
                break
            total += n
    return total


def move_range(path, fd, start, end, dst):
    """Copy bytes [start, end) of path to offset dst <= start via fd. Uses
    several concurrent large O_DIRECT reads where supported, since small
    buffered sequential reads are latency-bound on network disks."""
    read_fd, direct = open_for_copy(path, start)
    offsets = iter(range(start - start % DIRECT_ALIGN if direct else start,
                         end, COPY_BUFFER))
    buffers = [mmap.mmap(-1, COPY_BUFFER) for _ in range(COPY_DEPTH)]  # page-aligned
    try:
        with ThreadPoolExecutor(COPY_DEPTH) as pool:
            pending = deque()
            for buf in buffers:
                offset = next(offsets, None)
                if offset is not None:
                    pending.append((offset, buf, pool.submit(read_full, read_fd, buf, offset)))
            while pending:
                offset, buf, future = pending.popleft()
                skip = max(start - offset, 0)
                n = min(future.result(), end - offset)
                if n <= skip:
                    raise IOError(f'Unexpected end of {path}')
                with memoryview(buf) as view, view[skip:n] as data:
                    written = 0
                    while written < len(data):
                        with data[written:] as rest:
                            written += os.pwrite(fd, rest, dst + written)
                dst += n - skip
                offset = next(offsets, None)
                if offset is not None:
                    pending.append((offset, buf, pool.submit(read_full, read_fd, buf, offset)))
    finally:
        os.close(read_fd)
        for buf in buffers:
            buf.close()


def pwrite_all(fd, data, offset):
    with memoryview(data) as view:
        written = 0
        while written < len(view):
            with view[written:] as rest:
                written += os.pwrite(fd, rest, offset + written)


def assemble(path, header, chunk_outputs, align):
    """Turn path, holding the chunk BAMs at the given (offset, end) byte ranges
    with reserved space between them, into a single BAM with the given
    uncompressed header. Gaps are removed with FALLOC_FL_COLLAPSE_RANGE where
    possible, which needs no data copying; otherwise the record blocks are
    moved down in place."""
    header_blocks = encode_header(header, align)
    bodies = []  # (start, end, last block offset, last block data)
    with open(path, 'rb') as fh:
        for offset, out_end in chunk_outputs:
            body_end = out_end - len(BGZF_EOF)
            fh.seek(body_end)
            if fh.read(len(BGZF_EOF)) != BGZF_EOF:
                raise ValueError('Chunk output does not end with a BGZF EOF block')
            body_start = header_end(fh, offset)
            if body_start is None:
                raise ValueError('Chunk header does not end on a BGZF block boundary')
            if body_end > body_start:
                bodies.append((body_start, body_end) +
                              tail_bgzf_blocks(fh, body_start, body_end))
    if not bodies:
        with open(path, 'r+b') as fh:
            fh.write(header_blocks + BGZF_EOF)
            fh.truncate()
        return

    # Byte ranges to keep, in order: the final header, the chunks' record
    # blocks, and the last chunk's EOF block
    if len(header_blocks) > bodies[0][0]:
        raise RuntimeError('Readgroup header exceeded its reserved space')
    segments = [[0, len(header_blocks)]] + \
        [[start, end] for start, end, _, _ in bodies]
    segments[-1][1] += len(BGZF_EOF)

    fd = os.open(path, os.O_RDWR)
    try:
        pwrite_all(fd, header_blocks, 0)
        os.ftruncate(fd, segments[-1][1])
        if all(start % align == 0 for start, _ in segments):
            try:
                gap = segments[1][0] - segments[0][1]
                collapse_range(fd, segments[0][1], gap)
                for segment in segments[1:]:
                    segment[0] -= gap
                    segment[1] -= gap
                segments[0][1] = segments[1][1]
                del segments[1]
                # segments[i] now holds chunk i's record blocks. Working back
                # from the end, re-encode each chunk's last record block so it
                # ends on a filesystem block boundary, then collapse the gap
                # after it.
                for i in range(len(segments) - 2, -1, -1):
                    block_offset = bodies[i][2] - gap
                    tail = stored_to_boundary(block_offset, bodies[i][3], align)
                    gap_start = block_offset + len(tail)
                    gap_end = segments[i + 1][0]
                    if gap_start >= gap_end:
                        break
                    pwrite_all(fd, tail, block_offset)
                    segments[i][1] = gap_start
                    collapse_range(fd, gap_start, gap_end - gap_start)
                    segments[i][1] += segments[i + 1][1] - segments[i + 1][0]
                    del segments[i + 1]
            except (OSError, AttributeError, ValueError):
                pass
        # Fall back to moving the remaining segments down in place
        dst = 0
        for start, end in segments:
            if start != dst:
                move_range(path, fd, start, end, dst)
            dst += end - start
        os.ftruncate(fd, dst)
    finally:
        os.close(fd)


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

    # Final BAM
    final_bam = in_bam_stem + '.readgroupsadded.bam'

    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as in_bam_fh:
        header_dict = in_bam_fh.header.to_dict()

    # The chunks are written straight into one sparse file, each at an offset
    # that leaves room for its output, then the gaps are removed. Finished
    # chunks are recorded in a checkpoint file, so a rerun in the same
    # directory, e.g. after preemption, only redoes unfinished chunks.
    temp_bam = final_bam + '.partial'
    checkpoint = final_bam + '.checkpoint'
    in_stat = Path(in_bam).stat()
    identity = dict(version=CHECKPOINT_VERSION, name=Path(in_bam).name,
                    size=in_stat.st_size, mtime_ns=in_stat.st_mtime_ns)
    resumed = load_checkpoint(checkpoint, temp_bam, identity)
    if resumed is None:
        # Drop any old checkpoint before temp_bam is truncated
        Path(checkpoint).unlink(missing_ok=True)
        layout = plan_layout(in_bam, temp_bam, args.proc, args.threads)
        layout['identity'] = identity
        done = {}
    else:
        layout, done = resumed
        print(f"Resuming: {len(done)} of {len(layout['chunks'])} chunks already done",
              file=sys.stderr)
    write_checkpoint(checkpoint, layout, done)
    chunks, offsets, limits = layout['chunks'], layout['offsets'], layout['limits']

    def finish(index, result):
        done[index] = result
        append_checkpoint(checkpoint, dict(chunk=index, readgroups=result[0],
                                           out_end=result[1]))

    jobs = {index: (in_bam, temp_bam, offsets[index], limits[index], start, end,
                    args.threads)
            for index, (start, end) in enumerate(chunks) if index not in done}
    if args.proc == 1:
        for index, job in jobs.items():
            finish(index, process_chunk(*job))
    elif jobs:
        with ProcessPoolExecutor(max_workers=min(args.proc, len(jobs))) as executor:
            futures = {executor.submit(process_chunk, *job): index
                       for index, job in jobs.items()}
            # Checkpoint chunks as they finish. On a failure, stop starting new
            # chunks, but still checkpoint the ones already running.
            error = None
            for future in as_completed(futures):
                if future.cancelled():
                    continue
                try:
                    finish(futures[future], future.result())
                except Exception as exc:
                    if error is None:
                        error = exc
                        for other in futures:
                            other.cancel()
            if error is not None:
                raise error
    results = [done[index] for index in range(len(chunks))]

    readgroups = {}
    for readgroup_ids, _ in results:
        for readgroup_id in readgroup_ids:
            readgroups.setdefault(readgroup_id, ReadGroup(readgroup_id))

    # Add readgroups to header
    header_dict['RG'] = [readgroup._asdict() for readgroup in readgroups.values()]
    header_bam = final_bam + '.header'
    with pysam.AlignmentFile(header_bam, 'wb', header=header_dict):
        pass
    header = bam_header_bytes(header_bam)
    Path(header_bam).unlink()

    # Assembly rewrites temp_bam in place, so a run interrupted from here on
    # can't be resumed
    append_checkpoint(checkpoint, dict(assembling=True))
    assemble(temp_bam, header,
             [(offset, out_end) for offset, (_, out_end) in zip(offsets, results)],
             layout['align'])
    os.replace(temp_bam, final_bam)
    Path(checkpoint).unlink()
