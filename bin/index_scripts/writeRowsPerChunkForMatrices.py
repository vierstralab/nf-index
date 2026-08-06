# Port of cpp to python

import sys
import numpy as np

def load_order_map(path, what):
    m = {}
    with open(path) as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if "\t" in line:
                sys.exit(f"Error: Tab-delimited data on line {i} of {path}.\n"
                         f"Contents should be one column of IDs defining output {what} order.")
            m[line] = i
    return m

def main():
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} sampleIDs chunkIDorder chunkData8col outBinary")

    samp_to_col  = load_order_map(sys.argv[1], "column")
    chunk_to_row = load_order_map(sys.argv[2], "row")
    ncols = len(samp_to_col)

    ro, co = [], []
    seen_rows = set()
    with open(sys.argv[3]) as f:
        for linenum, line in enumerate(f, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                sys.exit(f"Error: line {linenum} of {sys.argv[3]} has {len(fields)} fields; expected 8.")
            samp, chunk = fields[3], fields[7]
            r = chunk_to_row.get(chunk)
            if r is None:
                sys.exit(f"Error: Chunk ID {chunk} (line {linenum}) not in ordered chunkIDs.")
            c = samp_to_col.get(samp)
            if c is None:
                sys.exit(f'Error: Sample "{samp}" (line {linenum}) not in sample IDs.')
            ro.append(r); co.append(c - 1)
            seen_rows.add(r)

    sorted_rows = sorted(seen_rows)
    row_to_out = {r: i for i, r in enumerate(sorted_rows)}
    nrows = len(sorted_rows)

    ro = np.fromiter((row_to_out[r] for r in ro), dtype=np.int64, count=len(ro))
    co = np.asarray(co, dtype=np.int64)

    B = np.zeros((nrows, ncols), dtype=bool)
    B[ro, co] = True

    np.savetxt(sys.argv[4], B.astype(np.int8), fmt="%d", delimiter="\t")

if __name__ == "__main__":
    main()