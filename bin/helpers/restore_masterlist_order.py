import sys
import numpy as np

np_array = np.load(sys.argv[1])
orig_idx = np.atleast_1d(np.loadtxt(sys.argv[2], dtype=int))
n_dhs = int(sys.argv[3])

restored = np.zeros(n_dhs, dtype=np_array.dtype)
restored[orig_idx] = np_array
np.save(sys.argv[4], restored)
