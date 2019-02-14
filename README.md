# Dynamic Time Warping
This repository implements dynamic time warping (or equivalently Levenshtein distance) between two multi-dimensional sequences and returns the time warped path. The weight for moving along diagonal path can be specified. Also a band width value can be specified so that the path is only looked for within that width along the diagonal. For the path, it returns the distances along the path as well as corresponding aligned indices of both the sequences. The path does not align the entire sequences, it starts from the end of both sequences and goes till one of the sequence reaches the beginning or get used up entirely (more probably the shorter sequence). The total DTW cost can be found by summing up the individual path distances.

## Python bindings
There are also python bindings for the compiled library using ctypes. The module PyDTW can be imported as it and used. As of now, the precompiled library is only for OS X, you can however compile it for linux and place the library inside the module to make it work on Linux.

Usage
```python
import numpy as np
from PyDTW import DTW_Extractor
e = DTW_Extractor()
seq1 = np.random.random([100, 20]) # sequence 1 is of length 100 and dimension 20
seq2 = np.random.random([130, 20]) # sequence 2 is of length 100 and dimension 20

distances, indices1, indices2 = e.compute_dtw(seq1, seq2, w_diag=1.1, band_win=5)
# the cost of moving along diagonal is 1.1 times moving up or to side.
# the path is found in window of [diagonal - 5, diagonal + 5]
```

#### TODO:
- Create a makefile
- Compile and add library for linux
