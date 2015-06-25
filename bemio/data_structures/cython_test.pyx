import cython
cimport cython

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef _process(np.ndarray[DTYPE_t, ndim=2] array):

    cdef unsigned int rows = array.shape[0]
    cdef unsigned int cols = array.shape[1]
    cdef unsigned int row
    cdef np.ndarray[DTYPE_t, ndim=2] out = np.zeros((rows, cols))

    for row in range(0, rows):
        for col in range(0, cols):
            for row2 in range(0, rows):
                out[row, col] += array[row2, col] - array[row, col]

    return out

def main():
    cdef np.ndarray[DTYPE_t, ndim=2] data
    cdef np.ndarray[DTYPE_t, ndim=2] out
    data = np.load('data.npy')
    out = _process(data)
    np.save('vialoop.npy', out)