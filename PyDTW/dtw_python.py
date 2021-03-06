"""
dtw_python.cpp
libDTW

Author: Gagandeep Singh, 13-02-2019
"""

import ctypes
import numpy as np
import platform
import librosa
import os

class DTW_Extractor(object):

    class DTW_path(ctypes.Structure):
        """
        Struct for DTW output. Contains DTW path distances as float array and the length of the path
        """
        _fields_ = [('distances', ctypes.POINTER(ctypes.c_float)),
                    ('seq1_indices', ctypes.POINTER(ctypes.c_int)),
                    ('seq2_indices', ctypes.POINTER(ctypes.c_int)),
                    ('length', ctypes.c_int)]

    def __init__(self):
        super(DTW_Extractor, self).__init__()
        if platform.system() == 'Linux':
            self.libmfcc = ctypes.CDLL(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'libDTW.so'))
        elif platform.system() == 'Darwin':
            self.libmfcc = ctypes.CDLL(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'libDTW.dylib'))
        else:
            raise ValueError('Platform: %s not supported'%(platform.system()))

        self.compute_dtw_func = self.libmfcc.dtw_path
        self.compute_dtw_func.argtypes = [ctypes.POINTER(ctypes.c_float), 
                                          ctypes.POINTER(ctypes.c_float),
                                          ctypes.c_int, 
                                          ctypes.c_int,
                                          ctypes.c_int,
                                          ctypes.c_float, 
                                          ctypes.c_int]
        self.compute_dtw_func.restype = DTW_Extractor.DTW_path

        self.deallocate_dtw_func = self.libmfcc.deallocate_DTW
        self.deallocate_dtw_func.argtypes = [ctypes.c_void_p]


    def compute_dtw(self, seq1, seq2, w_diag, band_win):
        """
        Compute DTW given two sequences
        :param seq1: np.ndarray, [len_seq1, dim] first sequence
        :param seq2: np.ndarray, [len_seq2, dim] second sequence
        :param w_diag: float, weight on moving diagonally
        :param band_win: int, find path within band_win (on both sides) of the diagonal. -1 will look unrestricted
        :return: (np.ndarray,)*3, [dtw_path_len,] distances of the DTW path, indices of seq1, indices of seq2
        """
        seq1_len = seq1.shape[0]
        dim1 = seq1.shape[1]
        seq1 = seq1.flatten().astype(dtype='float32')
        seq1_p = np.ctypeslib.as_ctypes(seq1)

        seq2_len = seq2.shape[0]
        dim2 = seq2.shape[1]
        assert dim1 == dim2, 'dimension of the two sequences does not match'
        seq2 = seq2.flatten().astype(dtype='float32')
        seq2_p = np.ctypeslib.as_ctypes(seq2)
        
        result = self.compute_dtw_func(seq1_p, seq2_p, seq1_len, seq2_len, dim1, w_diag, band_win)
        temp = np.ctypeslib.as_array(result.distances, (result.length,))
        dtw_dist = np.copy(temp)
        temp = np.ctypeslib.as_array(result.seq1_indices, (result.length,))
        seq1_indices = np.copy(temp)
        temp = np.ctypeslib.as_array(result.seq2_indices, (result.length,))
        seq2_indices = np.copy(temp)
        self.deallocate_dtw_func(result.distances)
        self.deallocate_dtw_func(result.seq1_indices)
        self.deallocate_dtw_func(result.seq2_indices)
        return dtw_dist, seq1_indices, seq2_indices


if __name__ == '__main__':
    """Some testing code"""
    ext = DTW_Extractor()
    t = np.random.random((1500))
    a = t[:1000].reshape((-1,20))
    b = t[40:1040].reshape((-1,20))
    distances, indices_a, indices_b = ext.compute_dtw(a, b, 1., -1)
