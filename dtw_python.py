import ctypes
import numpy as np
import platform
import librosa

class DTW_Extractor(object):

    class DTW_path(ctypes.Structure):
        """
        Struct for DTW output. Contains DTW path distances as float array and the length of the path
        """
        _fields_ = [('distances', ctypes.POINTER(ctypes.c_float)),
                    ('length', ctypes.c_int)]

    def __init__(self, config_path):
        super(DTW_Extractor, self).__init__()
        if platform.system() == 'Linux':
            self.libmfcc = ctypes.CDLL('./libMfccExtractor.so')
        elif platform.system() == 'Darwin':
            self.libmfcc = ctypes.CDLL('./libDTW.dylib')
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
        self.deallocate_dtw_func.argtypes = [DTW_Extractor.DTW_path]


    def compute_dtw(self, seq1, seq2, w_diag, band_win):
        """
        Compute DTW givene the audio data of one utterance (usually recorded utterance)
        and mfccs of the second utterance
        :param audio_data: np.ndarray, one dimensional audio signal
        :param mfccs: np.ndarray, two dimensional mfcc data
        :return: np.ndarray, distances of the DTW path
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
        self.deallocate_dtw_func(result)
        return dtw_dist


if __name__ == '__main__':
    """Some testing code"""
    ext = DTW_Extractor('mfcc_config')
    t, _ = librosa.core.load('orig1.wav', sr=None)
    a = t[:1000].reshape((-1,20))
    b = t[1000:2000].reshape((-1,20))
    dtw = ext.compute_dtw(a, b, 1., -1)
    print(dtw)
