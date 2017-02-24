import numpy as np

####################################################################################################
# Exercise 1: DFT

def dft_matrix(n: int) -> np.ndarray:
    """
    Construct DFT matrix of size n.

    Arguments:
    n: size of DFT matrix

    Return:
    F: DFT matrix of size n

    Forbidden:
    - numpy.fft.*
    """
    # TODO: initialize matrix with proper size
    F = np.zeros((n, n), dtype='complex128')

    # TODO: create principal term for DFT matrix
    w = np.exp(-2* np.pi * 1j/n)
    # TODO: fill matrix with values
    for i in range (0,n):
        for j in range (0,n):
            F[i][j] = w**(i*j)
    return F


def is_unitary(matrix: np.ndarray) -> bool:
    """
    Check if the passed in matrix of size (n times n) is unitary up to a normalization factor of 1/n.

    Arguments:
    matrix: the matrix which is checked

    Return:
    unitary: True if the matrix is unitary
    """
    unitary = True
    
    # TODO: check that F is unitary, if not return false
    a, b = matrix.shape
    if a != b:
        unitary = False
        
    if np.allclose(matrix.dot(np.transpose(np.conj(matrix))),a*np.eye(a)) != True :
        unitary = False
        
    return unitary


def create_harmonics(n: int = 128) -> (list, list):
    """
    Create delta impulse signals and perform the fourier transform on each signal.

    Arguments:
    n: the length of each signal

    Return:
    sigs: list of np.ndarrays that store the delta impulse signals
    fsigs: list of np.ndarrays with the fourier transforms of the signals
    """

    # list to store input signals to DFT
    sigs = []
    # Fourier-transformed signals
    fsigs = []

    # TODO: create signals and extract harmonics out of DFT matrix
    F = dft_matrix(n)
    for i in range(0,n):
       signal = np.zeros(n)
       signal[i] = 1
       sigs.append(signal)
       fsigs.append(F.dot(signal))
    return sigs, fsigs


####################################################################################################
# Exercise 2: FFT

def shuffle_bit_reversed_order(data: np.ndarray) -> np.ndarray:
    """
    Shuffle elements of data using bit reversal of list index.

    Arguments:
    data: data to be transformed (shape=(n,), dtype='float64')

    Return:
    data: shuffled data array
    """

    # TODO: implement shuffling by reversing index bits
    size = data.size
    shuffleddata = np.zeros(size, dtype = 'complex128')
    for i in range(0, size):
        pos = int(bin(i)[2:].zfill(int(np.log2(size)))[::-1], 2)
        shuffleddata[pos] = data[i]
    
    return shuffleddata
    
    return data


def fft(data: np.ndarray) -> np.ndarray:
    """
    Perform real-valued discrete Fourier transform of data using fast Fourier transform.

    Arguments:
    data: data to be transformed (shape=(n,), dtype='float64')

    Return:
    fdata: Fourier transformed data

    Note:
    This is not an optimized implementation but one to demonstrate the essential ideas
    of the fast Fourier transform.

    Forbidden:
    - numpy.fft.*
    """

    fdata = np.asarray(data, dtype='complex128')
    n = fdata.size

    # check if input length is power of two
    if not n > 0 or (n & (n - 1)) != 0:
        raise ValueError

    # TODO: first step of FFT: shuffle data
    fdata = shuffle_bit_reversed_order(fdata)

    # TODO: second step, recursively merge transforms
    count = 1
    while n > count:
        for i in range(0, count):
            omega = np.complex(np.exp((-1j * 2.0 * i * np.pi)/(2*count)))
            for j in range(i,n,2*count):
                tmp = omega * fdata[j+count]
                fdata[j+count] = fdata[j] - tmp
                fdata[j] = fdata[j] + tmp
        count = count*2
    return fdata


def generate_tone(f: float = 261.626, num_samples: int = 44100) -> np.ndarray:
    """
    Generate tone with frequency f (default mid C: f = 261.626 Hz) and return the signal.

    Arguments:
    f: frequency of the tone

    Return:
    data: the generated signal
    """

    # sampling range
    x_min = 0.0
    x_max = 2.0 * np.pi
    data = np.zeros(num_samples)

    # TODO: Generate sine wave with proper frequency
    x = (x_max - x_min)/(num_samples - 1)
    for i in range(0, num_samples):
        data[i] = np.sin(f*x*i)

    return data


def low_pass_filter(adata: np.ndarray, bandlimit: int = 5000) -> np.ndarray:
    """
    Filter high frequencies above bandlimit.

    Arguments:
    adata: data to be filtered
    bandlimit: bandlimit above which to cut off frequencies

    Return:
    adata_filtered: filtered data
    """

    # TODO: compute Fourier transform of input data
    adata = fft(adata)
    # TODO: set high frequencies above bandlimit to zero, make sure the almost symmetry of the transform is respected
    adata[bandlimit+1:adata.size-bandlimit] = 0
    # TODO: compute inverse transform and extract real component
    var = np.conjugate(fft(np.conjugate(adata)))
    adata_filtered = np.real(1/adata.size * var)
    
    return adata_filtered


if __name__ == '__main__':
    print("All requested functions for the assignment have to be implemented in this file and uploaded to the "
          "server for the grading.\nTo test your implemented functions you can "
          "implement/run tests in the file tests.py (> python3 -v test.py [Tests.<test_function>]).")
