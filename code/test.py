from scipy import weave, rand, product
from scipy.weave import converters
from scipy import zeros, array, arange
import time
from mpi4py import MPI

def writeBinaryArray(A, filename):
    """write a binary array in a file"""
    wrapping_code = """
    using namespace blitz;
    cout << sizeof(A) << endl;
    blitz::ofstream fout(filename.c_str(), blitz::ios::binary);
    fout.write((char *)(A.data()), A.size()*8);
    cout << sum(A) << endl;
    fout.close();
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<blitz/array.h>'],
                 compiler = 'gcc')
    return

def readBinaryArray(N, filename):
    """read a binary array in a file"""
    wrapping_code = """
    using namespace blitz;
    Array<complex<float>, 2> B(N, 2);
    blitz::ifstream fin(filename.c_str(), blitz::ios::binary);
    fin.read((char *)(B.data()), B.size()*8);
    cout << sum(B) << endl;
    fin.close();
    """
    weave.inline(wrapping_code,
                 ['N', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<blitz/array.h>'],
                 compiler = 'gcc')

def createZ_sparse_MLFMA(chunkNumber):
    wrapping_code = """
    Z_sparse_MLFMA Z_near;
    Z_near.setZ_sparse_MLFMAFromFile(chunkNumber);
    """
    weave.inline(wrapping_code,
                 ['chunkNumber'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<blitz/array.h>','"Z_sparse_MLFMA.h"'],
                 compiler = 'gcc')


if __name__=="__main__":
    #N = int(4)
    #A = (rand(N, 2) + 1.j * rand(N, 2)).astype('F')
    #print product(A.shape)*A.itemsize/1024.**2
    #t0 = time.time()
    #filename = "A.txt"
    ##writeBinaryArray(A, filename)
    #time_writeBinaryArray = time.time() - t0
    #t0 = time.time()
    ##readBinaryArray(N, filename)
    #time_readBinaryArray = time.time() - t0
    #print "time for writing Binary Array =", time_writeBinaryArray
    #print "time for reading Binary Array =", time_readBinaryArray

    #createZ_sparse_MLFMA(0)
    sendbuf = rand(3,3).tolist()
    recvbuf = zeros(3).tolist()
    #MPI.Init()
    status = MPI.Status()
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    print "           my ID is", my_id, "in", num_proc, "processes"
    if (my_id == 1):
        tag = 0
        print "sendbuf = ", sendbuf
        MPI.COMM_WORLD.Send([sendbuf, MPI.FLOAT], 0, tag)
    else:
        tag = my_id
        recvbuf = MPI.COMM_WORLD.Recv(1, tag, status)
        print "recvbuf = ", recvbuf
    MPI.Finalize()


