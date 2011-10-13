import os, sys
from scipy import zeros, real, imag, resize
from scipy import weave
from scipy.weave import converters

def readIntFromDisk(filename):
    f = open(filename, 'r')
    out = f.readlines();
    f.close()
    return int(out[0])

def readFloatFromDisk(filename):
    f = open(filename, 'r')
    out = f.readlines();
    f.close()
    return float(out[0])

def writeScalarToDisk(x, filename):
    f = open(filename, 'w')
    if isinstance(x, complex):
        string = '(' + str(real(x)) + ',' + str(imag(x)) + ')'
        f.write(string)
    else:
        f.write(str(x))
    f.close()

def writeASCIIBlitzArrayToDisk(A, filename):
    wrapping_code = """
    using namespace blitz;
    blitz::ofstream ofs(filename.c_str());
    if (!ofs.is_open())
    {
        cerr << "ReadWriteBlitzArray.py::writeASCIIBlitzArrayToDisk: Unable to write to file: " << filename << endl;
        exit(1);
    }
    ofs.precision(18);
    ofs << A << endl;
    ofs.close();
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])

def readASCIIBlitzIntArray2DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    f.close()
    # we need to take into account the new Blitz++ ASCII format
    # for expressing an array's dimensions
    indLines_expr, indColumns_expr = data.split('x')[0], data.split('x')[1]
    if '(' not in indLines_expr and ',' not in indLines_expr and ')' not in indLines_expr:
        Nl, Nc = int(indLines_expr), int(indColumns_expr)
    else:
        indLines_expr = indLines_expr.replace('(','')
        indLines_expr = indLines_expr.replace(')','')
        indLines_expr = indLines_expr.replace('\n','')
        indColumns_expr = indColumns_expr.replace('(','')
        indColumns_expr = indColumns_expr.replace(')','')
        indColumns_expr = indColumns_expr.replace('\n','')
        startIndLines, endIndLines = int(indLines_expr.split(',')[0]), int(indLines_expr.split(',')[1])
        startIndColumns, endIndColumns = int(indColumns_expr.split(',')[0]), int(indColumns_expr.split(',')[1])
        Nl = endIndLines - startIndLines + 1
        Nc = endIndColumns - startIndColumns + 1
    A = zeros((Nl, Nc), 'i')
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str());
    if (!ifs.is_open())
    {
        cerr << "ReadWriteBlitzArray.py::readASCIIBlitzIntArray2DFromDisk: Unable to read from file: " << filename << endl;
        exit(1);
    }
    ifs.precision(18);
    blitz::Array<int, 2> B;
    ifs >> B;
    ifs.close();
    A = B;
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return A

def readASCIIBlitzComplexFloatArray2DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    f.close()
    # we need to take into account the new Blitz++ ASCII format
    # for expressing an array's dimensions
    indLines_expr, indColumns_expr = data.split('x')[0], data.split('x')[1]
    if '(' not in indLines_expr and ',' not in indLines_expr and ')' not in indLines_expr:
        Nl, Nc = int(indLines_expr), int(indColumns_expr)
    else:
        indLines_expr = indLines_expr.replace('(','')
        indLines_expr = indLines_expr.replace(')','')
        indLines_expr = indLines_expr.replace('\n','')
        indColumns_expr = indColumns_expr.replace('(','')
        indColumns_expr = indColumns_expr.replace(')','')
        indColumns_expr = indColumns_expr.replace('\n','')
        startIndLines, endIndLines = int(indLines_expr.split(',')[0]), int(indLines_expr.split(',')[1])
        startIndColumns, endIndColumns = int(indColumns_expr.split(',')[0]), int(indColumns_expr.split(',')[1])
        Nl = endIndLines - startIndLines + 1
        Nc = endIndColumns - startIndColumns + 1
    A = zeros((Nl, Nc), 'F')
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str());
    if (!ifs.is_open())
    {
        cerr << "ReadWriteBlitzArray.py::readASCIIBlitzComplexFloatArray2DFromDisk: Unable to read from file: " << filename << endl;
        exit(1);
    }
    ifs.precision(18);
    blitz::Array<std::complex<float>, 2> B;
    ifs >> B;
    ifs.close();
    A = B;
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])

    return A

def readASCIIBlitzFloatArray2DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    f.close()
    # we need to take into account the new Blitz++ ASCII format
    # for expressing an array's dimensions
    indLines_expr, indColumns_expr = data.split('x')[0], data.split('x')[1]
    if '(' not in indLines_expr and ',' not in indLines_expr and ')' not in indLines_expr:
        Nl, Nc = int(indLines_expr), int(indColumns_expr)
    else:
        indLines_expr = indLines_expr.replace('(','')
        indLines_expr = indLines_expr.replace(')','')
        indLines_expr = indLines_expr.replace('\n','')
        indColumns_expr = indColumns_expr.replace('(','')
        indColumns_expr = indColumns_expr.replace(')','')
        indColumns_expr = indColumns_expr.replace('\n','')
        startIndLines, endIndLines = int(indLines_expr.split(',')[0]), int(indLines_expr.split(',')[1])
        startIndColumns, endIndColumns = int(indColumns_expr.split(',')[0]), int(indColumns_expr.split(',')[1])
        Nl = endIndLines - startIndLines + 1
        Nc = endIndColumns - startIndColumns + 1
    A = zeros((Nl, Nc), 'f')
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str());
    if (!ifs.is_open())
    {
        cerr << "ReadWriteBlitzArray.py::readASCIIBlitzFloatArray2DFromDisk: Unable to read from file: " << filename << endl;
        exit(1);
    }
    ifs.precision(18);
    blitz::Array<float, 2> B;
    ifs >> B;
    ifs.close();
    A = B;
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return A

def readASCIIBlitzFloatArray1DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    f.close()
    # we need to take into account the new Blitz++ ASCII format
    # for expressing an array's dimensions
    indLines_expr = data.split()[0]
    if '(' not in indLines_expr and ',' not in indLines_expr and ')' not in indLines_expr:
        Nl = int(indLines_expr)
    else:
        indLines_expr = indLines_expr.replace('(','')
        indLines_expr = indLines_expr.replace(')','')
        indLines_expr = indLines_expr.replace('\n','')
        startIndLines, endIndLines = int(indLines_expr.split(',')[0]), int(indLines_expr.split(',')[1])
        Nl = endIndLines - startIndLines + 1
    A = zeros(Nl, 'f')
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str());
    if (!ifs.is_open())
    {
        cerr << "ReadWriteBlitzArray.py::readASCIIBlitzFloatArray1DFromDisk: Unable to read from file: " << filename << endl;
        exit(1);
    }
    ifs.precision(18);
    blitz::Array<float, 1> B;
    ifs >> B;
    ifs.close();
    A = B;
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return A

def readASCIIBlitzIntArray1DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    f.close()
    # we need to take into account the new Blitz++ ASCII format
    # for expressing an array's dimensions
    indLines_expr = data.split()[0]
    if '(' not in indLines_expr and ',' not in indLines_expr and ')' not in indLines_expr:
        Nl = int(indLines_expr)
    else:
        indLines_expr = indLines_expr.replace('(','')
        indLines_expr = indLines_expr.replace(')','')
        indLines_expr = indLines_expr.replace('\n','')
        startIndLines, endIndLines = int(indLines_expr.split(',')[0]), int(indLines_expr.split(',')[1])
        Nl = endIndLines - startIndLines + 1
    A = zeros(Nl, 'i')
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str());
    if (!ifs.is_open())
    {
        cerr << "ReadWriteBlitzArray.py::readASCIIBlitzFloatArray1DFromDisk: Unable to read from file: " << filename << endl;
        exit(1);
    }
    blitz::Array<int, 1> B;
    ifs >> B;
    ifs.close();
    A = B;
    """
    weave.inline(wrapping_code,
                 ['A', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w', '-fPIC'])
    return A

def writeBlitzArrayToDisk(A, filename):
    sizeOfItem = A.itemsize
    wrapping_code = """
    using namespace blitz;
    blitz::ofstream fout(filename.c_str(), blitz::ios::binary);
    fout.write((char *)(A.data()), A.size()*sizeOfItem);
    fout.close();
    """
    weave.inline(wrapping_code,
                 ['A', 'sizeOfItem', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])

def readBlitzArrayFromDisk(filename, Nl, Nc, ELEM_TYPE):
    A = zeros((Nl, Nc), ELEM_TYPE)
    sizeOfItem = A.itemsize
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str(), blitz::ios::binary);
    if (! ifs.is_open()) { 
      cout << "ReadWriteBlitzArray.py::readBlitzArrayFromDisk: error opening " << filename << endl;
      exit (1);
    }
    ifs.read((char *)(A.data()), A.size()*sizeOfItem);
    ifs.close();
    """
    weave.inline(wrapping_code,
                 ['A', 'sizeOfItem', 'filename'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>','<fstream>','<sstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return A

def read1DBlitzArrayFromDisk(filename, ELEM_TYPE):
    B = zeros(1, ELEM_TYPE)
    sizeOfItem = B.itemsize
    sizeOfA = zeros(1, 'i')
    wrapping_code = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str(), blitz::ios::binary);
    if (! ifs.is_open()) { 
      cout << "ReadWriteBlitzArray.py::read1DBlitzArrayFromDisk: error opening " << filename << endl;
      exit (1);
    }
    ifs.seekg (0, blitz::ios::end);
    int length = ifs.tellg();
    ifs.close();
    sizeOfA(0) = length/sizeOfItem;
    """
    weave.inline(wrapping_code,
                 ['sizeOfA','sizeOfItem', 'filename'],
                 type_converters = converters.blitz,
                 headers = ['<iostream>','<fstream>','<sstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    # resizing A
    A = zeros(sizeOfA[0], ELEM_TYPE)
    wrapping_code2 = """
    using namespace blitz;
    blitz::ifstream ifs(filename.c_str(), blitz::ios::binary);
    ifs.read((char *)(A.data()), A.size()*sizeOfItem);
    ifs.close();
    """
    weave.inline(wrapping_code2,
                 ['A', 'sizeOfItem', 'filename'],
                 type_converters = converters.blitz,
                 headers = ['<iostream>','<fstream>','<sstream>','<complex>','<string>','<blitz/array.h>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return A

def writeToDisk_chunk_of_Z_sparse(path, name, Z, q_array, rowIndexToColumnIndexes, RWG_numbers, chunkNumber):
    """this function writes to disk the chunks of Z sparse and the corresponding indexes arrays, each with a number"""
    chunkNumberString = str(chunkNumber)
    writeBlitzArrayToDisk(Z, os.path.join(path, name) + str(chunkNumber) + '.txt')
    writeBlitzArrayToDisk(q_array, os.path.join(path, 'q_array') + str(chunkNumber) + '.txt')
    writeBlitzArrayToDisk(rowIndexToColumnIndexes, os.path.join(path, 'rowIndexToColumnIndexes') + str(chunkNumber) + '.txt')
    writeBlitzArrayToDisk(RWG_numbers, os.path.join(path, 'RWG_numbers') + str(chunkNumber) + '.txt')
    # now we write the scalar values
    N_RWG_File = os.path.join(path, 'N_RWG') + str(chunkNumber) + '.txt'
    writeScalarToDisk(RWG_numbers.shape[0], N_RWG_File)
    N_near_File = os.path.join(path, 'N_near') + str(chunkNumber) + '.txt'
    writeScalarToDisk(Z.shape[0], N_near_File)
    N_q_array_File = os.path.join(path, 'N_q_array') + str(chunkNumber) + '.txt'
    writeScalarToDisk(q_array.shape[0], N_q_array_File)


if __name__=="__main__":
    x = 4.0 + 1.j * 3.3
    filename = './tmp/x.txt'
    writeScalarToDisk(x, filename)

    
    
