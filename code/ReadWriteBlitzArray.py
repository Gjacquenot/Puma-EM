import os, sys
import numpy
from scipy import zeros, real, imag, array

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
    f = open(filename, 'w')
    dimensions = A.shape 
    N_dim = len(dimensions)
    string = '(0,' + str(dimensions[0]-1) + ')'
    for i in dimensions[1:]:
        string += ' x (0,' + str(i-1) + ')'
    string += '\n'
    f.write(string)
    f.write('[ ')
    if 'complex' not in A.dtype.name:
        if N_dim==1:  # 1-D arrays
            for i in range(dimensions[0]):
                string = str(A[i])
                f.write(string + ' ')
            f.write(']\n')

        elif N_dim==2: # 2-D arrays
            for i in range(dimensions[0]):
                for j in range(dimensions[1]):
                    string = str(A[i, j])
                    f.write(string + ' ')
                if i<dimensions[0]-1:
                    f.write('\n')
            f.write(']\n')

    else: # array elements are complex
        if N_dim==1: # 1-D arrays
            for i in range(dimensions[0]):
                string1 = str(real(A[i]))
                string2 = str(imag(A[i]))
                string = '(' + string1 + ',' + string2 + ') '
                f.write(string)
            f.write(']\n')

        elif N_dim==2: # 2-D arrays
            for i in range(dimensions[0]):
                for j in range(dimensions[1]):
                    string1 = str(real(A[i,j]))
                    string2 = str(imag(A[i,j])) 
                    string = '(' + string1 + ',' + string2 + ') '
                    f.write(string)
                if i<dimensions[0]-1:
                    f.write('\n')
            f.write(']\n')
    f.close()


def get2DArrayDimensions(data):
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
    return Nl, Nc

def readASCIIBlitzIntArray2DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    Nl, Nc = get2DArrayDimensions(data)
    A = zeros((Nl, Nc), 'i')
    if Nl==1:
        data = f.readline()
        A[0] = list(map(int, data.split()[1:-1]))
    else:
        for i in range(Nl):
            data = f.readline()
            if i==0:
                # we need to jump over the '[' character
                A[0] = list(map(int, data.split()[1:]))
            elif i<Nl-1:
                A[i] = list(map(int, data.split()))
            else:
                # we need to ignore the last ']' character
                A[Nl-1] = list(map(int, data.split()[:-1]))
    f.close()
    return A

def readASCIIBlitzComplexFloatArray2DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    Nl, Nc = get2DArrayDimensions(data)
    A = zeros((Nl, Nc), 'F')
    if Nl==1:
        data = f.readline()
        line = data[1:-1].split()
        for j in range(Nc):
            temp2 = line[j][1:-1].split(',')
            A[0, j] = float(temp2[0]) + 1.j * float(temp2[1])
    else:
        for i in range(Nl):
            data = f.readline()
            if i==0:
                # we need to jump over the '[' character
                line = data[1:].split()
                for j in range(Nc):
                    temp2 = line[j][1:-1].split(',')
                    A[0, j] = float(temp2[0]) + 1.j * float(temp2[1])
            elif i<Nl-1:
                line = data.split()
                for j in range(Nc):
                    temp2 = line[j][1:-1].split(',')
                    A[i, j] = float(temp2[0]) + 1.j * float(temp2[1])
            else:
                # we need to ignore the last ']' character
                line = data[:-1].split()
                for j in range(Nc):
                    temp2 = line[j][1:-1].split(',')
                    A[i, j] = float(temp2[0]) + 1.j * float(temp2[1])
    f.close()
    return A

def readASCIIBlitzFloatArray2DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline()
    Nl, Nc = get2DArrayDimensions(data)
    A = zeros((Nl, Nc), 'f')
    if Nl==1:
        data = f.readline()
        A[0] = array(list(map(float, data.split()[1:-1])))
    else:
        for i in range(Nl):
            data = f.readline()
            if i==0:
                # we need to jump over the '[' character
                A[0] = list(map(float, data.split()[1:]))
            elif i<Nl-1:
                A[i] = list(map(float, data.split()))
            else:
                A[Nl-1] = list(map(float, data.split()[:-1]))
    f.close()
    return A

def readASCIIBlitzFloatArray1DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline() # we skip the dimensions
    data = f.readline()
    A = array(list(map(float, data.split()[1:-1])))
    f.close()
    return A

def readASCIIBlitzIntArray1DFromDisk(filename):
    f = open(filename, 'r')
    data = f.readline() # we skip the dimensions
    data = f.readline()
    A = array(list(map(int, data.split()[1:-1])))
    f.close()
    return A

def writeBlitzArrayToDisk(A, filename):
    output_file = open(filename, 'wb')
    A.tofile(output_file)
    output_file.close()

def readBlitzArrayFromDisk(filename, Nl, Nc, ELEM_TYPE):
    input_file = open(filename, 'rb')
    A = numpy.fromfile(file=input_file, dtype=ELEM_TYPE).reshape((Nl, Nc))
    input_file.close()
    return A

def read1DBlitzArrayFromDisk(filename, ELEM_TYPE):
    input_file = open(filename, 'rb')
    A = numpy.fromfile(file=input_file, dtype=ELEM_TYPE)
    input_file.close()
    return A


if __name__=="__main__":
    x = 4.0 + 1.j * 3.3
#    filename = './tmp/x.txt'
#    writeScalarToDisk(x, filename)



