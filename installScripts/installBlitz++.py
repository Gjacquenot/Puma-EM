import os

import os

def installBlitz():
    URL = "http://ftp.debian.org/debian/pool/main/b/blitz++/"
    TARGET = "blitz++_0.9.orig"
    os.system("wget " + URL + TARGET + ".tar.gz")
    os.system("tar xvzf " + TARGET + ".tar.gz")
    os.chdir("blitz++-0.9.orig")
    os.system("./configure")
    os.system("make lib")
    os.system("make install")

if __name__=="__main__":
    installBlitz()



