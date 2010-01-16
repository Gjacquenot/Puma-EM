import os

def installGMSH():
    URL = "http://www.geuz.org/gmsh/bin/Linux/"
    TARGET = "gmsh-2.4.2"
    os.system("wget " + URL + TARGET + "-Linux.tgz")
    os.system("tar xzf " + TARGET + "-Linux.tgz")
    os.system("mv " + os.path.join(TARGET, "gmsh") + " /usr/bin")    
    os.system("rm -rf " + TARGET + "-Linux.tgz")

if __name__=="__main__":
    installGMSH()


