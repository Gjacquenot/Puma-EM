import os

def installGMSH():
    URL = "http://www.geuz.org/gmsh/bin/Linux/"
    TARGET = "gmsh-2.5.0-Linux"
    os.system("wget " + URL + TARGET + ".tgz")
    os.system("tar xzf " + TARGET + ".tgz")
    os.system("mv " + os.path.join(TARGET, "bin/gmsh") + " /usr/bin")    
    os.system("rm -rf " + TARGET + ".tgz")

if __name__=="__main__":
    installGMSH()


