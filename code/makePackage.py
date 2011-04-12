import os
import string

def constructArborescense(path):
    """this function construct the arborescence of a directory, i.e. it seeks for all its sub-directories"""
    arborescence = []
    for root, direct, files in os.walk(path):
        if '.git' not in root and 'unused' not in root:
            if files != []:
                arborescence.append(root)
    return arborescence

def copySources(SHORT_LICENCE, COPYRIGHT, CONTACT, pathToCopyFrom, pathToCopyTo):
    short_licence = SHORT_LICENCE

    listSourcesToCopy = os.listdir(pathToCopyFrom)
    for sourceFileName in listSourcesToCopy:
        if not os.path.isdir(os.path.join(pathToCopyFrom, sourceFileName)):
            f = open(os.path.join(pathToCopyFrom, sourceFileName), 'r')
            sourceFile = f.readlines()
            f.close()
            newSourceFile = []
            COPY = False
            if '.cpp' in sourceFileName or '.h' in sourceFileName or '.geo' in sourceFileName:
                newSourceFile.append('/**********************************************************************\n')
                newSourceFile.append(' *\n')
                newSourceFile.append(' * ')
                newSourceFile.append(sourceFileName)
                newSourceFile.append('\n')
                newSourceFile.append(' *\n')
                newSourceFile.append(' * ')
                newSourceFile.append(COPYRIGHT)
                newSourceFile.append('\n')
                newSourceFile.append(' *\n')
                for line in short_licence:
                    newLine = ' * ' + line
                    newSourceFile.append(newLine)
                newSourceFile.append(' *\n')
                newSourceFile.append(' * ')
                newSourceFile.append(CONTACT)
                newSourceFile.append('\n')
                newSourceFile.append(' *\n')
                newSourceFile.append(' **********************************************************************/\n')
                newSourceFile.append('\n')

                COPY = True

            elif '.py' in sourceFileName and 'makePackage' not in sourceFileName or 'makefile' in sourceFileName:
                newSourceFile.append('#######################################################################\n')
                newSourceFile.append('##\n')
                newSourceFile.append('## ')
                newSourceFile.append(sourceFileName)
                newSourceFile.append('\n')
                newSourceFile.append('##\n')
                newSourceFile.append('## ')
                newSourceFile.append(COPYRIGHT)
                newSourceFile.append('\n')
                newSourceFile.append('##\n')
                for line in short_licence:
                    newLine = '## ' + line
                    newSourceFile.append(newLine)
                newSourceFile.append('##\n')
                newSourceFile.append('## ')
                newSourceFile.append(CONTACT)
                newSourceFile.append('\n')
                newSourceFile.append('##\n')
                newSourceFile.append('#######################################################################\n')
                newSourceFile.append('\n')

                COPY = True

            elif '.f' in sourceFileName or '.sh' in sourceFileName or 'COPYING' in sourceFileName or 'REFERENCES' in sourceFileName or 'GUIDE' in sourceFileName or 'hostfile' in sourceFileName or 'repo' in sourceFileName:
                COPY = True

            elif 'Doxyfile' in sourceFileName or 'README' in sourceFileName or 'EXAMPLES' in sourceFileName:
                COPY = True

            else:
                COPY = False

            if COPY:
                for line in sourceFile:
                    newSourceFile.append(line)
                f = open(os.path.join(pathToCopyTo, sourceFileName), 'w')
                for line in newSourceFile:
                    f.write(line)
                f.close()
                if '.sh' in sourceFileName:
                    os.system("chmod ug+rwx " + os.path.join(pathToCopyTo, sourceFileName))

if __name__=="__main__":
    RELEASE = '0.0.0'
    MAIN_DIR = '.'
    SHORT_LICENCE = ['This file is part of Puma-EM.\n', '\n', 'Puma-EM is free software: you can redistribute it and/or modify\n', 'it under the terms of the GNU General Public License as published by\n', 'the Free Software Foundation, either version 3 of the License, or\n', '(at your option) any later version.\n', '\n', 'Puma-EM is distributed in the hope that it will be useful,\n', 'but WITHOUT ANY WARRANTY; without even the implied warranty of\n', 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n', 'GNU General Public License for more details.\n', '\n', 'You should have received a copy of the GNU General Public License\n', 'along with Puma-EM.  If not, see <http://www.gnu.org/licenses/>.\n']
    COPYRIGHT = "Copyright (C) 2011 Idesbald Van den Bosch"
    CONTACT = "Suggestions/bugs : <vandenbosch.idesbald@gmail.com>"

    arborescence = constructArborescense(MAIN_DIR)
    print arborescence

    PACKAGE_DIR = os.path.join(MAIN_DIR, 'Puma-EM')
    for pathToCopyFrom in arborescence:
        if pathToCopyFrom=='.':
            pathToCopyTo = PACKAGE_DIR
            os.mkdir(pathToCopyTo)
        else:
            pathToCopyFrom = pathToCopyFrom[2:]
            pathToCopyTo = os.path.join(PACKAGE_DIR, pathToCopyFrom)
            os.makedirs(pathToCopyTo)
        copySources(SHORT_LICENCE, COPYRIGHT, CONTACT, pathToCopyFrom, pathToCopyTo)

