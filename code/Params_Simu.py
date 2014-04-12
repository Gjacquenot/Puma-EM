import string

class Params_Simu:

    def display(self):
        """
        Function that displays the values of the attributes of the class.
        @return: No returned value
        """

        attributes = dir(self)
        MAX_LENGTH = 0
        for att in list(range(0,len(attributes))):
            this_at = attributes[att]
            if this_at[0] != '_' and this_at != 'create_attribute' and this_at != 'display' and this_at != 'saveTo':
                if len(this_at) > MAX_LENGTH:
                    MAX_LENGTH = len(this_at)
        for att in list(range(0,len(attributes))):
            this_at = attributes[att]
            if this_at[0] != '_' and this_at != 'create_attribute' and this_at != 'display' and this_at != 'saveTo':
                baba = this_at.rjust(MAX_LENGTH) + " = " + str(getattr(self, this_at))
                print(baba)

    def saveTo(self, filename):
        """
        Function that saves the values of the attributes of the class to a file.
        @return: No returned value
        """

        fileToSaveTo = open(filename, 'w')
        fileToSaveTo.write("# parameters of the simulation\n")
        attributes = dir(self)
        for att in list(range(0,len(attributes))):
            this_at = attributes[att]
            if this_at[0] != '_' and this_at != 'create_attribute' and this_at != 'display' and this_at != 'saveTo':
                this_at_value = str(getattr(self, this_at))
                fileToSaveTo.write(this_at + " = " + this_at_value + '\n')
        fileToSaveTo.close()

