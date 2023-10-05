# These were used to convert notes to code and check formatting
# This class can be used on a .pyeqn file to confirm which equations are valid
import re, os


class PyEqnValidator:
    def __init__(self, infile: str):
        self.infile = infile
        self.eqn_regex = "\d{1,2}-\d{1,2}(?:\w)?"

    def reveal_blank_eqn_names(self):
        ix = 1
        with open(os.getcwd() +'/' + self.infile) as file:
            for l in file.readlines():
                # Use regular expression to find patterns like "1-2a"
                # print(l)
                if eqn_number := re.compile(self.eqn_regex).findall(l):
                    print(eqn_number)
                # Use regular expression to find patterns like "1-2"
                if x := re.compile(self.eqn_regex).findall(l):
                    eqn_number = x[0]
                # Check if the line contains "=" but not followed by ":"
                if "=" in l and ":" not in l.split("=")[0]:
                    self.parse_eqn(l.strip().split("#")[0])

    def parse_eqn(self, l: str):
        try:
            # Try to evaluate the expression in the line
            eval(l)
        except SyntaxError as se:
            # If there's a syntax error, extract the tokens and print the line
            to_tokens = se.text.split("=")[0]
            y = re.compile("\w[A-Za-z0-9\-\_]+").findall(to_tokens)
            print(l)
        except NameError as ne:
            # If there's a name error, ignore it
            pass

if __name__ == '__main__':
    infile = 'pyeqn_examples/simple.pyeqn'
    pev = PyEqnValidator(infile)
    pev.reveal_blank_eqn_names()