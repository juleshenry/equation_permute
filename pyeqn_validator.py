# These were used to convert notes to code and check formatting
# This class can be used on a .pyeqn file to confirm which equations are valid
import re, os


class PyEqnValidator:
    def __init__(self, infile: str):
        self.infile = infile

    def reveal_blank_eqn_names(self):
        ix = 1
        for o in os.listdir(os.getcwd() + self.infile):
            if o[0] == "_":
                continue
            with open(os.getcwd() + self.infile) as file:
                for l in file.readlines():
                    if x := re.compile("\d{1,2}-\d{1,2}\w").findall(l):
                        eqn_number = x[0]
                        if len(l) < 10:
                            ix += 1
                            # print(ix, l.strip(), "needs name!")

    def see_which_notes_are_valid_Python(self):
        for o in os.listdir(os.getcwd() + "/chapters"):
            if o[0] == "_":
                continue
            with open(os.getcwd() + "/chapters/" + o) as file:
                eqn_number = ""
                for l in file.readlines():
                    if x := re.compile("\d{1,2}-\d{1,2}").findall(l):
                        eqn_number = x[0]
                    if "=" in l and ":" not in l.split("=")[0]:
                        self.parse_eqn(l.strip().split("#")[0])

    def parse_eqn(self, l: str):
        try:
            eval(l)
        except SyntaxError as se:
            to_tokens = se.text.split("=")[0]
            y = re.compile("\w[A-Za-z0-9\-\_]+").findall(to_tokens)
            for o in y:
                print(o, "=", 1.1)
            print(l)
        except NameError as ne:
            pass
