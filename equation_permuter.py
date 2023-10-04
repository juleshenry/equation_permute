import os
import re
import time
from sympy import Symbol, solve, log
import subprocess

TAB = "    "
TYPE = ": float"
STD = "WRITE"
OUTFILE = "vakyume2.py"
_DEBUG = True
ROOT = 'tests'

def debug(s:str):
    if _DEBUG:
        print("[DEBUG]",s)

def stdout(s):
    if STD == "WRITE":
        with open(OUTFILE, "a+") as o:
            o.write(s + "\n")
    else:
        print(s)


class Solver:

    def valid_toke(s, t):
        # print(t)
        valid = t.isidentifier() | t.split('**')[0].strip().isidentifier()
        # print(t.split('**')[0].strip().isidentifier(),'~',valid)
        return valid

    def clean_t(s,t):
        d = t.strip().replace("(", "").replace(")", "")
        o = d.split('**')
        if len(o)>1 and o[1]:
            d = f"{t.split('**')[0]} ** {''.join(o[1:])}" 
            # print(d)
        return d

    def get_tokes(s, eqn):
        """
        Assumes equations are symbol-separated by spaces
        Anything like math.sin(a**2) must be treated via conversion to Sympy syntax
        B**2 is treated a token by space separation.
        `clean toke` cleans v**2 to v ** 2 for processing by Symbol
        """
        tokes = set()
        for t in eqn.split(" "):
            clean = s.clean_t(t)
            # print(clean,clean.isidentifier())
            if s.valid_toke(clean) and t not in {"ln", "log"}:
                tokes.add(clean)
            # else:
            #     print("invalid toke",clean)
        tokes = list(tokes)
        # print('tokes',tokes)
        return tokes

    def permute(self, eqn: str, eqn_n: str):
        tokens = self.get_tokes(eqn)
        normal_form = eqn.split("=")[1].strip() + " - " + eqn.split("=")[0].strip()
        for token in tokens:
            # print("investigating",t)
            args = sorted(filter(lambda x: x != token, tokens))
            args = list(map(lambda x:x.split('**')[0].strip(), args))
            typed_args = str(f"{TYPE}, ").join(args)
            if typed_args:
                typed_args += TYPE
            print(typed_args)
            token = token.split('**')[0]
            stdout(f"{TAB}@staticmethod")
            stdout(f'{TAB}def eqn_{eqn_n.replace("-","_")}__{token}({typed_args}):')
            stdout(f"{TAB}# {eqn.strip().replace('#','')}")
            # try:
            print('*'*88)
            print(normal_form)
            print(Symbol(token))
            solns = solve(normal_form, (token))
            print(solns)
            print('*'*88)
            # print('NORM',normal_form)
            if not len(solns):
                print('FAILED on:',fr'{normal_form}',rf'{token}',sep='\n')
                stdout(f"{TAB*2}pass # unable  to solve")
                continue
            stdout(TAB * 2 + "result = []")
            for soln in solns:
                stdout(f"{TAB*2}{token} = {soln}")
                stdout(f"{TAB*2}result.append({token})")
            stdout(TAB * 2 + f"return {token}")
            # except:
            # stdout(f"{TAB*2}pass #NotImplementedError")

    def analyze(s, i):
        root_dir = os.getcwd() + f'//{ROOT}//'
        get = list(filter(lambda x: i in x, os.listdir(root_dir)))[0]
        with open(root_dir + get) as file:
            eqn_number = ""
            for l in file.readlines():
                if x := re.compile("\d{1,2}-\d{1,2}\w{,2}").findall(l):
                    eqn_number = x[0]
                if " = " in l:
                    # print("[DEBUG]",eqn_number)
                    s.permute(l, eqn_number)


class SetupMethods:
    # These were used to convert my notes to code and check formatting
    # Universally, they are the 1st step in the process
    @staticmethod
    def reveal_blank_eqn_names():
        ix = 1
        for o in os.listdir(os.getcwd() + "/chapters"):
            if o[0] == "_":
                continue
            with open(os.getcwd() + "/chapters/" + o) as s:
                for l in s.readlines():
                    if x := re.compile("\d{1,2}-\d{1,2}\w").findall(l):
                        eqn_number = x[0]
                        if len(l) < 10:
                            ix += 1
                            # print(ix, l.strip(), "needs name!")

    def see_which_notes_are_valid_Python(s):
        for o in os.listdir(os.getcwd() + "/chapters"):
            if o[0] == "_":
                continue
            with open(os.getcwd() + "/chapters/" + o) as s:
                eqn_number = ""
                for l in s.readlines():
                    if x := re.compile("\d{1,2}-\d{1,2}").findall(l):
                        eqn_number = x[0]
                    if "=" in l and ":" not in l.split("=")[0]:
                        s.parse_eqn(l.strip().split("#")[0])

    def parse_eqn(s, l):
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


def get_class_name(mods: str):
    # ends in .pyeqn, length 6
    return "".join(x[0].upper() + x[1:] for x in mods)[:-6]

if __name__ == "__main__":
    X = Solver()
    for module in sorted(os.listdir(os.getcwd() + f"/{ROOT}")):
        if not module.startswith('quad'):
            continue
        chap, mods = module.split("_")[0], module.split("_")[1:]
        cls_name = get_class_name(mods)
        stdout(f"class {cls_name}:")
        print(cls_name)
        X.analyze(chap)
    # Run black on the file
    subprocess.run(['black', OUTFILE])