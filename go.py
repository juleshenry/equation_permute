import argparse
from equation_permuter import main

def parse_arguments():
    parser = argparse.ArgumentParser(description="Custom Constants Parser")
    parser.add_argument("--tab", default="    ", help="Custom indentation string")
    parser.add_argument("--typehint", default=": float", help="Default type annotation")
    parser.add_argument("--infile", default="in_file.py", help="Input file name")
    parser.add_argument("--outfile", default="out_file.py", help="Output file name")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main(args.tab, args.typehint, args.infile, args.outfile)