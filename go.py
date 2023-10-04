import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Custom Constants Parser")

    parser.add_argument("--tab", default="    ", help="Custom indentation string")
    parser.add_argument("--type", default=": float", help="Default type annotation")
    parser.add_argument(
        "--std", choices=["WRITE", "PRINT"], default="WRITE", help="Output mode"
    )
    parser.add_argument("--outfile", default="vakyume.py", help="Output file name")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    print(f"TAB: {args.tab}")
    print(f"TYPE: {args.type}")
    print(f"STD: {args.std}")
    print(f"OUTFILE: {args.outfile}")
