import argparse
parser = argparse.ArgumentParser(description="Searches alignment scoring matrix in file for highest match")
parser.add_argument( "-f", "--file", help="alignment score table file")
parser.add_argument( "-s", "--search", help="search for highest match of given individual")
parser.add_argument( "-a", "--all", help="print full content", action="store_true")
args = parser.parse_args()

print(f"searching {args.f}")
if args.s:
    print(f"printing best match for {args.s}")
if args.a:
    print("printing the entire table!!!")
