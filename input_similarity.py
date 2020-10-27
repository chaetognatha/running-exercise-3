import argparse
parser = argparse.ArgumentParser(description="This is a description of what this program does")
parser.add_argument("input", help="name of input file")
parser.add_argument("name", help="the name to search")
parser.add_argument( "-ALL", help="prints all fields", action="store_true")
args = parser.parse_args()

input_file = args.input
if args.ALL:
    similarity_input = "ALL"
else:
    similarity_input = args.name

with open(input_file, 'r') as matrix:  # open the file

    if similarity_input == 'ALL':
        for element in matrix:
            print(element)
    else:
        head_list = matrix.readline().strip().split("\t")

        item_list = []
        for row in matrix:
            if similarity_input.capitalize() in row.capitalize():
                my_split = row.split("\t")
                my_max = max(my_split[1:])
                my_index = my_split.index(my_max)
                best_name = head_list[my_index]
                print(f"The best match is {best_name} with a score of {my_max}")
