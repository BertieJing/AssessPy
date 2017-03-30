# import
import pickle
import argparse
# code
parser = argparse.ArgumentParser(
    description='Script to change repeatmask coordinate based on dictionary of coordinates chromosome')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-c", "--coordinate", help="original coordinate",
                    type=str, required=True)
parser.add_argument("-b", "--bed", help="bed file of the repeatmask result",
                    type=str, required=True)
parser.add_argument("-o", "--output", help="output file",
                    type=str, required=True)

args = parser.parse_args()
original_coordinate = args.coordinate
repeat_mask_bed = args.bed
output = args.output


class ChangeBedCoordinate:
    def __init__(self, original_coordinate_file, bed_file, output):
        self.original_coordinate_dictionary = pickle.load(open(original_coordinate_file, 'rb'))
        self.bed_file = bed_file
        self.output = output

    def change_coordinate(self):
        with open(self.bed_file, 'r') as data_in, open(self.output, 'w') as data_out:
            for line in data_in:
                line_split = line.split()
                offset = self.original_coordinate_dictionary[line_split[0]][0]
                line_split[1] = offset + int(line_split[1])
                line_split[2] = offset + int(line_split[2])
                data_out.write("\t".join(str(v) for v in line_split)+"\n")
        print('Done')

if __name__ == "__main__":
    a = ChangeBedCoordinate(original_coordinate_file=original_coordinate,
                            bed_file=repeat_mask_bed,
                            output=output)
    a.change_coordinate()