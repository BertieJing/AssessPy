# import
import os
import pickle
import argparse
# code
parser = argparse.ArgumentParser(
    description='Script to change dictionary format to bed file ')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-c", "--coordinate", help="dictionary of  coordinate",
                    type=str, required=True)
parser.add_argument("-o", "--output", help="output file",
                    type=str, required=True)
parser.add_argument("-ch", "--chromosome", help="chromosome to be written in bed",
                    type=str, required=True)

args = parser.parse_args()
original_coordinate = args.coordinate
output_result = args.output
chromosome = args.chromosome


class DictionaryToBed:
    def __init__(self, coordinate_dictionary, chromosome, output=None):
        self.coordinate_dictionary = coordinate_dictionary
        self.chromosome = chromosome
        if output:
            self.output = output
        else:
            directory = os.path.dirname(coordinate_dictionary)
            basename = os.path.basename(coordinate_dictionary).split('.', 1)[0]
            self.output = os.path.join(directory, basename+'.bed')

    def from_dictionary_to_bed(self):
        coordinate_dictionary = pickle.load(open(self.coordinate_dictionary, 'rb'))
        with open(self.output, 'w') as data_out:
            for k, v in coordinate_dictionary.iteritems():
                i = 1
                for element in v:
                    data_out.write("{}\t{}\t{}\t{}\n".format(self.chromosome, element[0], element[1], k + '.' + str(i)))
                    i += 1

if __name__ == "__main__":
    dictionary_to_bed = DictionaryToBed(coordinate_dictionary=original_coordinate, chromosome=chromosome, output=output_result)
    dictionary_to_bed.from_dictionary_to_bed()
    print('Done')