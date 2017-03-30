# import
import argparse
import pickle
import os
from AssessUtil import AssessUtil
# code
parser = argparse.ArgumentParser(
    description='Script to extract coordinate from chromosome')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-c", "--coordinate", help="original coordinate",
                    type=str, required=True)
parser.add_argument("-b", "--blatResult", help="blat file",
                    type=str, required=True)

args = parser.parse_args()
original_coordinate = args.coordinate
blat_align_result = args.blatResult


class GetCoordinates:
    def __init__(self, original_read_coordinate, blat_result):
        self.original_read_coordinate = pickle.load(open(original_read_coordinate, 'rb'))
        self.blat_result = blat_result

    def get_corrected_read_coordinate(self):
        corrected_read_dictionary = {}
        with open(self.blat_result, 'r') as input_data:
            for line in input_data:
                line_split = line.split()
                cor_seq_name = line_split[9]
                original_seq_name = line_split[13]
                begin = int(line_split[15])
                end = int(line_split[16])
                if cor_seq_name.split('.', 1)[0] == original_seq_name:
                    # convert to original chromosome coordinate
                    begin = begin + self.original_read_coordinate[original_seq_name][0]
                    end = end + self.original_read_coordinate[original_seq_name][0]
                    if original_seq_name in corrected_read_dictionary:
                        corrected_read_dictionary[original_seq_name].append([begin, end])
                    else:
                        corrected_read_dictionary[original_seq_name] = []
                        corrected_read_dictionary[original_seq_name].append([begin, end])

        return AssessUtil.sort_dictionary_values(unsorted_dictionary=corrected_read_dictionary)

    def get_uncorrected_coordinate(self):
        intervals_of_uncorrected_reads = AssessUtil.uncorrected_interval(
            corrected_coordinate=self.get_corrected_read_coordinate(), original_seq_length=self.original_read_coordinate)
        return intervals_of_uncorrected_reads

if __name__ == "__main__":
    get_coordinate = GetCoordinates(original_read_coordinate=original_coordinate, blat_result=blat_align_result)
    corrected_coordinate_dictionary = get_coordinate.get_corrected_read_coordinate()
    uncorrected_coordinate_dictionary = get_coordinate.get_uncorrected_coordinate()
    directory = os.path.dirname(blat_align_result)
    pickle.dump(corrected_coordinate_dictionary, open(directory + "/corrected_read_coordinate.p", 'wb'))
    pickle.dump(uncorrected_coordinate_dictionary, open(directory + "/uncorrected_read_coordinate.p", 'wb'))