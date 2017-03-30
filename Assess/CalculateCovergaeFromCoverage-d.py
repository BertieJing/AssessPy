# import

import argparse
# code
parser = argparse.ArgumentParser(
    description='Script to change repeatmask coordinate based on dictionary of coordinates chromosome')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-i", "--input", help="read coverage using bedtools coverage -d",
                    type=str, required=True)
parser.add_argument("-o", "--output", help="output file",
                    type=str, required=True)

args = parser.parse_args()
input_file = args.input
output = args.output


class CalculateCoverage:
    def __init__(self, input_coverage, output_average_coverage):
        self.input_coverage = input_coverage
        self.output_average_coverage = output_average_coverage
        self.coverage_dictionary = {}

    def read_coverage(self):
        with open(self.input_coverage, 'r') as data_in:
            for line in data_in:
                line_split = line.strip().split("\t")
                begin = line_split[1]
                end = line_split[2]
                name = begin + ":" + end
                coverage = line_split[-1]
                if name in self.coverage_dictionary:
                    self.coverage_dictionary[name] = [self.coverage_dictionary[name][0] + int(coverage),
                                                      self.coverage_dictionary[name][1] + 1]
                else:
                    self.coverage_dictionary[name] = [int(coverage), int(1)]

    def write_average_coverage(self):
        with open(self.output_average_coverage, 'w') as data_out:
            for k, v in self.coverage_dictionary.iteritems():
                data_out.write("{}\t{}\n".format(k, v[0]/v[1]))
        print('Done')


if __name__ == "__main__":
    cal = CalculateCoverage(input_coverage=input_file, output_average_coverage=output)
    cal.read_coverage()
    cal.write_average_coverage()
