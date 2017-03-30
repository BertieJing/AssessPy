""" Script to extract uncorrected reads from proovread full read correction file
"""

from AssessUtil import AssessUtil
import argparse


parser = argparse.ArgumentParser(description='Script to extract uncorrected reads'
                                             ' from proovread full read correction file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')

parser.add_argument("-c", "--correctedReads", help="trimmed splitted proovread file",
                    type=str, required=True)
parser.add_argument("-f", "--fullLength", help="Proovread untrimmed unsplitted file",
                    type=str, required=True)
parser.add_argument("-o", "--outputFile", help="A file that contains uncorrected PacBio reads",
                    type=str, required=True)
parser.add_argument("-m", "--minimumLength", help="The minimum length of uncorrected PacBio reads default 50",
                    type=str, required=False)

args = parser.parse_args()

corrected = args.correctedReads
full_reads = args.fullLength
output = args.outputFile
minimum = args.minimumLength
if not minimum:
    minimum = 50


if __name__ == '__main__':
    AssessUtil.extract_uncorrected_reads_proovread(trimmed_corrected_reads=corrected,
                                                   full_corrected_reads=full_reads,
                                                   output_file=output, threshold=minimum)