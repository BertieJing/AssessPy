""" Script to calculate complexity and GC content of read
"""

from AssessUtil import AssessUtil
import argparse
import time
from pyfaidx import Fasta


parser = argparse.ArgumentParser(description='Script to calculate complexity and GC content of read')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')

parser.add_argument("-i", "--input", help="Fatsa reads file",
                    type=str, required=True)
parser.add_argument("-o", "--outputFile", help="output reads file",
                    type=str, required=True)


args = parser.parse_args()

input_file = args.input
output_file = args.outputFile

if __name__ == '__main__':
    with Fasta(input_file) as data_in, open(output_file, 'w') as data_out:
        data_out.write("read_name\tGC\tcomplexity\n")
        for line in data_in:
            seq_name = line.name
            seq = str(line)
            gc = AssessUtil.gc_content(seq)
            complexity = AssessUtil.complexity(seq)
            data_out.write("{}\t{}\t{}\n".format(seq_name, gc, complexity))
    print('Done {} {}'.format(time.strftime("%d/%m/%Y"), time.strftime("%I:%M:%S")))