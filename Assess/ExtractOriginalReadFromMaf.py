""" Script to extract original read  from maf file
"""
import argparse
from AssessUtil import AssessUtil

parser = argparse.ArgumentParser(
    description='Script to extract original read  from maf file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-m", "--mafFile", help="A maf file that contains the alignment of original and artificial data",
                    type=str, required=True)
parser.add_argument("-o", "--outputFile", help="A file that contains assessment of the corrected PacBio reads",
                    type=str, required=True)

args = parser.parse_args()
maf = args.mafFile
result = args.outputFile

if __name__ == "__main__":
    original_coordinate = AssessUtil.extract_original_read_maf(maf_file=maf, maf_fasta=result)
    print("Dictionary with the coordinate in: {}".format(original_coordinate))
