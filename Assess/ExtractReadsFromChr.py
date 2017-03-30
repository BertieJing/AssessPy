# import
import os
import uuid
import pickle
import ntpath
from pyfaidx import Fasta
from AssessUtil import AssessUtil
from Blat import Blat

# code


class ExtractReadsFromChr:
    def __init__(self, maf_file, working_dir, corrected_file=None):
        self.maf = maf_file
        if corrected_file:
            self.corrected_file = corrected_file
        else:
            self.corrected_file = ""
        self.work_dir = working_dir
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)
        self.original_seq_file = self.work_dir + "/original_seq.fa"
        # this returns file name for a dictionary contains read and coordinate
        self.original_coordinate = AssessUtil.extract_original_read_maf(maf_file=self.maf, maf_fasta=self.original_seq_file)
        self.original_coordinate = pickle.load(open(self.original_coordinate, "rb"))
        self.original_sequence = Fasta(self.original_seq_file)

    def get_corrected_coordinate(self):
        corrected_reads_coordinate = {}
        with Fasta(self.corrected_file) as data_in:
            for line in data_in:
                corrected_seq_name = line.name
                original_seq_name = line.name.split('.', 1)[0]
                original_sequence_bases = self.original_sequence[original_seq_name]

                # create two files to use in blat align
                unique_filename = str(uuid.uuid4())
                db_file = self.work_dir + "/" + unique_filename + "_ref.fa"
                query_file = self.work_dir + "/" + unique_filename + "_query.fa"

                with open(db_file, 'w') as ref_data, open(query_file, 'w') as query_data:
                    ref_data.write(">{}\n{}".format(original_seq_name, original_sequence_bases))
                    query_data.write(">{}\n{}".format(corrected_seq_name, str(line)))

                blat = Blat(db=db_file, query=query_file, align_type="psl", header=False)
                blat.run_blat()
                align_result = blat.parse_output()

                offset = int(self.original_coordinate[original_seq_name][0])
                begin = offset + int(align_result[15])
                end = offset + int(align_result[16])
                if original_seq_name in corrected_reads_coordinate:
                    corrected_reads_coordinate[original_seq_name].append([begin, end])
                else:
                    corrected_reads_coordinate[original_seq_name] = []
                    corrected_reads_coordinate[original_seq_name].append([begin, end])
                os.remove(db_file)
                os.remove(query_file)
        return corrected_reads_coordinate

    def get_corrected_uncorrected_coordinate(self):
        # run blat for the maf and corrected db
        corrected_basename = ntpath.basename(self.corrected_file).split('.', 1)[0]
        my_psl_output = self.work_dir + "/" + corrected_basename + ".psl"
        blat = Blat(db=self.original_sequence, query=self.corrected_file, original_coordinate=self.original_coordinate,
                    align_type="psl", header=False, output=my_psl_output)
        blat.run_blat()
        corrected_read_coordinate = blat.get_corrected_read_coordinate()
        uncorrected_read_coordinate = self.get_uncorrected_coordinate(corrected_read_coordinate)
        return corrected_read_coordinate, uncorrected_read_coordinate

    def get_uncorrected_coordinate(self, dictionary_of_corrected_reads):
        intervals_of_uncorrected_reads = AssessUtil.uncorrected_interval(
            corrected_coordinate=dictionary_of_corrected_reads, original_seq_length=self.original_coordinate)
        return intervals_of_uncorrected_reads

if __name__ == "__main__":
    a = ExtractReadsFromChr(maf_file='/home/medhat/Downloads/yeast4_pacbio.fa_0001.maf',
                            corrected_file='/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/yeast.trimmed.fa',
                            working_dir='/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/test_coordent'
                            )
    corrected, uncorrected = a.get_corrected_uncorrected_coordinate()



