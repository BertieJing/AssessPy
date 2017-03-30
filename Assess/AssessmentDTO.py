# import
from AssessUtil import AssessUtil
from ReadFasta import ReadFasta
# class definition


class AssessmentDTO:
    def __init__(self, corrected_reads, maf_file, corrected_file_type, align_type, log=None):
        self.corrected_reads = corrected_reads
        self.maf_file = maf_file
        if log is not None:
            self.log = log
        self.corrected_file_type = corrected_file_type
        self.align = align_type

    def get_maf_dictionary(self):
        maf_dict = AssessUtil()
        return maf_dict.read_maf(self.maf_file)

    def get_maf_sign_dictionary(self):
        maf_dict = AssessUtil()
        return maf_dict.read_maf_sign(self.maf_file)

    def get_lsc_dictionary(self):
        return ReadFasta.read_lsc_result(self.corrected_reads)

    def get_proovread_dictionary(self):
        return ReadFasta.read_proovread_result(self.corrected_reads)

    def get_lordec_dictionary(self):
        return ReadFasta.read_lordec_result(self.corrected_reads)

    def get_pbcr_dictionary(self):
        try:
            if self.log is not None:
                read_fasta = ReadFasta()
                return read_fasta.read_pbcr_result(self.corrected_reads, self.log)
        except AttributeError as pbcr_error:
            print(pbcr_error.message)

    def get_general_dictionary(self):
        return ReadFasta.read_general_corrected(self.corrected_reads)

    def get_corrected_file(self):
        if self.align.lower() == "local":
            if self.corrected_file_type.lower() == "lordec":
                return self.get_lordec_dictionary()
            elif self.corrected_file_type.lower() == "proovread":
                return self.get_proovread_dictionary()
            elif self.corrected_file_type.lower == "lsc":
                return self.get_lsc_dictionary()
            else:
                return self.get_pbcr_dictionary()
        else:
            if self.corrected_file_type.lower() == "lsc":
                return self.get_lsc_dictionary()
            else:
                return self.get_general_dictionary()
