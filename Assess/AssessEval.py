# import

from AssessmentDTO import AssessmentDTO
from AssessAlign import AssessAlign
import time

# class definition


class AssessEval:
    def __init__(self):
        pass

    @classmethod
    def eval_mafft_code(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001_test.maf")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.mafft_align(v[1], v[0], k, k+"_ref")
        end = time.clock()
        print("overall time for code mafft: {} ".format(end - start))

    @classmethod
    def eval_mafft_command(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001_test.maf", corrected_file_type="")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.mafft_align_command_line(v[1], v[0], k, k+"_ref",
                                                 directory="/home/medhat/Desktop/Desktop_files/asses_test/exampl/")
        end = time.clock()

        print("overall time for command local mafft: {} ".format(end - start))

    @classmethod
    def eval_mafft_code_global(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001.maf")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.mafft_align(v[1], v[0], k, k+"_ref", "global")
        end = time.clock()
        print("overall time for code mafft global alignment: {} ".format(end - start))

    @classmethod
    def eval_mafft_command_global(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001.maf", corrected_file_type="")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.mafft_align_command_line(v[1], v[0], k, k+"_ref",  "global")
        end = time.clock()
        print("overall time for command mafft global alignment: {} ".format(end - start))

    @classmethod
    def eval_fogsaa_command_global(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001_test.maf")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.fogsaa_align_command_line(str(v[1]).strip(), str(v[0]).strip(), k, k+"_ref",
                                                  directory="/home/medhat/mytest/")
        end = time.clock()
        print("overall time for command fogsaa global alignment: {} ".format(end - start))

    @classmethod
    def eval_needle_global(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001_test.maf", corrected_file_type="")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.needle_align_code(str(v[1]).strip(), str(v[0]).strip())
        end = time.clock()
        print("overall time for command needle global alignment: {} ".format(end - start))

    @classmethod
    def eval_water_local(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001_test.maf", corrected_file_type="")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.water_align_code(str(v[1]).strip(), str(v[0]).strip())
        end = time.clock()
        print("overall time for command water local alignment: {} ".format(end - start))

    @classmethod
    def eval_water_local_command(cls):
        assessDto = AssessmentDTO("/home/medhat/Desktop/Desktop_files/asses_test/exampl/lol.txt",
                                  "/home/medhat/Downloads/ecoli_pacbio.fa_0001_test.maf", corrected_file_type="")
        maf_dictionary = assessDto.get_maf_dictionary()
        start = time.clock()
        for k, v in maf_dictionary.iteritems():
            AssessAlign.water_align_command(str(v[1]).strip(), str(v[0]).strip())
        end = time.clock()
        print("overall time for command water local alignment: {} ".format(end - start))

if __name__ == "__main__":
    # AssessEval.eval_mafft_code()

    # AssessEval.eval_mafft_code_global()
    # AssessEval.eval_mafft_command_global()
    # AssessEval.eval_fogsaa_command_global()
    # AssessEval.eval_needle_global()
    # AssessEval.eval_water_local_command()
    # AssessEval.eval_mafft_command()
    AssessEval.eval_water_local()
