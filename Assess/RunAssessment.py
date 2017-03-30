# import
import os, sys
from multiprocessing import Process
from AssessmentDTO import AssessmentDTO
from AssessUtil import AssessUtil
from AssessAlign import AssessAlign
import time
# class


class RunAssessment:
    def __init__(self, corrected_fasta, maf_file, corrected_file_type, assess_result_file, log_file=None, work_dir="./"
                 , number_of_threads=1, align="local"):
        print("\n######\nRunning Info:\n")
        print("Corrected File: {}\nMaf File: {}\nSoftware used: {}\nWorking directory: {}\nLog file: {}"
              "\nNumber of threads running: {}\nAlign Method: {}\n#####\n".format(corrected_fasta, maf_file,
                                                                                  corrected_file_type, work_dir,
                                                                                  log_file, number_of_threads, align))
        self.number_of_threads = number_of_threads
        self.assess_result_file = assess_result_file
        self.align = align
        self.corrected_file_type = corrected_file_type
        if corrected_file_type.lower() == "pbcr":
            if log_file is not None:
                self.log_file = log_file
            else:
                sys.exit('to assess corrected result form PBcR you need log file')
        assess_dto = AssessmentDTO(corrected_fasta, maf_file, corrected_file_type, align, log_file)
        print("#############\n Reading Maf file \n#############")
        self.maf_dictionary = assess_dto.get_maf_sign_dictionary()

        # create working directories
        root_directory = work_dir
        try:
            os.makedirs(root_directory+"/result/assess")
            os.makedirs(root_directory+"/result/run")
            os.makedirs(root_directory+"/result/miss")

            # declare directories
            self.assess_dir = root_directory+"/result/assess/"
            self.run_dir = root_directory+"/result/run/"
            self.miss_dir = root_directory+"/result/miss/"
        except OSError:
            if not os.path.isdir(root_directory):
                raise
        print("#############\n Reading corrected file \n#############")
        if corrected_file_type.lower() == "pbcr":
            self.corrected_dictionary, self.lost_data = assess_dto.get_corrected_file()
            with open(self.miss_dir + self.assess_result_file + "_missed.txt", 'w') as missed:
                missed.write("amount of missed data in case of pbcr = {}".format(self.lost_data))
        else:
            self.corrected_dictionary = assess_dto.get_corrected_file()

    def correct_assess(self):
        start = time.clock()
        start_wall_time = time.time()
        assess_subdirectory = AssessUtil.chunks(self.corrected_dictionary,
                                                len(self.corrected_dictionary)/int(self.number_of_threads))
        # loop list and run threading
        my_process = []
        assess_files = []
        miss_files = []
        i = 1
        print("#############\n Initialize assessment  \n#############")
        for dic in assess_subdirectory:
            assess_files.append(str(i)+"_"+self.assess_result_file)
            miss_files.append(str(i)+self.assess_result_file + "_missed.txt")
            my_process.append(Process(target=self.run_assess, args=(dic, str(i)+"_"+self.assess_result_file,)))
            i += 1
        for pr in my_process:
            pr.start()
        for pr in my_process:
            pr.join()
        print("#############\n Begin assessment  \n#############")
        with open(self.assess_dir + self.assess_result_file, 'a') as outfile:
            outfile.write("sequence\tsubsequence\tseqlength\tsubLength\tsubsimilarity\tsizeOfAllSub\n")
            for file_name in assess_files:
                with open(self.assess_dir+file_name) as infile:
                    for line in infile:
                        outfile.write(line)
        end = time.clock()
        end_wall_time = time.time()
        print("#################\n process finished taking {} processor seconds\n and {} wall time\n#################"
              .format(end-start, end_wall_time - start_wall_time))

    def run_assess(self, correct_dictionary, assess_file):
        if self.corrected_file_type.lower() == "lsc" and self.align.lower() == "global":
            self.assess_lsc(correct_dictionary, assess_file)
        else:
            if self.align.lower() == "local":
                self.assess_other(correct_dictionary, assess_file)
            else:
                self.assess_other_global(correct_dictionary, assess_file)

    def assess_other(self, correct_dictionary, assess_file):
        fo = open(self.assess_dir + assess_file, 'w')
        print("#############\n writing to File {} \n#############".format(assess_file))
        for k in correct_dictionary.keys():
            value = correct_dictionary[k]
            size_of_all_subreads = int(value[0][0])
            original_seq = self.maf_dictionary[k][1]
            sequence_strand = self.maf_dictionary[k][0]
            original_sequence_length = len(original_seq)
            for element in value[1:]:
                seq_length = element[0]
                seq_name = element[1]
                if self.corrected_file_type.lower() == "pbcr":
                    seq_value = str(AssessUtil.revers_complement(element[3])) if int(sequence_strand) < 0 else element[3]
                else:
                    seq_value = str(AssessUtil.revers_complement(element[2])) if int(sequence_strand) < 0 else element[2]
                # run local alignment
                corrected_local_similarity = AssessAlign.water_align_code(seq_value, original_seq)
                fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(k, seq_name, original_sequence_length, seq_length,
                                                         corrected_local_similarity, size_of_all_subreads))
        fo.close()
        print("#############\n closing to File {} \n#############".format(assess_file))

    def assess_other_global(self, correct_dictionary, assess_file):
        fo = open(self.assess_dir + assess_file, 'w')
        print("#############\n writing to File {} \n#############".format(assess_file))
        for k in correct_dictionary.keys():
            original_seq = self.maf_dictionary[k][1]
            sequence_strand = self.maf_dictionary[k][0]
            seq_value = str(AssessUtil.revers_complement(correct_dictionary[k])) if int(sequence_strand) < 0 else correct_dictionary[k]
            size_of_all_subreads = len(seq_value)
            original_sequence_length = len(original_seq)
            corrected_global_similarity = AssessAlign.needle_align_code(seq_value, original_seq)
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(k, k, original_sequence_length, size_of_all_subreads,
                                                         corrected_global_similarity, size_of_all_subreads))
        fo.close()
        print("#############\n closing to File {} \n#############".format(assess_file))

    def assess_lsc(self, correct_dictionary, assess_file):
        fo = open(self.assess_dir + assess_file, 'w')
        print("#############\n writing to File {} \n#############".format(assess_file))
        for k in correct_dictionary.keys():
            value = correct_dictionary[k]
            original_seq = self.maf_dictionary[k][1]
            sequence_strand = self.maf_dictionary[k][0]
            original_sequence_length = len(original_seq)
            size_of_all_subreads = original_sequence_length
            seq_value = AssessUtil.revers_complement(value) if int(sequence_strand) < 0 else value
            # run local alignment
            corrected_global_similarity = AssessAlign.needle_align_code(seq_value, original_seq)
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(k, k, original_sequence_length, len(value),
                                                       corrected_global_similarity, size_of_all_subreads))
        fo.close()
        print("#############\n closing to File {} \n#############".format(assess_file))

    @staticmethod
    def assess_maf(maf, threads, directory, output, align="global"):
        print("\n######\nRunning Info:\n")
        print("Maf: {}\nThreads: {}\nDirectory: {}\nOutput: {}\nAlign: {}\n#########\n\n".format(maf, threads,
                                                                                                 directory, output,
                                                                                                 align))
        start = time.clock()
        start_wall_time = time.time()
        my_process = []
        assess_files = []
        i = 1
        maf_dictionary = AssessUtil.read_maf(maf)
        os.makedirs(directory+"/result/assess")

        # declare directories
        assess_dir = directory+"/result/assess/"
        assess_subdirectory = AssessUtil.chunks(maf_dictionary,
                                                len(maf_dictionary)/int(threads))
        print("#############\n Initialize assessment  \n#############")
        for dic in assess_subdirectory:
            assess_files.append(str(i) + "_" + output)
            print("#############\n Begin assessment batch {} \n#############".format(i))
            my_process.append(Process(target=RunAssessment.run_maf_assess, args=(dic, str(i)+"_"+output, assess_dir,)))
            i += 1
        for pr in my_process:
            pr.start()
        for pr in my_process:
            pr.join()

        with open(assess_dir + output, 'a') as outfile:
            outfile.write("seq_name\toriginal_seq_len\tartificial_seq_length\tglobal_align_value\tgc\tcomplexity\n")
            for file_name in assess_files:
                with open(assess_dir+file_name) as infile:
                    for line in infile:
                        outfile.write(line)
        end = time.clock()
        end_wall_time = time.time()
        print("#################\n process finished taking {} processor seconds\n and {} wall time\n#################"
              .format(end-start, end_wall_time - start_wall_time))

    @classmethod
    def run_maf_assess(cls, dictionary, output, assess_dir, align="global"):
        print("#############\n writing to File {} \n#############".format(output))
        fo = open(assess_dir+output, 'w')
        assess_util = AssessUtil()
        for k, v in dictionary.iteritems():
            seq_name = k
            original_seq = v[0]
            artificial_seq = v[1]
            gc_content = assess_util.gc(original_seq)
            complexity = assess_util.complexity(original_seq)
            corrected_global_similarity = AssessAlign.needle_align_code(artificial_seq, original_seq)
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(seq_name, len(original_seq), len(artificial_seq),
                                                       corrected_global_similarity, gc_content,
                                                       complexity))

        print("#############\n closing to File {} \n#############".format(output))
        fo.close()

if __name__ == "__main__":
    print("running main")