# import
import os
import pickle
import ntpath
from collections import defaultdict
from pyfaidx import Fasta
from itertools import islice
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.lcc import lcc_simp
# class definition


class AssessUtil:
    def __init__(self):
        pass

    @classmethod
    def read_maf(cls, maf_input_file):
        """
        read maf file
        maf_input_file is file
        """
        maf_dictionary_reads = defaultdict(list)
        with open(maf_input_file, 'r') as input_maf_file:
            for line in input_maf_file:
                if line.strip().startswith("a"):
                    ref_line = next(input_maf_file).strip().split()
                    art_line = next(input_maf_file).strip().split()
                    ref_seq = ref_line[6].replace("-", "")
                    art_seq = art_line[6].replace("-", "")

                    art_strand = "+1" if art_line[4] == '+' else "-1"
                    art_seq_name = art_line[1]
                    maf_dictionary_reads[art_seq_name] = [ref_seq, art_seq, art_strand]
        return maf_dictionary_reads

    @classmethod
    def read_maf_align(cls, maf_input_file):
        """
        read maf file
        maf_input_file is file
        """
        maf_dictionary_reads = defaultdict(list)
        with open(maf_input_file, 'r') as input_maf_file:
            for line in input_maf_file:
                if line.strip().startswith("a"):
                    ref_line = next(input_maf_file).strip().split()
                    art_line = next(input_maf_file).strip().split()
                    ref_seq = ref_line[6]
                    art_seq = art_line[6]
                    art_strand = "+1" if art_line[4] == '+' else "-1"
                    art_seq_name = art_line[1]
                    maf_dictionary_reads[art_seq_name] = [ref_seq, art_seq, art_strand]
        return maf_dictionary_reads

    @staticmethod
    def extract_original_read_maf(maf_file, maf_fasta):
        coordinate = {}
        with open(maf_file, 'r') as input_data, open(maf_fasta, 'w') as data_out:
            for line in input_data:
                # if new alignment
                if line.startswith('a'):
                    ref_line = next(input_data).strip().split()
                    begin = int(ref_line[2])
                    end = begin + int(ref_line[3])
                    ref_seq = ref_line[6].replace('-', '')
                    seq_name = next(input_data).strip().split()[1]
                    seq_length = len(ref_seq)
                    coordinate[seq_name] = [begin, seq_length]  # also it is better to use sequence end
                    data_out.write(">{} {}:{}\n{}\n".format(seq_name, begin, end, ref_seq))
        out_dir = os.path.dirname(maf_fasta)
        base_name = ntpath.basename(maf_fasta).split('.', 1)[0]
        file_name = out_dir+"/"+base_name+".dictionary"
        pickle.dump(coordinate, open(file_name, "wb"))
        return file_name

    @classmethod
    def read_maf_sign(cls, maf_input_file):
        """
        read maf file maf_input_file
        """
        maf_dictionary_reads = defaultdict(list)
        with open(maf_input_file, 'r') as input_maf_file:
            for line in input_maf_file:
                if line.strip().startswith("a"):
                    ref_line = next(input_maf_file).strip().split()
                    art_line = next(input_maf_file).strip().split()
                    ref_seq = ref_line[6].replace("-", "")
                    art_strand = "+1" if art_line[4] == '+' else "-1"
                    art_seq_name = art_line[1]
                    maf_dictionary_reads[art_seq_name] = [art_strand, ref_seq]
        return maf_dictionary_reads

    @staticmethod
    def read_repeat_mask_file(self, repeat_file):
        """
        read masked reads
        """
        if repeat_file:
            try:
                repeat_data = {}

                with open(repeat_file, 'r') as data_in:
                    for line in data_in:
                        seq_name, start, end, repeat_type = line.split()
                        start = int(start)
                        end = int(end)
                        if seq_name in repeat_data:
                            if len(repeat_data[seq_name]) < start:
                                repeat_data[seq_name][len(repeat_data[seq_name]):start] = \
                                    'N' * (start - len(repeat_data[seq_name])-1)
                            if repeat_type == 'Simple_repeat':
                                repeat_data[seq_name][start:end] = 'S' * (end - start+1)
                            else:
                                repeat_data[seq_name][start:end] = 'H' * (end - start+1)
                        else:
                            repeat_data[seq_name] = []
                            if start > 0:
                                    # -1 cause the alignment is 1 based
                                    repeat_data[seq_name][0:(start-1)] = 'N' * (start-1)
                            if repeat_type == 'Simple_repeat':
                                repeat_data[seq_name][start:end] = 'S' * (end - start+1)
                            else:
                                repeat_data[seq_name][start:end] = 'H' * (end - start+1)

                    return repeat_data
            except IOError as e:
                print(e.message)

    def gc(self, seq):
        """ Return the GC content of seq as a float """
        g = seq.count('G')
        g += seq.count('g')
        c = seq.count('C')
        c += seq.count('c')
        return (g + c) / float(len(seq))

    @staticmethod
    def gc_content(seq):
        gc = sum(seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        try:
            return gc * 100.0 / len(seq)
        except ZeroDivisionError:
            return 0.0

    @staticmethod
    def chunks(data, SIZE=10000):
        it = iter(data)
        dictionary_list=[]
        for i in xrange(0, len(data), int(SIZE)):
            dictionary_list.append({k:data[k] for k in islice(it, SIZE)})
        return dictionary_list

    @staticmethod
    def revers_complement(seq):
        my_dna = Seq(seq, generic_dna)
        return my_dna.reverse_complement()

    @staticmethod
    def complexity(seq):
        return lcc_simp(seq)

    @staticmethod
    def get_number_of_unique_reads_in_lordec(lordec_corrected_file):
        unique_reads = set()
        with Fasta(lordec_corrected_file) as data_in:
            for line in data_in:
                if line.name.count('_') > 1:
                    unique_reads.add("_".join(line.name.split("_", 2)[:2]))
                else:
                    unique_reads.add(line)
        print("Number of reads that was corrected is: {}".format(len(unique_reads)))

    @staticmethod
    def convert_proovread_corrected_read_to_python(corrected_file):
        corrected_reads_coordinate = {}
        with open(corrected_file, 'r') as input_correct:
            for line in input_correct:
                if line.startswith('>'):
                    seq_end = 0
                    end = 0
                    line_split = line.replace("SUBSTR:", "").replace('>', '').split()
                    if line_split[1].startswith('HPL'):
                        del(line_split[1])
                    if line_split[-1].startswith('SI'):
                        del(line_split[-1])
                    seq = line_split[0].split(".", 1)[0]
                    if len(line_split) > 2:
                        if "," in line_split[1]:
                            begin, end = line_split[1].split(',')
                        else:
                            begin = line_split[1]
                        # get the substring from last element of the list which contains the coordinate
                        last_element = line_split[-1]
                        if "," in last_element:
                            seq_begin, seq_end = last_element.split(',')
                            seq_begin = int(seq_begin) + int(begin)
                            seq_end = seq_begin + int(seq_end)
                        else:
                            seq_begin = int(last_element) + int(begin)
                            seq_end = seq_begin + int(end)
                    else:
                        if "," in line_split[1]:
                            seq_begin, seq_end = line_split[1].split(',')
                            seq_end = int(seq_begin) + int(seq_end)
                        else:
                            seq_begin = line_split[1]

                    if seq in corrected_reads_coordinate:
                        corrected_reads_coordinate[seq].append([int(seq_begin), seq_end])
                    else:
                        corrected_reads_coordinate[seq] = []
                        corrected_reads_coordinate[seq].append([int(seq_begin), seq_end])
        return AssessUtil.sort_dictionary_values(corrected_reads_coordinate)

    @staticmethod
    def remove_string(target_list, query, replace_value):
        list_length = len(target_list)
        for i in range(0, list_length):
            target_list[i] = target_list[i].replace(query, replace_value)
        return target_list

    @staticmethod
    def sort_dictionary_values(unsorted_dictionary):
        sorted_dictionary = {}
        for k, v in unsorted_dictionary.iteritems():
            sorted_dictionary[k] = sorted(v, key=lambda x: x[0])
        return sorted_dictionary

    @staticmethod
    def uncorrected_interval(corrected_coordinate, original_seq_length=None):
        uncorrected_dict_coordinate = {}
        if original_seq_length:
            for k, v in corrected_coordinate.iteritems():
                uncorrected_dict_coordinate[k] = AssessUtil.list_complement_with_seq_len(v, original_seq_length[k])
        else:
            for k, v in corrected_coordinate.iteritems():
                uncorrected_dict_coordinate[k] = AssessUtil.list_complement(v)
        return uncorrected_dict_coordinate

    @staticmethod
    def list_complement(my_list):
        complement_list = []
        first = True
        last = len(my_list)
        i = 0
        last_element = 0
        for element in my_list:
            i += 1
            if first:
                if int(element[0]) == 0:
                    # complement_list.append([0, element[0]])
                    last_element = element[1]+1
                else:
                    complement_list.append([0, element[0]-1])
                    last_element = element[1]+1
                first = False
            else:
                complement_list.append([last_element, element[0]-1])
                last_element = element[1]+1
            if last == i:
                complement_list.append([last_element, 0])
        if len(my_list) == 1:
            complement_list.append([last_element, 0])
        return complement_list

    @staticmethod
    def list_complement_with_seq_len(my_list, begin_length):
        complement_list = []
        first = True
        last = len(my_list)
        i = 0
        last_element = 0
        for element in my_list:
            i += 1
            if first:
                if int(element[0]) == 0:
                    # complement_list.append([0, element[0]])
                    last_element = element[1]+1
                else:
                    if element[0] > begin_length[0]:
                        complement_list.append([begin_length[0], element[0]-1])
                    last_element = element[1]+1
                first = False
            else:
                if last_element < element[0]-1:
                    complement_list.append([last_element, element[0]-1])
                last_element = element[1]+1
            if last == i:
                if last_element < (int(begin_length[0]) + int(begin_length[1])):
                    complement_list.append([last_element, (int(begin_length[0]) + int(begin_length[1]))])
        if (last == 1) and (last_element < (int(begin_length[0]) + int(begin_length[1]))):
            complement_list.append([last_element, (int(begin_length[0]) + int(begin_length[1]))])
        return complement_list

    @staticmethod
    def extract_proovread_corrected_reads_by_coordinate(full_corrected_proovread_file, read_coordinate,
                                                        uncorrected_reads, threshold):
        my_full_reads = Fasta(full_corrected_proovread_file)
        assessUtil = AssessUtil()
        with open(uncorrected_reads, 'w') as output:
            for k, v in read_coordinate.iteritems():
                    i = 1
                    for element in v:
                        my_new_element = assessUtil.get_element_coordinate(element, len(my_full_reads[k]), threshold)
                        if my_new_element is not None:
                            if element[1] == 0 or element[1] < len(my_full_reads[k]):
                                output.write(">{}.{} | {}\n{}\n".format(k, i, my_full_reads[k][
                                                                            my_new_element[0]:my_new_element[
                                                                                1]].longname,
                                                                      my_full_reads[k][
                                                                      my_new_element[0]:my_new_element[1]]
                                                                      ))
                            i += 1

    @staticmethod
    def extract_uncorrected_reads_proovread(trimmed_corrected_reads, full_corrected_reads, output_file, threshold=50):
        corrected_reads_python_coordinate = AssessUtil.convert_proovread_corrected_read_to_python(
            trimmed_corrected_reads)
        uncorrected_reads_python_coordinate = AssessUtil.uncorrected_interval(corrected_reads_python_coordinate)
        AssessUtil.extract_proovread_corrected_reads_by_coordinate(full_corrected_proovread_file=full_corrected_reads,
                                                                   read_coordinate=uncorrected_reads_python_coordinate,
                                                                   uncorrected_reads=output_file, threshold=threshold)
        print("Done")

    @classmethod
    def get_element_coordinate(cls, element, length_of_read, threshold):
        if element[0] < element[1]:
            if element[1] - element[0] >= threshold:
                return element
            else:
                return None
        elif element[1] == 0:
            if element[0] >= length_of_read:
                return None
            else:
                if element[1] - element[0] >= threshold:
                    return [element[0], length_of_read]
                else:
                    return None

if __name__ == "__main__":
    print(AssessUtil.complexity("AggutttAAAAAcAA"))


