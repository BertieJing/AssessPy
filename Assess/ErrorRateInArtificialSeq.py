# import
from AssessUtil import AssessUtil

# code


class ErrorRate:
    def __init__(self, maf_file, blat_alignment_file, corrected_error_rate, uncorrected_error_rate):
        self.corrected_error_rate = corrected_error_rate
        self.uncorrected_error_rate = uncorrected_error_rate
        self.maf_file = maf_file
        self.maf_dictionary = AssessUtil.read_maf_align(self.maf_file)
        self.blat_align = blat_alignment_file

    def get_correction_coordinate(self):
        corrected_read_dictionary = {}
        target_length = {}
        with open(self.blat_align, 'r') as input_data:
            for line in input_data:
                line_split = line.split()
                cor_seq_name = line_split[9]
                original_seq_name = line_split[13]
                if cor_seq_name.split('.', 1)[0] == original_seq_name:
                    begin = int(line_split[15])
                    end = int(line_split[16])
                    if original_seq_name in corrected_read_dictionary:
                        corrected_read_dictionary[original_seq_name].append([begin, end])
                    else:
                        corrected_read_dictionary[original_seq_name] = []
                        corrected_read_dictionary[original_seq_name].append([begin, end])
                        target_length[original_seq_name] = [0, int(line_split[14])]

        return AssessUtil.sort_dictionary_values(unsorted_dictionary=corrected_read_dictionary), target_length

    @classmethod
    def get_uncorrected_coordinate(cls, corrected_coordinate_dict, original_coordinate_dict):
        return AssessUtil.uncorrected_interval(corrected_coordinate=corrected_coordinate_dict, original_seq_length=original_coordinate_dict)

    def get_error_comparison(self):
        corrected_coordination, original_coordination = self.get_correction_coordinate()
        uncorrected_coordinate = self.get_uncorrected_coordinate(corrected_coordinate_dict=corrected_coordination,
                                                                 original_coordinate_dict=original_coordination)
        with open(self.uncorrected_error_rate, 'w') as uncorrected_error_out:
            for seq_name, seq_coordinate_list in uncorrected_coordinate.iteritems():
                original_seq = self.maf_dictionary[seq_name][0]
                corrected_seq_coordinate_list = self.get_original_coordinate(seq_coordinate_list, original_seq)
                for element in corrected_seq_coordinate_list:
                    start = element[0]
                    end = element[1]
                    sub_original_seq = self.maf_dictionary[seq_name][0][start:end]
                    sub_artificial_seq = self.maf_dictionary[seq_name][1][start:end]
                    # print("art -> {}\nref -> {}\n".format(sub_artificial_seq, sub_original_seq))
                    mismatch, indel = self.get_align_value(sub_original_seq, sub_artificial_seq)
                    uncorrected_error_out.write("{}_{}:{}\t{}\t{}\t{}\n".format(seq_name, start, end, (end-start), mismatch, indel))
        with open(self.corrected_error_rate, 'w') as corrected_error_out:
            for seq_name, seq_coordinate_list in corrected_coordination.iteritems():
                original_seq = self.maf_dictionary[seq_name][0]
                corrected_seq_coordinate_list = self.get_original_coordinate(seq_coordinate_list, original_seq)
                for element in corrected_seq_coordinate_list:
                    start = element[0]
                    end = element[1]
                    sub_original_seq = self.maf_dictionary[seq_name][0][start:end]
                    sub_artificial_seq = self.maf_dictionary[seq_name][1][start:end]
                    # print("art -> {}\nref -> {}\n".format(sub_artificial_seq, sub_original_seq))
                    mismatch, indel = self.get_align_value(sub_original_seq, sub_artificial_seq)
                    corrected_error_out.write("{}_{}:{}\t{}\t{}\t{}\n".format(seq_name, start, end, (end-start), mismatch, indel))


    @classmethod
    def get_align_value(cls, query, target):
        indel = mismatch = 0
        for q, t in zip(query, target):
            if q != "-" and t != "-":
                if q != t: mismatch += 1
            else:
                indel += 1
        return mismatch, indel

    @classmethod
    def get_original_coordinate(cls, coordinate_list, original_seq):
        corrected_coordinate_list = []
        for element in coordinate_list:
            element_start = element[0]
            element_end = element[1] - element_start
            i = begin = 0
            for character in original_seq:
                if character.isalpha():
                    i += 1
                begin += 1
                if i == element_start:
                    break
            corrected_coordinate_list.append([begin, begin+element_end])
        return corrected_coordinate_list


if __name__ == "__main__":
    a = ErrorRate("/home/medhat/Downloads/yeast4_pacbio.fa_0001.maf",
                  "/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/yeast_blat_result.psl"
                  , "/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/corrected_error_rate.txt",
                  "/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_error_rate.txt")
    # corrected_coordinate, original_coordinate = a.get_correction_coordinate()
    # uncorrected_coordinate = a.get_uncorrected_coordinate(corrected_coordinate_dict=corrected_coordinate, original_coordinate_dict=original_coordinate)
    a.get_error_comparison()