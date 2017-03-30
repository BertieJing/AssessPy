# import

try:
    from pyfaidx import Fasta
    import re
except ImportError as e:
    print(e.message)


# class definition


class ReadFasta:

    def __init__(self):
        pass

    @classmethod
    def read_general_corrected(cls, corrected_file):
        corrected_dictionary = {}
        with Fasta(corrected_file) as data_in:
            for element in data_in:
                corrected_dictionary[element.name] = str(element)
        return corrected_dictionary

    # read lordec corrected fasta file
    @staticmethod
    def read_lordec_result(lordec_file):
        lordec_dictionary = {}
        if lordec_file:
            try:
                with Fasta(lordec_file) as data_in:
                    for element in data_in:
                        seq_length = len(str(element))

                        if element.name.count("_") == 2:
                            sub_seq_name = element.name.replace(">", "").strip()
                            seq_name = element.name.replace(">", "").strip().rsplit("_", 1)[0]
                        else:
                            seq_name = sub_seq_name = element.name.replace(">", "").strip()

                        if len(str(element)) >= 200:
                            if seq_name in lordec_dictionary:
                                # append
                                lordec_dictionary[seq_name].append([seq_length, sub_seq_name, str(element)])
                                lordec_dictionary[seq_name][0][0] += seq_length
                            else:
                                # add
                                lordec_dictionary[seq_name] = []
                                lordec_dictionary[seq_name].append([seq_length])
                                lordec_dictionary[seq_name].append([seq_length, sub_seq_name, str(element)])

            except AttributeError as error:
                print(error)
            # print(len(lordec_dictionary))
            return lordec_dictionary
        else:
            print("file not found ")

    # read Lsc file
    @staticmethod
    def read_lsc_result(lsc_file):
        lsc_dictionary = {}
        with Fasta(lsc_file) as data_in:
            for element in data_in:
                lsc_dictionary[element.name.split("|")[0]] = str(element)
        return lsc_dictionary

    # read proovread corrected file
    @staticmethod
    def read_proovread_result(proovread_file):
        proovread_dictionary = {}
        if proovread_file:
            try:
                with Fasta(proovread_file) as data_in:
                    for element in data_in:
                        seq_length = len(str(element))

                        if element.name.find(".") > -1:
                            sub_seq_name = element.name.replace(">", "").strip()
                            seq_name = element.name.replace(">", "").strip().split(".", 1)[0]
                        else:
                            seq_name = sub_seq_name = element.name.replace(">", "").strip()

                        if len(str(element)) >= 200:
                            if seq_name in proovread_dictionary:
                                # append
                                proovread_dictionary[seq_name].append([seq_length, sub_seq_name, str(element)])
                                proovread_dictionary[seq_name][0][0] += seq_length
                            else:
                                # add
                                proovread_dictionary[seq_name] = []
                                proovread_dictionary[seq_name].append([seq_length])
                                proovread_dictionary[seq_name].append([seq_length, sub_seq_name, str(element)])
                return proovread_dictionary
            except AttributeError as proovread_error:
                print(proovread_error.message)

    # read pbcr file
    def read_pbcr_result(self, pbcr_file, log_file):
        old_new_dic = self.extract_original_new_name(log_file)
        pbcr_fasta_read = {}
        number_of_lost_data = 0
        with Fasta(pbcr_file) as data_in:
            for element in data_in:
                sequence_length = len(str(element))
                if sequence_length >= 200:
                    seq_name = element.name
                    name, new_name, start, end = self.extract_numbers(seq_name)
                    if name in old_new_dic.keys():
                            old_name = old_new_dic[name]

                            if old_name in pbcr_fasta_read:
                                pbcr_fasta_read[old_name].append([sequence_length, new_name, "{}-{}".format(start, end),
                                                                  str(element)])
                                # adding the new length.
                                pbcr_fasta_read[old_name][0][0] += sequence_length
                            else:
                                pbcr_fasta_read[old_name] = []
                                pbcr_fasta_read[old_name].append([sequence_length])
                                pbcr_fasta_read[old_name].append([sequence_length, new_name, "{}-{}".format(start, end),
                                                                  str(element)])
                    else:
                        number_of_lost_data += 1
        return pbcr_fasta_read, number_of_lost_data

    def extract_original_new_name(self, input_file):
        if input_file:
            mapping_new_to_original = {}
            try:
                with open(input_file, 'r') as input_data:
                    # skipping the title line
                    input_data.next()
                    for line in input_data:
                        line = line.strip()
                        line_split = line.split()
                        original_name = line_split[0]
                        new_name = self.extract_numbers_pbcr(line_split[1])
                        mapping_new_to_original[new_name] = original_name
                return mapping_new_to_original
            except IOError as details:
                print(details)

    def extract_numbers_pbcr(self, string):
        if string:
            return re.findall('\d+', string)[0]

    def extract_numbers(slef, string):
        # I am using negative number to get coordinate so if the name of the file contains number wont effect
        #  the end result.
        if string:
            numbers = (re.findall('\d+', string))
            return numbers[-4], "{}_{}".format(numbers[-4], numbers[-3]), numbers[-2], numbers[-1]

if __name__ == "__main__":
    s = ReadFasta()
    # d = s.read_proovread_result("/data/test_data_from_server_maize/proovread/ecoli/ecoli-backup/ecoli.trimmed.fa")
    d = s.read_lsc_result("/data/test_data_from_server_maize/LSC/ecoli/output/corrected_LR.fa")
    print(len(d))

    print(d["S1_11366|0.75"])


