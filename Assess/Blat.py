# import
import subprocess
import os
import uuid
from AssessUtil import AssessUtil


# code


class Blat:
    def __init__(self, db, query, original_coordinate, align_type=None, output=None, header=True, **kwargs):
        self.db = db
        self.query = query
        self.original_coordinate = original_coordinate
        if output:
            self.output = output
        else:
            unique_filename = uuid.uuid4()
            self.output = os.path.dirname(self.query) + "/" + str(unique_filename) + ".result"
        if align_type:
            self.output_type = align_type
        else:
            self.output_type = "blast8"
        if header:
            self.header = ""
        else:
            self.header = "-noHead"

        self.step = 12
        self.repeat = 256
        self.min_score = 0
        self.min_identity = 0
        for key in kwargs:
            if key == "step":
                self.step = kwargs[key]
            elif key == "repeat":
                self.repeat = kwargs[key]
            elif key == "min_score":
                self.min_score = kwargs[key]
            elif key == "min_identity":
                self.min_identity = kwargs[key]

    def run_blat(self):
        command = "blat {} {} -stepSize={} -repMatch={} -minScore={} -minIdentity={} -out={} {} {}".format(self.db,
                                                                                                           self.query,
                                                                                                           self.step,
                                                                                                           self.repeat,
                                                                                                           self.min_score,
                                                                                                           self.min_identity,
                                                                                                           self.output_type,
                                                                                                           self.header,
                                                                                                           self.output)
        try:
            blat_align = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        except Exception as e:
            print(e.message)
        self.stdout, self.stderr = blat_align.communicate()
        blat_align.wait()
        blat_align.stdout.close()
        print(self.stdout)
        if self.stderr:
            raise IOError('blat failed')
        return True

    def get_corrected_read_coordinate(self):
        corrected_read_dictionary = {}
        print(os.path.isfile(self.output))
        with open(self.output, 'r') as input_data:
            for line in input_data:
                line_split = line.split()
                cor_seq_name = line_split[9]
                original_seq_name = line_split[13]
                begin = line_split[15]
                end = line_split[16]
                if cor_seq_name.split('.', 1)[0] == original_seq_name:
                    begin += self.original_coordinate[original_seq_name][0]
                    end += self.original_coordinate[original_seq_name][0]
                if original_seq_name in corrected_read_dictionary:
                    corrected_read_dictionary[original_seq_name].append([begin, end])
                else:
                    corrected_read_dictionary[original_seq_name] = []
                    corrected_read_dictionary[original_seq_name].append([begin, end])

        return AssessUtil.sort_dictionary_values(unsorted_dictionary=corrected_read_dictionary)

# psl file  http://www.ensembl.org/info/website/upload/psl.html
# text = Query id, Subject id, % identity, alignment length, mismatches,gap openings,
#  q. start, q. end, s. start, s. end, e-value, bit score

# csv = "Query id","Subject id","% identity","alignment length","mismatches",
# "gap openings","q. start","q. end","s. start","s. end","e-value","bit score"
# tabs = Query id Subject id % identity alignment length mismatches gap openings
# q. start q. end s. start s. end e-value bit score

if __name__ == "__main__":
    a = Blat(db="/home/medhat/Downloads/test.fa", query="/home/medhat/Downloads/query.fa", original_coordinate="",
             align_type="psl", header=False)
    a.run_blat()

    # Blat.test_parse_out("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/test_coordent/out.psl")
