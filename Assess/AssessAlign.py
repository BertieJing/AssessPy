# import
try:
    from Bio.Align.Applications import MafftCommandline
    import fileinput, os, subprocess
    from StringIO import StringIO
    from Bio import AlignIO
    from datetime import datetime
    from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline
    import re
except ImportError as import_error:
    print("Problem in importing data see message --> \n{}\n".format(import_error.message))
# class


class AssessAlign:

    def __init__(self):
        pass

    @staticmethod
    def mafft_align(query_seq, target_seq, query_name, target_name, align_method="local", directory="./", quiet=False):
        # add time to file name to make it unique '20160809-144522_' 2016-08-09 14:45:22
        file_name = directory + datetime.now().strftime("%Y%m%d_%H%M%S_") + query_name + ".fasta"
        with open(file_name, 'w') as data_out:
            data_out.write(">{}\n{}\n>{}\n{}".format(target_name, target_seq, query_name, query_seq))
        if align_method == "local":
            mafft_cline = MafftCommandline(input=directory + file_name, nuc=True, localpair=True, maxiterate=1000,
                                            quiet=quiet)

        else:
            mafft_cline = MafftCommandline(input=directory + file_name, nuc=True, globalpair=True, maxiterate=1000,
                                           quiet=quiet)
        out, _ = mafft_cline()
        align = AlignIO.read(StringIO(out), "fasta")
        my_list = list(align)
        # target name, target seq, query name, query seq
        os.remove(file_name)
        return my_list[0].id, my_list[0].seq, my_list[1].id, my_list[1].seq

    @staticmethod
    def mafft_align_command_line(query_seq, target_seq, query_name, target_name, align_method="local",
                                 directory="./"):
        align_file = directory + datetime.now().strftime("%Y%m%d_%H%M%S_") + query_name + ".fasta"
        with open(align_file, 'w') as query_file:
            query_file.write(">{}\n{}\n>{}\n{}".format(target_name, target_seq, query_name, query_seq))
        if align_method == "local":
            mafft_align = subprocess.Popen(["mafft", "--treeout", "--thread", "-1", "--quiet", "--localpair", "--maxiterate", "1000"
                                               ,align_file], stdout=subprocess.PIPE)
        else:
            mafft_align = subprocess.Popen(["mafft", "--treeout", "--thread", "-1", "--quiet", "--globalpair", "--maxiterate", "1000"
                                               , align_file], stdout=subprocess.PIPE)
        mafft_align.communicate()
        mafft_align.wait()
        mafft_align.stdout.close()
        os.remove(align_file)
        # return my_list[0].id, my_list[0].seq, my_list[1].id, my_list[1].seq
        with open(align_file+".tree", 'r') as data_in:
            result = data_in.readlines()[2]
        os.remove(align_file+".tree")
        return float(re.findall(r'\d+\.\d+', result)[0])

    @staticmethod
    def fogsaa_align_command_line(query_seq, target_seq, query_name, target_name, directory="./", dna="1", gap="1",
                                  match="1", mismatch="-1", gap_opening="-10", gap_extension="-1"):
        query_file = directory + datetime.now().strftime("%Y%m%d_%H%M%S_") + query_name + ".txt"
        target_file = directory + datetime.now().strftime("%Y%m%d_%H%M%S_") + query_name + "_ref.txt"

        query_file_al = directory + datetime.now().strftime("%Y%m%d_%H%M%S_") + query_name + "_al.txt"
        target_file_al = directory + datetime.now().strftime("%Y%m%d_%H%M%S_") + query_name + "_al_ref.txt"

        with open(query_file, 'w') as q_file, open(target_file, 'w') as t_file:
            q_file.write(">{}\n{}".format(query_name, query_seq))
            t_file.write(">{}\n{}".format(target_name, target_seq))
        # print("****\n my code \n***")
        print("fogsaa_twik", query_file, target_file, dna, gap, match, mismatch,
                                         gap_opening, gap_extension, query_file_al, target_file_al)
        fogsaa_align = subprocess.Popen(["/home/medhat/source/FOGSAA/fogsaa2", query_file, target_file, dna, gap,
                                         match, mismatch, gap_opening, gap_extension, query_file_al, target_file_al],
                                        stdout=subprocess.PIPE)

        out, _ = fogsaa_align.communicate()
        print(out)
        fogsaa_align.wait()
        fogsaa_align.stdout.close()
        q_align = ""
        t_align = ""
        print(open(query_file_al).readlines())
        print(query_file_al, os.path.exists(query_file_al), os.path.isfile(query_file_al))
        if os.path.isfile(query_file_al):
            try:

                with open(query_file_al, 'r') as query_align_data, open(target_file_al, 'r') as target_align_data:
                    q_align = query_align_data.readlines()
                    t_align = target_align_data.readlines()

                os.remove(query_file)
                os.remove(query_file_al)
                os.remove(target_file)
                os.remove(target_file_al)
                if os.path.isfile(directory + "seq1.txt"):
                    os.remove(directory + "seq1.txt")
                    os.remove(directory + "seq2.txt")
            except Exception as e:
                print(e.message)
        else:
            print("this sequence {} was not aligned".format(query_name))
        return target_name, t_align, query_name, q_align

    @staticmethod
    def align_score(seq_one, seq_two):
        matches = 0
        mismatches = 0
        gap = 0
        for i in range(min(len(seq_one), len(seq_two))):
            if (seq_one[i].lower() == seq_two[i].lower()) and (seq_one[i] != "-" and seq_two[i] != "-"):
                matches += 1
            elif seq_one[i] == "-" or seq_two[i] == "-":
                gap += 1
            else:
                mismatches += 1
        return matches, mismatches, gap

    @staticmethod
    def needle_align_code(query_seq, target_seq):
        needle_cline = NeedleCommandline(asequence="asis:" + query_seq,
                                         bsequence="asis:" + target_seq,
                                         aformat="simple",
                                         gapopen=10,
                                         gapextend=0.5,
                                         outfile='stdout'
                                         )
        out_data, err = needle_cline()
        # print("data from stdout {}".format(out_data))
        out_split = out_data.split("\n")
        p = re.compile("\((.*)\)")
        return p.search(out_split[25]).group(1).replace("%", "")

    @staticmethod
    def water_align_code(query_seq, target_seq):
        needle_cline = WaterCommandline(asequence="asis:" + query_seq,
                                        bsequence="asis:" + target_seq,
                                        aformat="simple",
                                        gapopen=10,
                                        gapextend=0.5,
                                        outfile='stdout'
                                        )
        out_data, err = needle_cline()
        out_split = out_data.split("\n")
        p = re.compile("\((.*)\)")
        return p.search(out_split[25]).group(1).replace("%", "")

    @staticmethod
    def water_align_command(query_seq, target_seq):
        water_cline = subprocess.Popen(["/home/medhat/source/EMBOSS-test/EMBOSS-6.6.0/emboss/water",
                                        "-asequence=asis:" + query_seq,
                                        "-bsequence=asis:" + target_seq,
                                        "-aformat", "srspair",
                                        "-gapopen", "10",
                                        "-gapextend", "0.5",
                                        "-outfile", "stdout"
                                        ], stdout=subprocess.PIPE)

        out_data, _ = water_cline.communicate()
        water_cline.wait()
        # print(water_cline.returncode)
        water_cline.stdout.close()
        # print(out_data)
        out_split = out_data.split("\n")
        p = re.compile("\((.*)\)")
        return p.search(out_split[25]).group(1).replace("%", "")


if __name__ == "__main__":
    print("Put some code")




