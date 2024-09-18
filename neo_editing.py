import argparse
import subprocess
import pandas as pd
import os
from Bio.Seq import Seq

# Global variables:
SIZE_OF_SEQ = 88
MER_LENGHT = 15
VERSION = "diabetes"
FASTA_PATH = None
NETMHCPAN_PATH = None

def setting_Environment_Variables():
    """
    This method configures the environment variables to enable the subsequent execution of netMHCIIpan4.3.
    """
    os.environ["NETMHCIIpan"] = NETMHCPAN_PATH
    os.environ["PLATFORM"] = "Linux_x86_64"
    os.environ["NMHOME"] = NETMHCPAN_PATH
    os.environ["NetMHCIIpanPLAT"] = os.environ["NMHOME"] + "/" + os.environ["PLATFORM"]
    os.environ["NetMHCIIpanWWWPATH"] = "/services/NetMHCIIpan-4.3/tmp"
    os.environ["NetMHCIIpanWWWDIR"] = "/var/www/html/services/NetMHCIIpan-4.3/tmp"
    os.environ["TMPDIR"] = "/tmp"

def delete_files(folder_path):
    # iterate over all files in the folder
    for filename in os.listdir(folder_path):
        # construct the full file path
        file_path = os.path.join(folder_path, filename)
        
        # check if the file path is a regular file (not a folder or a symlink)
        if os.path.isfile(file_path):
            # delete the file
            os.remove(file_path)

def creat_FASTA_file(temp_path, list_of_seq, num):
    with open(os.path.join(temp_path, f"input_seq_{num}.fsa"), "w") as fasta_file:
        for i, sequence in enumerate(list_of_seq):
            fasta_file.write(f">{i}\n{str(Seq(sequence).translate())}\n")

def run_MHCpan2(temp_path, hla, num):
    output = subprocess.run(f'{os.path.join(NETMHCPAN_PATH,"Linux_x86_64", "bin", "NetMHCIIpan-4.3")} -f {os.path.join(temp_path, f"input_seq_{num}.fsa")} -a {hla} -length {MER_LENGHT}', shell=True, capture_output=True, text=True)
    return output.stdout

def create_results_file(sb_tuple_edit, wb_tuple_edit, output_edit, output_NO_edit=None, sb_tuple_NO_edit=None, wb_tuple_NO_edit=None):
    output_edit_list = output_edit.split("\n")
    output_NO_edit_list = output_NO_edit.split("\n")

    for result in output_edit_list:
        if result[-2:] == "WB":
            final_result = [value for value in result.split(" ") if value.strip()]
            wb_tuple_edit.append((int(final_result[6]), final_result[2], float(final_result[9]), final_result[1], final_result[0]))
        if result[-2:] == "SB":
            final_result = [value for value in result.split(" ") if value.strip()]
            sb_tuple_edit.append((int(final_result[6]), final_result[2], float(final_result[9]), final_result[1], final_result[0])) 

    for result in output_NO_edit_list:
        if result[-2:] == "WB":
            final_result = [value for value in result.split(" ") if value.strip()]
            wb_tuple_NO_edit.append((int(final_result[6]), final_result[1], final_result[0]))
        if result[-2:] == "SB":
            final_result = [value for value in result.split(" ") if value.strip()]
            sb_tuple_NO_edit.append((int(final_result[6]), final_result[1], final_result[0])) 

def main(args):
    if VERSION not in ["most_common", "diabetes"]:
        exit(0)
    list_of_HLA_most_common = ['HLA-DQA10602-DQB10602', 'HLA-DPA10201-DPB10201', 'DRB1_1501', 'HLA-DPA10401-DPB10401', 'HLA-DQA10302-DQB10302']
    list_of_HLA_diabetes = ['HLA-DQA10302-DQB10302', 'HLA-DQA10201-DQB10201', 'DRB1_0401', 'DRB1_0301']
    check_flag = True
    hla_dic = {"diabetes": list_of_HLA_diabetes, "most_common": list_of_HLA_most_common}
    setting_Environment_Variables()
    project_dir = args.project_dir
    results_dir = args.results_dir
    sup_dir = args.sup_dir
    os.makedirs(os.path.join(results_dir , VERSION), exist_ok=True)
    if os.path.exists(os.path.join(results_dir , VERSION, "sorted_by_HLA")):
        delete_files(os.path.join(results_dir , VERSION, "sorted_by_HLA"))
    else:
        os.mkdir(os.path.join(results_dir , VERSION, "sorted_by_HLA"))
    with open (os.path.join(results_dir , VERSION, "final_results.tsv"), 'w') as final_file:
        final_file.write("#NEO-ADAR-TIGEN 1.0\n")
        final_file.write(f"Gene_Name\tEdit-site\tNeoantigen\n")
        editing_file = pd.read_csv(os.path.join(project_dir , "NEW_orsh_sites_nature_AG.bed6.bed"), sep="\t", header=None)
        for _, row in editing_file.iterrows():
            mut_chr = row[0]
            mut_pos = row[2]
            mut_seq = row[3][1]
            mutation = f"{row[0]}:g.{row[2]}{row[3][0]}>{row[3][1]}"

            # Run a program that accepts as input a chromosome and location and returns the 60 bases around it:
            p = subprocess.run(f'sh {os.path.join(os.getcwd(), "src", "find_seq_of_mutation.sh")} {mut_chr} {mut_pos} {int(SIZE_OF_SEQ / 2)} {FASTA_PATH}', capture_output=True,text=True, shell=True)
            seq = p.stdout.replace("\n", "") 

            # We will run a program that will return the reading frame of the sequence:
            q = subprocess.run(f'sh {os.path.join(os.getcwd(), "src", "find_position.sh")} {mut_chr} {mut_pos} {os.path.join(sup_dir, "bed6_of_ensemble.bed")}', shell=True, capture_output=True,text=True)
            res_list = [res for res in q.stdout.split("\n") if res]

            if not res_list:
                final_file.write(f"{row[4]}\t{mutation}\tNo1\n")
                continue
            
            temp = seq[::-1]
            list_of_codons = []
            frame_init_pos = int(mut_pos) - int((SIZE_OF_SEQ/2) - 1)

            for res in res_list:
                row_s = [value for value in res.split("\t") if value.strip()] 
                strand = row_s[4]
                reading_frame = int(row_s[5])
                if strand == "+":
                    cds_start_pos = int(row_s[1])
                    cds_end_pos = int(row_s[2])
                    count = 0
                    count2 = 1
                else:
                    count = 1
                    count2 = 0
                    seq = temp
                    frame_init_pos = int(mut_pos) + int((SIZE_OF_SEQ/2) - 1)
                    cds_start_pos = int(row_s[2])
                    cds_end_pos = int(row_s[1])
                end = len(seq)
                flag = True
                if abs(int(mut_pos) - cds_start_pos) < int((SIZE_OF_SEQ/2) - 1): 
                    start = int((SIZE_OF_SEQ/2) - 1) - abs(int(mut_pos) - cds_start_pos) + reading_frame + count
                    frame_init_pos = cds_start_pos + reading_frame 
                    flag = False
                if abs(cds_end_pos - int(mut_pos)) < int(SIZE_OF_SEQ/2):
                    end = abs(cds_end_pos - frame_init_pos) 
                    if not flag:
                        end += start
                if flag:
                    if (reading_frame == 0):
                        seq_frame = abs(frame_init_pos - cds_start_pos + count) % 3 
                    elif(reading_frame == 1):
                        seq_frame = abs(frame_init_pos - 1 - cds_start_pos) % 3 
                    else:
                        seq_frame = abs(frame_init_pos - 2 - cds_start_pos - count) % 3

                    if seq_frame == 0:
                        start = seq_frame 
                    elif seq_frame == 1:
                        start = seq_frame + 1
                    elif seq_frame == 2:
                        start = seq_frame - 1

                codon = [seq[i:i+3] for i in range(start, end, 3)]
                # Remove duplicates:
                if codon not in list_of_codons:
                    list_of_codons.append(codon)

            # Creats mutatation in the seq for each isoform:
            mut_isoforms = []
            pos = int(SIZE_OF_SEQ/2) - start - count2 + count
            if strand == "-":
                mut_seq = str(Seq(mut_seq).complement())
                pos = pos - len(mut_seq)
            for codon in list_of_codons:
                # Convert the list of codons to a DNA sequence
                dna_seq = ''.join(codon)
                edit_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos + 1:]
                if strand == "-":
                    dna_seq = str(Seq(dna_seq).complement())
                    edit_dna_seq = str(Seq(edit_dna_seq).complement())
                if len(edit_dna_seq) < MER_LENGHT * 3:
                    continue
                mut_isoforms.append((edit_dna_seq, dna_seq))

            if mut_isoforms:
                sb_tuple_edit = []
                wb_tuple_edit = []
                creat_FASTA_file(sup_dir, [m[0] for m in mut_isoforms], 1)
                if check_flag:
                    sb_tuple_NO_edit = []
                    wb_tuple_NO_edit = []
                    creat_FASTA_file(sup_dir, [m[1] for m in mut_isoforms], 2)
                for hla in hla_dic[VERSION]:
                    output1 = run_MHCpan2(sup_dir, hla, 1)
                    if check_flag:
                        # Running the tool on the non-edited sequences for comparison:
                        output2 = run_MHCpan2(sup_dir, hla, 2)
                        create_results_file(sb_tuple_edit, wb_tuple_edit, output1, output2, sb_tuple_NO_edit, wb_tuple_NO_edit) 
                    else:
                        create_results_file(sb_tuple_edit, wb_tuple_edit, output1) 


                uniq_sb = list(set(sb_tuple_edit))
                uniq_wb = list(set(wb_tuple_edit))

                if check_flag:
                    wb_sb_NO_edit = list(set(sb_tuple_NO_edit)) + list(set(wb_tuple_NO_edit))

                if not uniq_sb and not uniq_wb:
                    final_file.write(f"{row[4]}\t{mutation}\tNo2\n")
                    continue

                first_flag = True

                if uniq_sb:
                    for sb in uniq_sb:
                        if check_flag and (sb[0],sb[3],sb[4]) in wb_sb_NO_edit:
                            continue
                        with open (os.path.join(results_dir , VERSION, "sorted_by_HLA", f"{sb[3]}.tsv"), 'a') as final:
                            final.write(f"{row[4]}\t{mutation}\t{strand}\tSB\t{sb[2]}\t{sb[1]}\t{row[6]}\n")
                        if first_flag:
                            first_flag = False
                            final_file.write(f"{row[4]}\t{mutation}\tSB\n")

                if uniq_wb:
                    for wb in uniq_wb:
                        if check_flag and (wb[0],wb[3],wb[4]) in wb_sb_NO_edit:
                            continue
                        with open (os.path.join(results_dir , VERSION, "sorted_by_HLA", f"{wb[3]}.tsv"), 'a') as final:
                            final.write(f"{row[4]}\t{mutation}\t{strand}\tWB\t{wb[2]}\t{wb[1]}\t{row[6]}\n")
                        if first_flag:
                            final_file.write(f"{row[4]}\t{mutation}\tWB\n")
                            first_flag = False

                if first_flag: 
                    final_file.write(f"{row[4]}\t{mutation}\tNo3\n")

if __name__ == "__main__":
    # Create parser
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='A tool for creating a neo-antigen by RNA editing')
    parser.add_argument('-i', '--project_dir', dest='project_dir', action='store', metavar='root_dir',
                        nargs='?', default=os.path.join(os.getcwd(), "testdata"), help='A folder containing the projects on which the tool will run')
    parser.add_argument('-o', '--results_dir', dest='results_dir', action='store', metavar='output_file', nargs='?',
                        default=os.path.join(os.getcwd(), "results"), help='The folder that will contain the result files')
    parser.add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path', nargs='?',
                        default=None, help='Log file, default is to create in input dir')
    parser.add_argument('-t','--sup_dir', dest='sup_dir', action='store' , metavar='sup_dir', nargs='?',
                        default=os.path.join(os.getcwd(), "sup"), help='A folder containing completions for the tool code')
    parser.add_argument('-p','--netmhc_path', dest='netmhc_path', action='store' , metavar='netmhc_path', required=True,
                        help='The path to the netMHCIIpan4.3 tool that needs to be downloaded locally')
    parser.add_argument('-f','--hg38_fa', dest='hg38_fa', action='store' , metavar='hg38_fa', required=True,
                         help='The path to the hg38 genom reference that needs to be downloaded locally')
    
    FASTA_PATH = parser.parse_args().hg38_fa
    NETMHCPAN_PATH = parser.parse_args().netmhc_path

    main(parser.parse_args())