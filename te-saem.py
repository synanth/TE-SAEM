##########################################################
#                                                        #
#  TE-SAEM                                               #
#    - A program to quantify TEs from bulk RNA-Seq data  #
#                                                        #
#  S. Texocaelum                                         #
#  Dawei Li Lab @ TTUHSC                                 #
#                                                        #
##########################################################

#!/bin/bash

import os
import pathlib
import argparse
import subprocess


## Load .ini ##
def load_ini(args, program_path):
    gtf_loc, genome_loc = None, None

    with open(args.ini, "r") as f:
        for line in f.readlines():
            buff = line.strip().split()
            if buff[0] == "genome_idx":
                if buff[-1][1:3] == "./":
                    genome_loc = program_path + buff[-1][3:-1]
            elif buff[0] == "gtf":
                if buff[-1][1:3] == "./":
                    gtf_loc = program_path + buff[-1][3:-1]
    if any([gtf_loc, genome_loc]) == None:
        print("Parameter .ini not loaded, please check its format")
        exit()
    return gtf_loc, genome_loc


def get_working_dir(args):
    working_dir = "/".join(str(args.out).split("/")[:-1]) + "/"
    if working_dir[0] != "/":
        working_dir = "./" + working_dir
    return working_dir
        

## Check existence of arguments ##
def check_ref(ref_path, ref, args):
    if ref == "genome_fa":
        error_str = "Genome fasta "
        file_path = ref_path + args.genome + ".fa"
    elif ref == "te_gtf":
        error_str = "Transposable element GTF "
        file_path = ref_path + args.genome + ".gtf"
    elif ref == "te_len":
        error_str = "Transposable element lengths file "
        file_path = ref_path + args.genome + ".len"
    elif ref == "idx":
        error_str = "Star index "
        file_path = ref_path + "star_" + args.genome + "_idx/"
    if not os.path.exists(file_path):
        print(error_str + "does not exist in references folder.")
        print("Please run setup to retrieve and parse genome information.")
        exit()


def check_read(read_path):
    if not os.path.exists(read_path):
        print(str(read_path) + " does not exist.")
        exit()


## Driver function ##
if __name__ == '__main__':
    ## Clear screen ##
    os.system("clear")

    ## Set up paths ##
    program_path = "/".join(os.path.abspath(__file__).split("/")[:-1]) + "/"
    scripts_path = program_path + "scripts/"
    

    ## Set up argument parser ##
    parser = argparse.ArgumentParser(prog="TE-SAEM",
             description="Quantification of TE expression from bulk RNA-Seq data through the utilization of a simulated annealing based expectation maximization algorithm",
             epilog="For help//queries please visit our github: https://github.com/synanth/TE-SAEM")


    parser.add_argument("-i", 
                        type = pathlib.Path,
                        default = program_path + "params.ini",
                        metavar = "ini",
                        dest = "ini",
                        help = ".ini file location (default: %(default)s)")
    parser.add_argument("-1", 
                        type = pathlib.Path,
                        required = True,
                        metavar = "read1",
                        dest = "read1",
                        help = "Location of read1 (default: %(default)s)")
    parser.add_argument("-2", 
                        type = pathlib.Path,
                        required = True,
                        metavar = "read2",
                        dest = "read2",
                        help = "Location of read2 (default: %(default)s)")
    parser.add_argument("-o",
                        type = pathlib.Path,
                        required = True,
                        default = "te-counts.csv",
                        metavar = "out",
                        dest = "out",
                        help = "Output count table location (default: %(default)s)")
    parser.add_argument("-c",
                        type = bool,
                        default = True,
                        metavar = "clean",
                        dest = "clean",
                        help = "Clean intermediate files (default: %(default)s)")
    parser.add_argument("-t",
                        type = int,
                        default = 1,
                        metavar = "threads",
                        dest = "threads",
                        help = "Number of threads to use (default: %(default)s)")


    ## Check parameters ##
    args = parser.parse_args()
    gtf_loc, genome_loc = load_ini(args, program_path)
    working_dir = get_working_dir(args)
    check_read(args.read1)
    check_read(args.read2)

    ## Check if initial setup has been performed ##

    ## Alignment ##
    align_call  = "python3 " + scripts_path + "align.py "
    align_call += "-1 " + str(args.read1) + " "
    align_call += "-2 " + str(args.read2) + " "
    align_call += "-t " + str(args.threads) + " "
    align_call += "-o " + working_dir + "star "
    align_call += "-s " + genome_loc

    print("Running external STAR command with the following parameters:\n")
    subprocess.call(align_call, shell=True)
    
    ## Parse post align files ##
    parse_call  = "python3 " + scripts_path + "parse_alignment.py "
    parse_call += "-s " + working_dir + "star.sam "
    parse_call += "-t " + str(args.threads) + " "
    parse_call += "-g " + gtf_loc

    print("Parsing alignment")
    subprocess.call(parse_call, shell=True)
 
    ## Run SAEM ##
    saem_call  = "python3 " + scripts_path + "saem.py "
    saem_call += "-d " + working_dir + " "
    saem_call += "-g " + gtf_loc + " " 
    saem_call += "-o " + str(args.out)

    print("Running SAEM")
    subprocess.call(saem_call, shell=True)
