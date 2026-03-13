#################################
##                             ##
## Aligns initial fasta files  ##
##                             ##
#################################

import subprocess
import argparse
import pathlib
import sys


def get_align_call(args):
    star_call  = "STAR --runThreadN " + str(args.threads) + " "
    star_call += "--genomeDir " + str(args.star_idx) + " "
    star_call += "--readFilesIn " + str(args.read1) + " " + str(args.read2) + " "
    star_call += "--outFilterMultimapNmax 100 "
    star_call += "--winAnchorMultimapNmax 100 "
    star_call += "--outFileNamePrefix " + str(args.out) + " " 
    star_call += " --outFilterScoreMinOverLread 0.3 "
    star_call += "--outFilterMatchNminOverLread 0.3 "
    star_call += "--outFilterMatchNmin 15 "
    star_call += "> " + str(args.out) + ".progress.log  2>&1" 
    return star_call


def clean(loc):
    old_sam = str(loc) + "Aligned.out.sam"
    new_sam = str(loc) + ".sam"
    old_log = str(loc) + "Log.final.out"
    new_log = str(loc) + ".final.log"
    mv_calls = ["mv " + old_sam + " " + new_sam, "mv " + old_log + " " + new_log]
    rm_files = ["Log.out", "Log.progress.out", "SJ.out.tab"]
    rm_calls = ["rm " + str(loc) + f for f in rm_files]

    for call in mv_calls:
        subprocess.call(call, shell=True)
    for call in rm_calls:
        subprocess.call(call, shell=True)


################
## driver fxn ##
################

if __name__ == '__main__':
    ## Set up argument parser ##
    parser = argparse.ArgumentParser(prog="TE-SAEM align",
             description="Alignment wrapper for TE-SAEM utilizing the STAR aligner",
             epilog="For help//queries please visit our github: https://github.com/synanth/TE-SAEM")

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
                        help = "Location of read1 (default: %(default)s)")
    parser.add_argument("-t",
                        type = int,
                        default = 1,
                        metavar = "threads",
                        dest = "threads",
                        help = "Number of threads to use (default: %(default)s)")
    parser.add_argument("-o",
                        type = str,
                        required = True,
                        metavar = "out",
                        dest = "out",
                        help = "Output SAM location (default: %(default)s)")
    parser.add_argument("-s",
                        type = pathlib.Path,
                        required = True,
                        metavar = "star_idx",
                        dest = "star_idx",
                        help = "Location of star index (default: %(default)s)")

    args = parser.parse_args()    
    
    # Run align call ##
    align_call = get_align_call(args)
    print(align_call + "\n")
    subprocess.run(align_call, shell=True)

    ## Clean file structure ##
    print("Cleaning align file structures")
    clean(args.out)
