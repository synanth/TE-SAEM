#################################################################
##                                                             ##
## This program parses STAR alignment to prepare data for SAEM ##
##                                                             ##
#################################################################

import subprocess
import argparse
import pathlib
import sys
import os


def sort_sam(base_loc, args):
    sam_loc = base_loc + "star.sam"
    new_loc = base_loc + "sorted.sam"
    
    sort_call  = "samtools sort " + sam_loc + " "
    sort_call += "-o " + new_loc + " "
    sort_call += "-@ " + str(args.threads)

    subprocess.call(sort_call, shell=True)


def sam_to_bed(base_loc):
    sam_loc = base_loc + "sorted.sam"
    bed_loc = base_loc + "aligned.bed"
    len_loc = base_loc + "read_lens.txt"
    seq_loc = base_loc + "read_seqs.txt"
    bed = []
    last_entry = ""
    read_lens = {}
    seqs = {}

    with open(sam_loc, "r") as f:
#        lines = f.readlines()
        for line in f:
            if line[0] == "@":
                continue
            buff = line.strip().split()
            
            name = buff[0].split("/")[0]
            chrom = buff[2]
            start = buff[3]
            end = str(int(buff[3]) + len(buff[9]))
            seq = buff[9]
            check_as = next((s for s in buff if "AS" in s), None)
            if check_as is not None:
                align_score = check_as.split(":")[-1]
            else:
                continue
            
            if int(end) < int(start) or "chr" not in chrom:
                continue
            curr_entry = "\t".join([chrom, start, end, name, align_score])
            if last_entry == curr_entry:
                continue
            last_entry = curr_entry
            bed += [last_entry]
            read_lens[name] = int(end) - int(start)
            seqs[name] = seq


    with open(bed_loc, "w") as f:
        for line in bed:
            f.write(line + "\n")
    with open(len_loc, "w") as f:
        for name, l in read_lens.items():
            f.write(name + "," + str(l) + "\n")
    with open(seq_loc, "w") as f:
        for name, seq in seqs.items():
            f.write(name + "," + seq + "\n")


def split_overlap(base_loc, overlap):
    bed_loc = base_loc + "/aligned.bed"
    
    unique_count_loc = base_loc + "unique.counts"
    multi_loc = base_loc + "multi.translation"
    unique_translate_loc = base_loc + "unique.translation"
    unique, unique_translation, multi = {}, {}, {}
    
    for k,v in overlap.items():
        if len(v) == 1:
            val = list(v.keys())[0]
            if val not in unique:
                unique[val] = 1
            else:
                unique[val] += 1
            unique_translation[k] = val
        else:
            multi[k] = [str(val) + "*" + key for key, val in v.items()]

    with open(unique_count_loc, "w") as f:
        for k,v in unique.items():
            f.write(k + "," + str(v) + "\n")
    with open(multi_loc, "w") as f:
        for k,v in multi.items():
            f.write(k +"," + ",".join(v) + "\n")


def find_overlap(base_loc, gtf_loc):
    bed_loc = base_loc + "aligned.bed"
    intersect_loc = base_loc + "intersect.out"
    log_loc = base_loc + "intersect.log"
    if "ervmap" or "telescope" in str(gtf_loc):
        print("wtaf")
        chr_loc = "/home/stexocae/li_lab/saem/refs/star_hs1_idx/chrNameLength.txt"
    else:
        chr_loc = "/".join(str(gtf_loc).split("/")[:-1]) + "/star_hs1_idx/chrNameLength.txt"

    overlap = {}
    bedtools_call  = "bedtools intersect -sorted -a " + bed_loc + " -b " + str(gtf_loc) + " "
    bedtools_call += "-wo -f 0.5 > " + intersect_loc + " -g " + chr_loc

#    bedtools_call += "2> " + log_loc

    print("Running external bedtools command with the following parameters:\n")
    print(bedtools_call + "\n")
    subprocess.call(bedtools_call, shell=True)

    with open(intersect_loc, "r") as f:
        for line in f.readlines():
            buff = line.strip().split()
            if buff[3] not in overlap:
                overlap[buff[3]] = {buff[16][1:-2] : buff[4]}
            else:
                overlap[buff[3]] |= {buff[16][1:-2] : buff[4]}
    
    split_overlap(base_loc, overlap) 


def clean(base_loc):
    rm_files = [base_loc + "aligned.bed", base_loc + "intersect.out", base_loc + "star.sam"]
    for file in rm_files:
        subprocess.call("rm " + file, shell=True)


## driver fxn ##
if __name__ == '__main__':

    ## Argparse ##
    parser = argparse.ArgumentParser(prog="TE-SAEM align parse",
             description = "Parser of aligned data for TE-SAEM",
             epilog = "For help//queries please visit our github: https://github.com/synanth/TE-SAEM")
    parser.add_argument("-s",
                        type = pathlib.Path,
                        required = True,
                        metavar = "sam",
                        dest = "sam",
                        help = "Location of SAM file (default: %(default)s)")
    parser.add_argument("-t",
                        type = int,
                        default = 1,
                        metavar = "threads",
                        dest = "threads",
                        help = "Number of threads to use (default: %(default)s)")
    parser.add_argument("-g",
                        type = pathlib.Path,
                        required = True,
                        metavar = "gtf",
                        dest = "gtf",
                        help = "Location of TE GTF file (default: %(default)s)")
    
    args = parser.parse_args()
    base_loc = "/".join(str(args.sam).split("/")[:-1]) + "/"
    
    ## Parse SAM file ##
    print("Sorting SAM file")
    sort_sam(base_loc, args)
    print("Converting SAM to BED")    
    sam_to_bed(base_loc)
    print("Finding overlap between reads and GTF annotation")
    find_overlap(base_loc, args.gtf)
    print("parse_done")
#    clean(base_loc)
