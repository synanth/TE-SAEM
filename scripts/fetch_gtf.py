import sys
import subprocess
import os
import pysam

def calc_gc(seq):
    return (seq.count("G") + seq.count("g") + seq.count("C") + seq.count("c"))/len(seq)


def convert_rmsk_to_gtf(rmsk_loc, gtf_loc, fa_loc):
    fa = pysam.FastaFile(fa_loc)
    skip = ["DNA", "DNA?", "Low_complexity", "LTR", "LTR?", "rRNA", "Satellite", "Simple_repeat", "scRNA", "snRNA", "srpRNA", "tRNA", "Unknown", "Unspecified"]
    gtf = []
    te_dup_number = {}
    with open(rmsk_loc, "r") as f:
        lines = f.readlines()
        for line in lines[3:]:
            buff = line.strip().split()
            if buff[10] in skip:
                continue
            chrm = buff[4]
            start = buff[5]
            end = buff[6]
            strand = buff[8]
            score = buff[0]
            if strand == "C":
                strand = "-"
            name = buff[9]
            class_fam = buff[10].split("/")
            te_class = class_fam[0]
            fam = class_fam[1]
            if name not in te_dup_number.keys():
                te_dup_number[name] = 1
                te_id = name
            else:
                te_id = name + "_dup" + str(te_dup_number[name])
                te_dup_number[name] += 1
            gc = calc_gc(fa.fetch(chrm, int(start), int(end)))
            last_col = 'gene_id "' + name + '"; transcript_id "' + te_id + '"; family_id "' + fam + '"; class_id "' +  te_class + '"; gene_name "' + name + ':TE"; locus "' + te_id + '"; gc_content "' + f"{gc:.3f}" + '";'
            gtf += [[chrm, "rmsk", "exon", start, end, score, strand, ".", last_col]]
    gtf = sorted(gtf, key = lambda x: (int(x[0][3:]) if x[0][3:].isnumeric() else ord(x[0][3:]), int(x[3])))

    with open(gtf_loc, "w") as f:
        for line in gtf:
            f.write("\t".join(line) + "\n")


def extract_lens(gtf_loc):
    len_loc = gtf_loc.replace("gtf", "len")
    lens = []
    with open(gtf_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split()
            lens += [int(buff[4]) - int(buff[3]) + 1]

    with open(len_loc, "w") as f:
        for l in lens:
            f.write(str(l) + "\n")


def check_annotation(loc):
    generate_annotation()
    quit()


def generate_annotation(rmsk_loc, gtf_loc, fa_loc):
    print("Downloading CHM13 repeatmasker output from UCSC")
    wget_call = "wget -q --show-progress -P " + refs_loc + " https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz"
    #subprocess.run(wget_call, shell=True)
    print("\nUnzipping repeatmasker output")
    gunzip_call = "gunzip " + rmsk_loc
    #subprocess.run(gunzip_call, shell=True)
    rmsk_loc = rmsk_loc[:-3]
    print("\nConverting repeatmasker output to gtf")
    convert_rmsk_to_gtf(rmsk_loc, gtf_loc, fa_loc)
    print("\nExtract TE lengths")
    extract_lens(gtf_loc)




if __name__ == '__main__':

    print(len(sys.argv))
    print(sys.argv)
    quit()












    chm13_rmsk_loc = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz"
    refs_loc = "/home/stexocae/li_lab/saem/refs/"
    rmsk_loc = refs_loc + "hs1.repeatMasker.out.gz"
    gtf_loc = refs_loc + "hs1.gtf"
    fa_loc = refs_loc + "hs1.fa"

    clear_call = "clear"
    subprocess.run(clear_call, shell=True)
    
    print("*****************************")
    print("***  Running SAEM setup   ***")
    print("*****************************\n")
    #generate_genome_idx(refs_loc)

    generate_annotation(rmsk_loc, gtf_loc, fa_loc)
    #convert_rmsk_to_gtf(rmsk_loc[:-3], gtf_loc, fa_loc)
    quit()
    print("Checking for .ini in home folder")
    params = check_ini()
    
    print("Checking for STAR index")
    check_genome_idx(params["genome_idx"])

    print("Checking for TE annotation file")
    check_annotation(gtf_loc)
