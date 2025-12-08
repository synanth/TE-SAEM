import sys
import math
import random


## em ##

def log_likelihood(theta, len_transcripts, multimapped_reads, readlen, avg_len, unique_counts):
    log_sum = 0
    for read, tes in multimapped_reads.items():
        log_sum += math.log(sum([theta[te] / max(len_transcripts[te]-read_lens[read]+1,1) for te in tes]))
    return log_sum


def init_abundance(len_transcripts, multimapped_reads, avg_len):
    tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    theta = {k:0 for k in tes}

    for te in tes:
        theta[te] = random.random()
#        theta[te] = (1/max(len_transcripts[te]-avg_len+1,1)+1e-12) 
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, len_transcripts, read_lens, avg_len):
    frac = {k:0 for k in theta}

    for read, tes in multimapped_reads.items():
        denom = sum([theta[te]/max(len_transcripts[te]-read_lens[read]+1,1) for te in tes])
        for te in tes:
            e_len = max(len_transcripts[te]-read_lens[read]+1,1)
            frac[te] += (theta[te] / e_len) / denom
    return frac


def m_step(frac, len_transcripts, read_lens, multimapped_reads, avg_len, unique_counts):
    theta = {k: 0 for k,v in frac.items()}
    denom = sum([frac[te] /max(len_transcripts[te]-avg_len+1,1) for te in theta]) + 1e-12

    for te in theta:
        theta[te] = (frac[te]/max(len_transcripts[te]-avg_len+1,1)+1e-12) / denom
    return theta


def em(len_transcripts, read_lens, multimapped_reads, unique_counts):
    threshold = 0.01
    avg_len = sum(read_lens.values())/len(read_lens)
    old_abundance = init_abundance(len_transcripts, multimapped_reads, avg_len)

    for i in range(1,10000):
#    for i in range(1,20):
        
        frac = e_step(old_abundance, multimapped_reads, len_transcripts, read_lens, avg_len)
        new_abundance = m_step(frac, len_transcripts, read_lens, multimapped_reads, avg_len, unique_counts)
        ll_new = log_likelihood(new_abundance, len_transcripts, multimapped_reads, read_lens, avg_len, unique_counts)
        ll_old = log_likelihood(old_abundance, len_transcripts, multimapped_reads, read_lens, avg_len, unique_counts)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)
        old_abundance = new_abundance

        if abs(diff) < threshold:
            return {k:v*len(multimapped_reads) for k,v in old_abundance.items()}
    return e_step(old_abundance, multimapped_reads, len_transcripts, read_lens, avg_len)
    return {k:v*len(multimapped_reads) for k,v in old_abundance.items()}



## driver fxn ##
if __name__ == '__main__':
    base_loc = "/home/stexocae/li_lab/saem/"
    gtf_loc = base_loc + "refs/hs1.gtf"
    unique_counts_loc = base_loc + "sim_data/star.unique.counts"
    multimapped_loc = base_loc + "sim_data/star.multi.translation"
    sam_loc = base_loc + "sim_data/alignment/segemehl.sam"

    len_transcripts = {}
    unique_counts = {}
    multimapped_reads = {}
    read_lens = {}
    e_lens = {}

    with open(sam_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "@":
                continue
            buff = line.strip().split()
            name = buff[0].split("/")[0]
            read_lens[name] = len(buff[9])
    
    with open(gtf_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split()
            len_transcripts[buff[11][1:-2]] = abs(int(buff[4]) - int(buff[3]))

    with open(unique_counts_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split(",")
            unique_counts[buff[0]] = int(buff[1])

    with open(multimapped_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split(",")
            multimapped_reads[buff[0]] = buff[1:]


    em_counts = em(len_transcripts, read_lens, multimapped_reads, unique_counts)
    non_zero = {k:int(v) for k,v in em_counts.items()}
    print(len(multimapped_reads))

    all_tes = {k:0 for k in len_transcripts}
    for te in all_tes:
        if te in unique_counts:
            all_tes[te] = unique_counts[te]
        if te in non_zero:
            all_tes[te] += non_zero[te]

    out_loc = base_loc + "manuscript/out/saem_em.out"
    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
