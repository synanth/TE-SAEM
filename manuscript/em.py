import sys
import math
import random


## em ##

def log_likelihood(theta, multimapped_reads, e_lens):
    log_sum = 0
    for read,tes in multimapped_reads.items():
        read_sum = 0
        for te in tes:
            read_sum += theta[te]/e_lens[te]
        log_sum += math.log(read_sum + 1e-300)
    return log_sum


def init_abundance(multimapped_reads, all_tes):
    theta = {k:0 for k in all_tes}

    for te in all_tes:
        theta[te] = random.random()
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, e_lens):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        frac_sum = 0
        for te in tes:
            frac[read][te] = theta[te]/e_lens[te]
            frac_sum += frac[read][te]
        for te in tes:
            frac[read][te] /= frac_sum
    return frac


def m_step(frac, multimapped_reads, all_tes):
    theta = {k: 0 for k in all_tes}
    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def em(multimapped_reads, e_lens):
    threshold = 0.01
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    
    old_abundance = init_abundance(multimapped_reads, all_tes)
    ll_old = log_likelihood(old_abundance, multimapped_reads, e_lens)

    for i in range(1,10000):
        frac = e_step(old_abundance, multimapped_reads, e_lens)
        new_abundance = m_step(frac, multimapped_reads, all_tes)
        ll_new = log_likelihood(new_abundance, multimapped_reads, e_lens)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)

        if abs(diff) < threshold:
            return {k:v*len(multimapped_reads) for k,v in old_abundance.items()}
        old_abundance = new_abundance
        ll_old = ll_new
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

    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    e_lens = {te:0 for te in all_tes}
    avg_len = int(sum(read_lens.values())/len(read_lens))
    for te in e_lens:
        e_lens[te] = max(len_transcripts[te] - avg_len+1, 1)

    em_counts = em(multimapped_reads, e_lens)
    non_zero = {k:int(v) for k,v in em_counts.items()}

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
