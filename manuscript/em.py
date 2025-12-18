import sys
import math
import random


## em ##
def gc_weight(gc, gc_bias, bin_width=.01):
    bin_idx = round(gc/bin_width) * bin_width
    if bin_idx not in gc_bias:
        bin_idx = min(gc_bias.keys(), key=lambda x: abs(x-bin_idx))
    return gc_bias[bin_idx]


def log_likelihood(theta, len_transcripts, multimapped_reads, read_lens, gc, bias):
    ll = 0.0
    for read, tes in multimapped_reads.items():
        xs = [
            math.log(theta[te] * gc_weight(gc[te], gc_bias)) -
            math.log(max(len_transcripts[te] - read_lens[read] + 1, 1))
            for te in tes
        ]
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    ll = ll/len(multimapped_reads)
    return ll


def init_abundance(len_transcripts, multimapped_reads, all_tes):
    theta = {k:0 for k in all_tes}
    avg_len = 120
    for te in all_tes:
        theta[te] = random.random()
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, len_transcripts, read_lens, gc, gc_bias):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            e_len = max(len_transcripts[te] - read_lens[read] + 1, 1)
            xs.append(math.log(theta[te] * gc_weight(gc[te], gc_bias)) - math.log(e_len))

        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes, xs):
            frac[read][te] = math.exp(x - m) / Z

    return frac


def m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes, unique_counts):
    theta = {k: 0 for k in all_tes}

    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]
    theta = {k:max(v, 1e-12) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def em(len_transcripts, read_lens, multimapped_reads, unique_counts, gc,gc_bias):
    threshold = 0.0001
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    old_abundance = init_abundance(len_transcripts, multimapped_reads, all_tes)
    ll_old = log_likelihood(old_abundance, len_transcripts, multimapped_reads, read_lens, gc, gc_bias)

    for i in range(1,10000):
        frac = e_step(old_abundance, multimapped_reads, len_transcripts, read_lens, gc, gc_bias)
        new_abundance = m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, len_transcripts, multimapped_reads, read_lens, gc, gc_bias)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)

        if abs(diff) < threshold:
            return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}
        old_abundance = new_abundance
        ll_old = ll_new
    return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}

def calc_gc_frac(seq):
    return (seq.count("g") + seq.count("G") + seq.count("c") + seq.count("C"))/len(seq)


def calc_gc_bias(unique_counts, unique_seq, bin_width=0.01):
    gc_bin_counts = {}

    for read, seq in unique_seq.items():
        gc_frac = calc_gc_frac(seq)
        bin_idx = round(gc_frac/bin_width) * bin_width
        if bin_idx not in gc_bin_counts:
            gc_bin_counts[bin_idx] = 1
        else:
            gc_bin_counts[bin_idx] += 1

    mean_count = sum(gc_bin_counts.values()) / len(gc_bin_counts)
    gc_bias = {k:v/mean_count for k,v in gc_bin_counts.items()}
    return gc_bias


## driver fxn ##
if __name__ == '__main__':
    base_loc = "/home/stexocae/li_lab/saem/"
    gtf_loc = base_loc + "refs/hs1.gtf"
    unique_counts_loc = base_loc + "sim_data/star.unique.counts"
    multimapped_loc = base_loc + "sim_data/star.multi.translation"
    sam_loc = base_loc + "sim_data/alignment/star.sam"
    assembly_loc = base_loc + "sim_data/assembly/to_contigs.sam"

    len_transcripts = {}
    unique_counts = {}
    unique_seq = {}
    multimapped_reads = {}
    read_lens = {}
    e_lens = {}
    gc = {}

    with open(sam_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "@":
                continue
            buff = line.strip().split()
            name = buff[0].split("/")[0]
            read_lens[name] = int(buff[8])
            if buff[-4] == "NH:i:1":
                unique_seq[name] = buff[9]
    
    with open(assembly_loc, "r") as f:
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
            name = buff[11][1:-2] 
            len_transcripts[name] = abs(int(buff[4]) - int(buff[3]))
            gc[name] = float(buff[-1][1:-2])

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


    gc_bias = calc_gc_bias(unique_counts, unique_seq)

    em_counts = em(len_transcripts, read_lens, multimapped_reads, unique_counts,gc, gc_bias)
    non_zero = {k:int(v) for k,v in em_counts.items() if v > 3}

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
