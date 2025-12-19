import sys
import math
import random


## em ##
def init_abundance(multimapped_reads, all_tes):
    theta = {k:0 for k in all_tes}
    for te in all_tes:
        theta[te] = random.random()
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def log_likelihood(theta, multimapped_reads, align_scores, e_lens, beta=1.0):
    ll = 0.0
    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            xs.append(math.log(theta[te]) + beta * align_scores[read][te] -  math.log(e_lens[read][te]))
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    return ll


def e_step(theta, multimapped_reads, align_scores, e_lens, beta = 1.0):
    frac = {k:{} for k in multimapped_reads}
    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            xs.append(math.log(theta[te]) + beta * align_scores[read][te] - math.log(e_lens[read][te]))

        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes, xs):
            frac[read][te] = math.exp(x - m) / Z

    return frac


def m_step(frac, multimapped_reads, all_tes, unique_counts):
    theta = {k: 0 for k in all_tes}

    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]

    theta = {k:max(v, 1e-12) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def em(e_lens, multimapped_reads, unique_counts, align_scores):
    threshold = 1e-2

    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))

    old_abundance = init_abundance(multimapped_reads, all_tes)
    ll_old = log_likelihood(old_abundance, multimapped_reads, align_scores, e_lens)

    for i in range(1,10000):
        frac = e_step(old_abundance, multimapped_reads, align_scores, e_lens)
        new_abundance = m_step(frac, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, multimapped_reads, align_scores, e_lens)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)

        if abs(diff) < threshold:
            return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}

        old_abundance = new_abundance
        ll_old = ll_new
    return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}




def calc_e_lens(len_transcripts, read_lens, multimapped_reads):
    e_lens = {}
    for read,tes in multimapped_reads.items():
        e_lens[read] = {}
        for te in tes:
            e_lens[read][te] = max(len_transcripts[te] - read_lens[read] + 1,1)
    return e_lens


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
    multi_seq = {}
    multimapped_reads = {}
    align_scores = {}
    read_lens = {}
    exclude = []


    with open(sam_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "@":
                continue
            buff = line.strip().split()
            name = buff[0].split("/")[0]
            read_lens[name] = len(buff[9])

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
            len_transcripts[name] = int(buff[4]) - int(buff[3]) +1
    
    with open(unique_counts_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split(",")
            unique_counts[buff[0]] = int(buff[1])
    
    with open(multimapped_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split(",")
            read = buff[0]
            te_names = [x.split("/")[1] for x in buff[1:]]
            te_scores = [x.split("/")[0] for x in buff[1:]]
            multimapped_reads[read] = te_names
            align_scores[read] = {te:float(te_scores[e]) for e, te in enumerate(te_names)}

    for key in exclude:
        multimapped_reads.pop(key, None)


    e_lens = calc_e_lens(len_transcripts, read_lens, multimapped_reads)

    em_counts = em(e_lens, multimapped_reads, unique_counts, align_scores)
    total_em = sum(em_counts.values())
    em_frac = {k:v/total_em for k,v in em_counts.items()}
    all_tes = {k:0 for k in len_transcripts}

    for te in all_tes:
        uc = unique_counts.get(te,0)
        frac = em_frac.get(te,0)
        if uc > 0 or frac > 5e-5:
            all_tes[te] += uc + int(em_counts.get(te,0))


    out_loc = base_loc + "manuscript/out/saem_em.out"
    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
