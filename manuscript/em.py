import sys
import math
import random


## em ##
def gc_weight(gc, gc_bias, gc_w=.01, min_w=0.5,max_w=2.0):
    g_bin = round(gc/gc_w) * gc_w
    if g_bin in gc_bias:
        return gc_bias[g_bin]
    bin_idx = min(gc_bias.keys(), key=lambda x: abs(x-g_bin))
    return max(min(gc_bias[bin_idx], max_w), min_w)


def log_likelihood(theta, len_transcripts, multimapped_reads, read_lens, gc_weights, align_scores, gc):
    ll = 0.0
    for read, tes in multimapped_reads.items():
        max_score = max(align_scores[read].values())
        tau = 1
        xs = [
            math.log(theta[te]) + (align_scores[read][te] - max_score)/tau -
            math.log(max(len_transcripts[te] - read_lens[read] + 1, 1))
            for te in tes
        ]
        xs.append(math.log(theta["_noise"]) + -5)
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    ll = ll/len(multimapped_reads)
    return ll


def init_abundance(len_transcripts, multimapped_reads, all_tes):
    theta = {k:0 for k in all_tes}
    random.seed(1)

    for te in all_tes:
        theta[te] = random.random()
    theta["_noise"] = 0.05
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, len_transcripts, read_lens, gc_weights, align_scores, gc):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        xs = []
        max_score = max(align_scores[read].values())
        tau = 1
        for te in tes:
            e_len = max(len_transcripts[te] - read_lens[read] + 1, 1)
            xs.append(math.log(theta[te]) + (align_scores[read][te]-max_score)/tau  - math.log(e_len))

        xs.append(math.log(theta["_noise"]) + -5)
        tes_plus = tes + ["_noise"]
        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes_plus, xs):
            frac[read][te] = math.exp(x - m) / Z

    return frac


def m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes, unique_counts):
    alpha = .3
    eps = 1e-12
    theta = {k: 0 for k in all_tes}
#    theta = {k: unique_counts.get(k,0) + (alpha-1) for k in all_tes}

    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]
        theta["_noise"] += frac[read]["_noise"]

    theta = {k:max(v, 1e-12) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def em(len_transcripts, read_lens, multimapped_reads, unique_counts, gc_weights, align_scores,gc):
    threshold = 1e-6
    noise = "_noise"
    log_p_noise = -5.0

    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    all_tes.append(noise)

    old_abundance = init_abundance(len_transcripts, multimapped_reads, all_tes)
    ll_old = log_likelihood(old_abundance, len_transcripts, multimapped_reads, read_lens, gc_weights, align_scores, gc)

    for i in range(1,10000):
        frac = e_step(old_abundance, multimapped_reads, len_transcripts, read_lens, gc_weights, align_scores, gc)
        new_abundance = m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, len_transcripts, multimapped_reads, read_lens, gc_weights, align_scores, gc)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)

        if abs(diff) < threshold:

            return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}
        old_abundance = new_abundance
        ll_old = ll_new
    return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}

def calc_gc_frac(seq):
    return (seq.count("g") + seq.count("G") + seq.count("c") + seq.count("C"))/len(seq)


def calc_gc_bias(unique_seq, gc_bin_size=0.01, pseudo_count=5):
    counts = {}

    for seq in unique_seq.values():
        gc_frac = calc_gc_frac(seq)
        g = round(gc_frac/gc_bin_size) * gc_bin_size
        counts[g] = counts.get(g,0) + 1
    for g in counts:
        counts[g] += pseudo_count
    counts_mean = sum(counts.values())/len(counts)

    return {k:v/counts_mean for k,v in counts.items()}


def norm_align_scores(align_scores, tau=.5):
    align_weight = {}

    for read, tes in align_scores.items():
        max_score = max(tes.values())
        weights = {te: math.exp((int(asv)-int(max_score))/tau) for te, asv in tes.items()}
        weight_sum = sum(weights.values())
        weights = {k:v/weight_sum for k,v in weights.items()}
        align_weight[read] = weights
        
    return align_weight


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
    e_lens = {}
    gc = {}

    with open(sam_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "@":
                continue
            buff = line.strip().split()
            name = buff[0].split("/")[0]
            read_lens[name] = len(buff[9])
            if buff[-4] == "NH:i:1":
                unique_seq[name] = buff[9]
                multi_seq[name] = buff[9]
            else:
                multi_seq[name] = buff[9]
    
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
            len_transcripts[name] = len(buff[9])
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
            read = buff[0]
            te_names = [x.split("/")[1] for x in buff[1:]]
            te_scores = [x.split("/")[0] for x in buff[1:]]
            multimapped_reads[read] = te_names
            align_scores[read] = {te:float(te_scores[e]) for e, te in enumerate(te_names)}


    gc_bias = calc_gc_bias(unique_seq)
    gc_weights = {read:math.log(gc_weight(calc_gc_frac(multi_seq[read]), gc_bias)) for read in multi_seq}
    gc = {k:math.log(gc_weight(v,gc_bias)) for k,v in gc.items()}

    em_counts = em(len_transcripts, read_lens, multimapped_reads, unique_counts, gc_weights, align_scores, gc)
    if "_noise" in em_counts:
        del em_counts["_noise"]
    total_em = sum(em_counts.values())
    em_frac = {k:v/total_em for k,v in em_counts.items()}
    all_tes = {k:0 for k in len_transcripts}

    for te in all_tes:
        uc = unique_counts.get(te,0)
        frac = em_frac.get(te,0)
        if uc > 0 or frac > 5e-4:
            all_tes[te] += uc + int(em_counts.get(te,0))

    


    out_loc = base_loc + "manuscript/out/saem_em.out"
    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
