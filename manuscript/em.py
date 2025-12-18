import sys
import math
import random


## em ##
def gc_weight(gc, gc_bias, l, gc_w=.01, min_w=0.5,max_w=2.0, l_w=25):
    g_bin = round(gc/gc_w) * gc_w
    l_bin = (l // l_w) / l_w
    if (g_bin, l_bin) in gc_bias:
        return(gc_bias[(g_bin, l_bin)])
    bin_idx = min(gc_bias.keys(), key=lambda x: abs(x[0]-g_bin)+abs(x[0]-l_bin))
    return max(min(gc_bias[bin_idx], max_w), min_w)


def log_likelihood(theta, len_transcripts, multimapped_reads, read_lens, gc_weights, align_scores):
    ll = 0.0
    for read, tes in multimapped_reads.items():
        xs = [
            math.log(theta[te]) + math.log(align_scores[read][te]) + gc_weights[read][te] -
            math.log(max(len_transcripts[te] - read_lens[read] + 1, 20))
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


def e_step(theta, multimapped_reads, len_transcripts, read_lens, gc_weights, align_scores):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            e_len = max(len_transcripts[te] - read_lens[read] + 1, 20)
            xs.append(math.log(theta[te]) + math.log(align_scores[read][te]) + gc_weights[read][te] - math.log(e_len))

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


def em(len_transcripts, read_lens, multimapped_reads, unique_counts, gc_weights, align_scores):
    threshold = 0.000001
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    old_abundance = init_abundance(len_transcripts, multimapped_reads, all_tes)
    ll_old = log_likelihood(old_abundance, len_transcripts, multimapped_reads, read_lens, gc_weights, align_scores)

    for i in range(1,10000):
        frac = e_step(old_abundance, multimapped_reads, len_transcripts, read_lens, gc_weights, align_scores)
        new_abundance = m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, len_transcripts, multimapped_reads, read_lens, gc_weights, align_scores)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)

        if abs(diff) < threshold:

            return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}
        old_abundance = new_abundance
        ll_old = ll_new
    return {k:v*len(multimapped_reads) for k,v in new_abundance.items()}

def calc_gc_frac(seq):
    return (seq.count("g") + seq.count("G") + seq.count("c") + seq.count("C"))/len(seq)


def calc_gc_bias(unique_counts, unique_seq, gc_w=0.01, len_w=25, pseudo=5):
    counts = {}
    total = 0

    for read, seq in unique_seq.items():
        gc_frac = calc_gc_frac(seq)
        g = round(gc_frac/gc_w)*gc_w
        l = (len(seq)//len_w) * len_w
        counts[(g,l)] = counts.get((g,l),0)+1
        total += 1
    g_set = set(k[0] for k in counts.keys())        
    l_set = set(k[0] for k in counts.keys())        
    for g in g_set:
        for l in l_set:
            counts[(g,l)] = counts.get((g,l), 0) + pseudo

    gc_bias = {k:v/total for k,v in counts.items()}
    return gc_bias


def norm_align_scores(align_scores, tau=5.0):
    align_weight = {}

    for read, tes in align_scores.items():
        max_score = max(tes.values())
        weights = {te: math.exp((int(asv)-int(max_score))/tau) for te, asv in tes.items()}
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
            align_scores[read] = {te:te_scores[e] for e, te in enumerate(te_names)}

    align_scores = norm_align_scores(align_scores)
    gc_bias = calc_gc_bias(unique_counts, unique_seq)
    gc_weights = {read:{te:math.log(gc_weight(gc[te], gc_bias, read_lens[read])) for te in tes} for read,tes in multimapped_reads.items()}

    em_counts = em(len_transcripts, read_lens, multimapped_reads, unique_counts, gc_weights, align_scores)
    total_em = sum(em_counts.values())
    em_frac = {k:v/total_em for k,v in em_counts.items()}
    all_tes = {k:0 for k in len_transcripts}

    for te in all_tes:
        uc = unique_counts.get(te,0)
        frac = em_frac.get(te,0)
        if uc > 0 or frac > 5e-3:
            all_tes[te] += uc + int(em_counts.get(te,0))

    


    out_loc = base_loc + "manuscript/out/saem_em.out"
    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
