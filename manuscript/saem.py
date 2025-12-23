import sys
import math
import random

## sa ##
def sa(old_abundance, temp):
    sa_abundance = old_abundance.copy()
    n_neighbors = max(1, min(round(len(sa_abundance)*temp),1000))
    tes = [k for k in sa_abundance.keys() if k != "_noise"]
    vals = sa_abundance.values()
    
    mean_val = sum(vals)/len(tes)
    step_size = math.sqrt(sum((v- mean_val)**2 for v in vals))


    for idx in random.sample(tes, n_neighbors):
        change = random.uniform(-step_size*temp, step_size*temp)
        sa_abundance[idx] = max(sa_abundance[idx] +change, 1e-100)
   
    sum_norm = sum(sa_abundance.values())
    sa_abundance["_noise"] = old_abundance["_noise"]
    sa_abundance= {k:v/sum_norm for k,v in sa_abundance.items()}
    return sa_abundance


def reduce_temp(temp, cooling_rate, i, ll_diff, c=2.0):
    ## testing geometric cooling ##
    new_temp = cooling_rate * temp
    ## testing log-geometric cooling ##
#    new_temp = (temp/math.log(i + c)) * (cooling_rate **i)
    ## cost aware cooling ##
#    new_temp = (temp/math.log(i+c+abs(ll_diff)))
    return new_temp 


def accept_sa(ll_old, ll_sa, temp, n):
    delta_ll = max(ll_sa -ll_old, -700)
    if delta_ll > 0:
        print(">")
        return True
    if abs(delta_ll) < 1e-6:
        return False
    print(-delta_ll/temp)
    accept_rate = math.exp(-delta_ll/temp)
    if random.random() < accept_rate:
        print(delta_ll, accept_rate)
        return True
    return False


## em ##
def gc_weight(gc, gc_bias, gc_w=0.1, min_w=0.5,max_w=2.0):
    g_bin = round(gc/gc_w) * gc_w

    if g_bin in gc_bias:
        return gc_bias[g_bin]
    bin_idx = min(gc_bias.keys(), key=lambda x: abs(x-g_bin))
    return max(min(gc_bias[bin_idx], max_w), min_w)


def log_likelihood(theta, multimapped_reads, gc_weights, align_scores, e_lens):
    ll = 0.0

    for read, tes in multimapped_reads.items():
        xs = [
            math.log(theta[te])
            + gc_weights[te] 
            + align_scores[read][te]
            - e_lens[read][te]
            for te in tes
        ]
        xs.append(math.log(theta["_noise"])
        + align_scores[read]["_noise"]
        - e_lens[read]["_noise"]
        )
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    ll = ll / len(multimapped_reads)
    return ll


def init_abundance(multimapped_reads, all_tes):
    theta = {k:random.random() for k in all_tes}
    theta["_noise"] = 1e-3
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, gc_weights, align_scores, e_lens):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            xs.append(math.log(theta[te]) 
            + gc_weights[te]
            + align_scores[read][te]
            - e_lens[read][te])

        xs.append(math.log(theta["_noise"])
        + align_scores[read]["_noise"]
        - e_lens[read]["_noise"]
        )
        tes_plus = tes + ["_noise"]
        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes_plus, xs):
            frac[read][te] = math.exp(x - m) / Z
    return frac


def m_step(frac, multimapped_reads, all_tes, unique_counts, alpha = 1):
    theta = {k: alpha*math.sqrt(unique_counts.get(k, 0)) for k in all_tes}
    theta["_noise"] = 1e-2
    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]
        theta["_noise"] += frac[read]["_noise"]
    theta = {k:max(v, 1e-9) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def theta_to_counts(frac, all_tes):
    counts = {k: 0 for k in all_tes}
    print(len(frac))
    rescue = []
    for read, tes in frac.items():
        vals = sorted(tes.items(), key = lambda x: x[1], reverse = True)
        if vals[0][1] - vals[1][1] > .25:
            counts[vals[0][0]] += 1
        else:
#            print(read, tes)
            rescue += [read]
    print(sum(counts.values()))
    return set(rescue), counts


def em(multimapped_reads, unique_counts, cooling_rate, gc_weights, align_scores, e_lens):
    threshold = 1e-6
    temp = 1.0
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    all_tes.append("_noise")
    old_abundance = init_abundance( multimapped_reads, all_tes)

    for i in range(1,10000):
        sa_abundance = sa(old_abundance, temp)
        ll_sa = log_likelihood(sa_abundance, multimapped_reads, gc_weights, align_scores, e_lens)
        ll_old = log_likelihood(old_abundance, multimapped_reads, gc_weights, align_scores, e_lens)
        if accept_sa(ll_old, ll_sa, temp, len(multimapped_reads)):
            print(i, ll_old, ll_sa, temp)
            old_abundance = sa_abundance
            temp = reduce_temp(temp, cooling_rate, i, (ll_sa - ll_old))
            continue
        frac = e_step(old_abundance, multimapped_reads, gc_weights, align_scores, e_lens)
        new_abundance = m_step(frac, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, multimapped_reads, gc_weights, align_scores, e_lens)
        diff = ll_new - ll_old
        print(i, diff, ll_old, ll_new)
        if abs(diff) < threshold and max(abs(new_abundance[k] - old_abundance[k]) for k in new_abundance) < 1e-4:
            return theta_to_counts(e_step(new_abundance, multimapped_reads, gc_weights, align_scores, e_lens), all_tes)
        temp = reduce_temp(temp, cooling_rate, i, diff)
        old_abundance = new_abundance
    return theta_to_counts(e_step(new_abundance, multimapped_reads, gc_weights, align_scores, e_lens), all_tes)


def gc_content(read):
    return (read.count("G") + read.count("C"))/max(1, len(read))


def calc_gc_bias(unique_seq, gc_bin_size=0.01, pseudo_count=5):
    counts = {}

    for seq in unique_seq.values():
        gc_frac = gc_content(seq)
        g = round(gc_frac/gc_bin_size) * gc_bin_size
        counts[g] = counts.get(g,0) + 1
    for g in counts:
        counts[g] += pseudo_count
    counts_mean = sum(counts.values())/len(counts)
    return {k:math.log(v/counts_mean) for k,v in counts.items()}


def norm_align_scores(align_scores):
    align_weight = {}
    for read, tes in align_scores.items():
        max_score = max(tes.values())
        weights = {te: int(asv)/int(max_score) for te, asv in tes.items()}
        weights["_noise"] = min(weights.values())/10
        align_weight[read] = weights
    return align_weight


def calc_e_lens(multimapped_reads, len_transcripts, read_lens):
    e_lens = {k:{} for k in multimapped_reads}
    len_max = max(len_transcripts.values())
    for read, tes in multimapped_reads.items():
        for te in tes:
            e_lens[read][te] = math.log(max(len_transcripts[te] - read_lens[read] + 1, 1))
    for read in multimapped_reads:
        e_lens[read]["_noise"] = math.log(max(len_max - read_lens[read] +1,1))
    return e_lens


## driver fxn ##
if __name__ == '__main__':

    cooling_rate = float(sys.argv[1])
    n_run = sys.argv[2]

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
    gc = {}

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
            name = buff[0]
            tes = [x.split("*")[1] for x in buff[1:]]
            scores = [x.split("*")[0] for x in buff[1:]]
            multimapped_reads[name] = list(set(tes))
            align_scores[name] = {te:scores[e] for e, te in enumerate(tes)}
    
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
    align_scores = norm_align_scores(align_scores)
    gc_bias = calc_gc_bias(unique_seq)
    gc_weights = {te:.2*gc_weight(gc[te], gc_bias) for te in len_transcripts}
    e_lens = calc_e_lens(multimapped_reads, len_transcripts, read_lens)
    em_counts = {}
    rescue, em_counts = em(multimapped_reads, unique_counts, cooling_rate, gc_weights, align_scores, e_lens)
    print(em_counts["_noise"])
    if "_noise" in em_counts:
        del em_counts["_noise"]
#    rescue_reads = {k:v for k,v in multimapped_reads.items() if k in rescue}
#    print(len(rescue))
#    toss, rescue_counts = em(rescue_reads, unique_counts, cooling_rate, gc_weights, align_scores, e_lens)

    all_tes = {k:unique_counts.get(k, 0) + em_counts.get(k, 0) for k in len_transcripts}

    out_loc = base_loc + "manuscript/out/saem_" + n_run + ".out"

    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
