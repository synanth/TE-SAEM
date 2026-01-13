import sys
import math
import random
import pathlib
import argparse


## sa ##
def sa(old_abundance, temp):
    sa_abundance = old_abundance.copy()
    tes = [k for k in sa_abundance.keys() if k != "_noise"]
    n_neighbors = random.randint(1, max(2, int(.1 * len(tes)**temp))) # added .1*

    for idx in random.sample(tes, n_neighbors):
        log_theta = math.log(sa_abundance[idx])
        delta = random.gauss(0, math.sqrt(temp)*.1* min(1, log_theta+ 10))
        sa_abundance[idx] = max(math.exp(log_theta+delta), 1e-100)
   
    sa_abundance["_noise"] = old_abundance["_noise"]
    sum_norm = sum(sa_abundance.values())
    sa_abundance= {k:v/sum_norm for k,v in sa_abundance.items()}
    return sa_abundance


def find_cooling_rate(n_tes, base_rate=.99, min_rate=.9, max_rate=.999):
    scale = math.log(n_tes+1)
    rate = 1 - (1-base_rate)/scale
    return min(max_rate, max(min_rate, rate))


def reduce_temp(temp, cooling_rate):
    new_temp = cooling_rate * temp
    return max(new_temp, 1e-6)


def accept_sa(ll_old, ll_sa,temp):
    if ll_sa - ll_old > 1e-6:
        return True, "ll_sa"
    delta = (ll_sa-ll_old)/temp
    if delta < -700 or ll_old - ll_sa < 1e-5:
        return False, "none"
    return random.random() < math.exp(delta), "accept"


## em ##
def log_likelihood(theta, multimapped_reads, weights, unique_counts):
    ll = 0.0
    log_theta = {k:math.log(v) for k,v in theta.items()}
    for read, tes in multimapped_reads.items():
        xs = [log_theta[te] + weights[read][te] for te in tes]
        xs.append(log_theta["_noise"] + weights[read]["_noise"])
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    for te in theta:
        if te == "_noise":
            continue
        ll += math.log(theta[te] + 1e-9)*unique_counts.get(te,0)
    return ll


def init_abundance(multimapped_reads, all_tes, unique_counts):
    theta = {k:1/len(all_tes) for k in all_tes}
    theta["_noise"] = 1e-2
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, weights):
    frac = {k:{} for k in multimapped_reads}
    log_theta = {k: math.log(v) for k,v in theta.items()}

    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            xs.append(log_theta[te] + weights[read][te])
        xs.append(log_theta["_noise"] + weights[read]["_noise"])

        tes_plus = tes + ["_noise"]
        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes_plus, xs):
            frac[read][te] = math.exp(x - m) / Z
    return frac


def m_step(frac, multimapped_reads, all_tes, unique_counts):
    theta = {k: unique_counts.get(k, 0) for k in all_tes} ## .1 
    theta["_noise"] = 1e-2
    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te] + .1
        theta["_noise"] += frac[read]["_noise"] + 1
    theta = {k:max(v, 1e-9) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def theta_to_counts(frac, all_tes, counts_threshold):
    counts = {k: 0 for k in all_tes}
    print(len(frac))
    rescue = []
    for read, tes in frac.items():
        vals = sorted(tes.items(), key = lambda x: x[1], reverse = True)
        if len(vals) == 1 or vals[0][1] / (vals[1][1] + 1e-9) > counts_threshold:
            counts[vals[0][0]] += 1
        else:
            rescue += [read]
    print(sum(counts.values()))
    return set(rescue), counts


def get_best(new_abundance, best_abundance, ll_new, ll_best):
    if ll_new > ll_best:
        return new_abundance, ll_new
    return best_abundance, ll_best


def em(multimapped_reads, unique_counts, weights, counts_threshold=1.05):
    threshold = 1e-4
    temp = 1.0
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    all_tes.append("_noise")
    old_abundance = init_abundance(multimapped_reads, all_tes, unique_counts)
    cooling_rate = find_cooling_rate(len(all_tes))
    best_abundance = old_abundance
    ll_old = log_likelihood(old_abundance, multimapped_reads, weights, unique_counts)
    ll_best = ll_old
    accepted = 0

    for i in range(1,10000):
        if i % 50 == 0:
            rate = accepted / 50
            accepted = 0
        #    print(str(i) + "\t" + str(rate) + "\t" + str(ll_old) + "\t" + str(temp))
        sa_abundance = sa(old_abundance, temp)
        ll_sa = log_likelihood(sa_abundance, multimapped_reads, weights, unique_counts)
        accept, lvl = accept_sa(ll_old, ll_sa, temp)
        
        if accept and temp > .05:
            old_abundance = sa_abundance
            ll_old = ll_sa
            best_abundance, ll_best = get_best(sa_abundance, best_abundance, ll_sa, ll_best)
            temp = reduce_temp(temp, cooling_rate)
            if lvl == "accept":
                accepted += 1
            continue
        elif temp < .05:
            old_abundance, ll_old = get_best(best_abundance, old_abundance, ll_best, ll_old)
            return theta_to_counts(e_step(best_abundance, multimapped_reads, weights), all_tes, counts_threshold)
            
        frac = e_step(old_abundance, multimapped_reads, weights)
        new_abundance = m_step(frac, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, multimapped_reads, weights, unique_counts)
        diff = ll_new - ll_old
        
        if abs(diff) < threshold and max(abs(new_abundance[k] - old_abundance[k]) for k in new_abundance) < 1e-4:
            return theta_to_counts(e_step(new_abundance, multimapped_reads, weights), all_tes, counts_threshold)
        temp = reduce_temp(temp, cooling_rate)
        old_abundance = new_abundance
        ll_old = ll_new
    return theta_to_counts(e_step(new_abundance, multimapped_reads, weights), all_tes, counts_threshold)


## LL model weight setup ##
def gc_content(seq):
    return (seq.count("g") + seq.count("G") + seq.count("c") + seq.count("C")) / max(1, len(seq))


def gc_bin(gc, bin_size):
    return round(gc / bin_size) * bin_size


def calc_gc_bias(unique_seqs, bin_size = 0.01, pseudo_count = 1):
    counts = {}
    for seq in unique_seqs.values():
        gc = gc_bin(gc_content(seq), bin_size)
        counts[gc] = counts.get(gc, 0) + 1
    for gc in counts:
        counts[gc] += pseudo_count
    mean = sum(counts.values()) / len(counts)
    return {k: math.log(v/mean) for k,v in counts.items()}


def build_mm_bias(multimapped_reads, transcript_gc, gc_bias):
    bias = {}
    for read, tes in multimapped_reads.items():
        bias[read] = {}

        for te in tes:
            g = transcript_gc[te]
            bin_idx = min(gc_bias, key= lambda x: abs(x-g))
            bias[read][te] = gc_bias[bin_idx]
        bias[read]["_noise"] = min(bias[read].values()) -1
    return bias


def norm_align_scores(align_scores):
    align_weight = {}
    for read, tes in align_scores.items():
        max_score = max([int(x) for x in tes.values()])
        weights = {te: math.log(int(asv)+1) - math.log(max_score + 1) for te, asv in tes.items()}
#        weights = {te: math.log(int(asv)+1) for te, asv in tes.items()} ## added max score
        weights["_noise"] = min(weights.values()) - .5
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



def combine_weights(multimapped_reads, align_scores, gc_bias, e_lens):
    weights = {}
    for read, tes in multimapped_reads.items():
        weights[read] = {}
        for te in tes:
            weights[read][te] = (align_scores[read][te] + gc_bias[read][te] - e_lens[read][te])
        weights[read]["_noise"] = (align_scores[read]["_noise"] + gc_bias[read]["_noise"] - e_lens[read]["_noise"])
    return weights


## load data ##
def parse_gtf(loc):
    gc, len_transcripts = {}, {}
    with open(loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split()
            name = buff[11][1:-2]
            len_transcripts[name] = abs(int(buff[4]) - int(buff[3]))
            gc[name] = float(buff[-1][1:-2])

    return gc, len_transcripts


def parse_unique(loc):
    unique_counts = {}
    with open(loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split(",")
            unique_counts[buff[0]] = int(buff[1])
    return unique_counts


def parse_multimapped(loc):
    multimapped_reads, align_scores = {}, {}
    with open(loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split(",")
            name = buff[0]
            tes = [x.split("*")[1] for x in buff[1:]]
            scores = [x.split("*")[0] for x in buff[1:]]
            multimapped_reads[name] = list(set(tes))
            align_scores[name] = {te:scores[e] for e, te in enumerate(tes)}
    return multimapped_reads, align_scores


def parse_sam(loc):
    unique_seq, read_lens = {}, {}
    with open(loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "@":
                continue
            buff = line.strip().split()
            name = buff[0].split("/")[0]
            read_lens[name] = len(buff[9])
            nh = next((s for s in buff if "NH" in s), None)
            if nh != None and int(nh[-1]) >= 1:
                unique_seq[name] = buff[9]
    return unique_seq, read_lens

def parse_args():
    parser = argparse.ArgumentParser(prog="TE-SAEM: SAEM",
             description = "Simulated annealing based expectation maximization",
             epilog="For help//queries please visit our github: https://github.com/synanth/TE-SAEM")
    parser.add_argument("-d",
                        type = pathlib.Path,
                        required = True,
                        metavar = "data",
                        dest = "data", 
                        help = "Input parsed SAM data location (default %(default)s)")
    parser.add_argument("-g",
                        type = pathlib.Path,
                        required = True,
                        metavar = "gtf",
                        dest = "gtf",
                        help = "TE GTF location (default %(default)s)")
    parser.add_argument("-o",
                        type = pathlib.Path,
                        required = True,
                        metavar = "out",
                        dest = "out",
                        help = "Counts output location (default %(default)s)")

    args = parser.parse_args()
    return args



## driver fxn ##
if __name__ == '__main__':
    n_run = "wut"

    base_loc = "/home/stexocae/li_lab/saem/"
    
    print("Prepping data for SAEM")
    args = parse_args()
    gc, len_transcripts = parse_gtf(args.gtf)
    unique_counts = parse_unique(str(args.data) + "/unique.counts") 
    multimapped_reads, align_scores = parse_multimapped(str(args.data) + "/multi.translation")
    unique_seq, read_lens = parse_sam(str(args.data) + "/star.sam")
    
    print("Calculating weights for log-likelihood")

    align_scores = norm_align_scores(align_scores)
    bias = calc_gc_bias(unique_seq)
    gc_bias = build_mm_bias(multimapped_reads, gc, bias)
    e_lens = calc_e_lens(multimapped_reads, len_transcripts, read_lens)
    weights = combine_weights(multimapped_reads, align_scores, gc_bias, e_lens)

    rescue_counts, em_counts = {}, {}
    print("Running initial SAEM")

    rescue, em_counts = em(multimapped_reads, unique_counts, weights)
    print(str(sum(em_counts.values())) + " reads classified.")
    if "_noise" in em_counts:
        del em_counts["_noise"]
    rescue_reads = {k:v for k,v in multimapped_reads.items() if k in rescue}
    print("Running rescue SAEM.")

    toss, rescue_counts = em(rescue_reads, unique_counts, weights)
    print(str(sum(rescue_counts.values())) + " reads rescued.")

    all_tes = {k:unique_counts.get(k, 0) + em_counts.get(k, 0) + rescue_counts.get(k,0) for k in len_transcripts}


    print("Writing count table to: " + str(args.out))
    with open(args.out, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
