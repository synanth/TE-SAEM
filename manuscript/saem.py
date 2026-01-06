import sys
import math
import random


## sa ##
def sa(old_abundance, temp, unique_counts):
    sa_abundance = old_abundance.copy()
    tes = [k for k in sa_abundance.keys() if unique_counts.get(k, 0) == 0]
    n_neighbors = 1
    
    for idx in random.sample(tes, n_neighbors):
        log_theta = math.log(sa_abundance[idx])
        delta = random.gauss(0, math.sqrt(temp) / len(tes))
        sa_abundance[idx] = max(math.exp(log_theta+delta), 1e-100)
   
    sum_norm = sum(sa_abundance.values())
    sa_abundance= {k:v/sum_norm for k,v in sa_abundance.items()}
    return sa_abundance


def find_cooling_rate(n_tes, base_rate=.99, min_rate=.9, max_rate=.999):
    scale = math.log(n_tes+1)
    rate = 1 - (1-base_rate)/scale
    return min(max_rate, max(min_rate, rate))


def reduce_temp(temp, cooling_rate):
    ## testing geometric cooling ##
    new_temp = cooling_rate * temp
    return max(new_temp, 1e-6)


def accept_sa(ll_old, ll_sa,temp):
    if ll_sa >= ll_old:
        return True, "ll_sa"
    delta = (ll_sa-ll_old)/temp
    if delta < -700:
        return False, "none"
    return random.random() < math.exp(delta), "accept"


## em ##
def log_likelihood(theta, multimapped_reads, align_scores, e_lens, unique_counts, all_tes):
    ll = 0.0
    for read, tes in multimapped_reads.items():
        xs = [
            math.log(theta[te])
            + align_scores[read][te]
            - e_lens[read][te]
            for te in tes
        ]
        
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    for te in all_tes:
        ll += unique_counts.get(te,0) * math.log(max(theta[te], 1e-300))
    return ll


def init_abundance(multimapped_reads, all_tes):
    theta = {k:1/len(all_tes) for k in all_tes}
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, align_scores, e_lens):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            xs.append(math.log(theta[te]) 
            + align_scores[read][te]
            - e_lens[read][te])
        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes, xs):
            frac[read][te] = math.exp(x - m) / Z
    return frac


def m_step(frac, multimapped_reads, all_tes, unique_counts):
    theta = {k: unique_counts.get(k, 0) for k in all_tes}
    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]
    theta = {k:max(v, 1e-9) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def theta_to_counts(frac, all_tes):
    counts = {k: 0 for k in all_tes}
    rescue = []
    
    for read, tes in frac.items():
        vals = sorted(tes.items(), key = lambda x: x[1], reverse = True)
        if vals[0][1] / (vals[1][1] + 1e-9) > 1.05:
            counts[vals[0][0]] += 1
        else:
            rescue += [read]
    return set(rescue), counts


def get_best(new_abundance, best_abundance, ll_new, ll_best):
    if ll_new > ll_best:
        return new_abundance, ll_new
    return best_abundance, ll_best
    

def em(multimapped_reads, unique_counts, align_scores, e_lens):
    threshold = 1e-6
    temp = 1.0
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    old_abundance = init_abundance(multimapped_reads, all_tes)
    cooling_rate = find_cooling_rate(len(all_tes))
    best_abundance = old_abundance
    ll_old = log_likelihood(old_abundance, multimapped_reads, align_scores, e_lens, unique_counts, all_tes)
    ll_best = ll_old

    for i in range(1,10000):
        sa_abundance = sa(old_abundance, temp, unique_counts)
        ll_sa = log_likelihood(sa_abundance, multimapped_reads, align_scores, e_lens, unique_counts, all_tes)
        accept, lvl = accept_sa(ll_old, ll_sa, temp)

        if accept and temp > .05:
            print(str(i) + "\t" + lvl + "\t" + str(ll_sa) +"\t" + str(ll_best) +"\t" + str(temp))
            old_abundance = sa_abundance
            ll_old = ll_sa
            best_abundance, ll_best = get_best(sa_abundance, best_abundance, ll_sa, ll_best)
            temp = reduce_temp(temp, cooling_rate)
            continue
        elif temp < .05:
            old_abundance, ll_old = get_best(old_abundance, best_abundance, ll_old, ll_best)
#            return theta_to_counts(e_step(old_abundance, multimapped_reads, align_scores, e_lens), all_tes)
            
            
        frac = e_step(old_abundance, multimapped_reads, align_scores, e_lens)
        new_abundance = m_step(frac, multimapped_reads, all_tes, unique_counts)
        ll_new = log_likelihood(new_abundance, multimapped_reads, align_scores, e_lens, unique_counts, all_tes)
        diff = ll_new - ll_old
        print(str(i) + "\treject\t" + str(ll_new) + "\t" + str(ll_best) + "\t" + str(temp))
        if abs(diff) < threshold:
            return theta_to_counts(e_step(new_abundance, multimapped_reads, align_scores, e_lens), all_tes)
        temp = reduce_temp(temp, cooling_rate)
        old_abundance = new_abundance
        ll_old = ll_new
    return theta_to_counts(e_step(new_abundance, multimapped_reads, align_scores, e_lens), all_tes)


## LL model setup ##
def norm_align_scores(align_scores):
    align_weight = {}
    for read, tes in align_scores.items():
        max_score = max(tes.values())
        weights = {te: int(asv)/int(max_score) for te, asv in tes.items()}
        weights = {te: math.log(int(asv)+1) for te, asv in tes.items()}
        align_weight[read] = weights
    return align_weight


def calc_e_lens(multimapped_reads, len_transcripts, read_lens):
    e_lens = {k:{} for k in multimapped_reads}
    len_max = max(len_transcripts.values())
    for read, tes in multimapped_reads.items():
        for te in tes:
            e_lens[read][te] = math.log(max(len_transcripts[te] - read_lens[read] + 1, 1))
    return e_lens




## driver fxn ##
if __name__ == '__main__':
    n_run = sys.argv[1]

    base_loc = "/home/stexocae/li_lab/saem/"
    gtf_loc = base_loc + "refs/hs1.gtf"
    unique_counts_loc = base_loc + "sim_data/star.unique.counts"
    multimapped_loc = base_loc + "sim_data/star.multi.translation"
    sam_loc = base_loc + "sim_data/alignment/star.sam"

    len_transcripts = {}
    unique_counts = {}
    unique_seq = {}
    multimapped_reads = {}
    align_scores = {}
    read_lens = {}

    with open(gtf_loc, "r") as f:
        lines = f.readlines()
        for line in lines:
            buff = line.strip().split()
            name = buff[11][1:-2]
            len_transcripts[name] = abs(int(buff[4]) - int(buff[3]))

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

    align_scores = norm_align_scores(align_scores)
    e_lens = calc_e_lens(multimapped_reads, len_transcripts, read_lens)
    rescue_counts, em_counts = {}, {}
    rescue, em_counts = em(multimapped_reads, unique_counts, align_scores, e_lens)


    rescue_reads = {k:v for k,v in multimapped_reads.items() if k in rescue}
    print(len(rescue))
    toss, rescue_counts = em(rescue_reads, unique_counts, align_scores, e_lens)

    all_tes = {k:unique_counts.get(k, 0) + em_counts.get(k, 0) + rescue_counts.get(k,0) for k in len_transcripts}

    out_loc = base_loc + "manuscript/out/saem_" + n_run + ".out"

    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")

