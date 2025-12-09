import sys
import math
import random

## sa ##
def sa(old_abundance, temp, step_size, neighborhood):
    sa_abundance = old_abundance.copy()
    if temp == 0:
        n_neighbors = 2
    else:
        n_neighbors = round(len(old_abundance.keys()) * (neighborhood))
    tes = list(sa_abundance.keys())
    for n in range(n_neighbors):
        idx = tes[random.randint(0, len(tes)-1)]
        change = random.uniform(-step_size, step_size)
        sa_abundance[idx] += change
        sa_abundance[idx] = max(sa_abundance[idx], 1e-300)
    sa_abundance = {k:v/sum(sa_abundance.values()) for k,v in sa_abundance.items()}
    return sa_abundance


def reduce_temp(i, temp, cooling_rate, cooling_schedule):
    if cooling_schedule < .25:
        return reduce_temp_linear(i, temp, cooling_rate)
    elif cooling_schedule < .5:
        return reduce_temp_log(i, temp, cooling_rate)
    elif cooling_schedule < .75:
        return reduce_temp_exp(i, temp, cooling_rate)
    return reduce_temp_mixed(i, temp, cooling_rate)

def reduce_temp_linear(i, temp, cooling_rate):
    temp -= cooling_rate
    return max(temp, 0.01)


def reduce_temp_log(i, temp, cooling_rate):
    temp /= math.log(1 + i * cooling_rate)
    return max(temp, 0.01)


def reduce_temp_exp(i, temp, cooling_rate):
    temp *= (cooling_rate*temp)
    return max(temp, 0.01)


def reduce_temp_mixed(i, temp, cooling_rate):
    if i < 20:
        temp *= (1-cooling_rate)
    elif i < 50:
        temp -= cooling_rate
    else:
        temp *= cooling_rate
    return max(temp, 0.01)



def accept_sa(ll_old, ll_sa, temp):
    accept = random.random()
#    print(temp, ll_sa - ll_old, math.exp((abs(ll_sa - ll_old)/temp))
    if ll_sa > ll_old:
        return True
    elif accept < math.exp((ll_sa + ll_old)/temp):
    #    print(accept, math.exp((ll_sa - ll_old)/temp))
        return True
    return False

## em ##

def log_likelihood(theta, len_transcripts, multimapped_reads, read_lens, avg_len):
    log_sum = 0
    for read, tes in multimapped_reads.items():
        if sum([theta[te] / max(len_transcripts[te]-read_lens[read]+1,1) for te in tes]) > 0:
            log_sum += math.log(sum([theta[te] / max(len_transcripts[te]-read_lens[read]+1,1) for te in tes]) + 1e-300)
    return log_sum


def init_abundance(len_transcripts, multimapped_reads):
    tes = set([x for sublist in multimapped_reads.values() for x in sublist])
    theta = {k:0 for k in tes}
    
    for te in tes:
        theta[te] = random.random()
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, len_transcripts, read_lens, avg_len):
    frac = {k:0 for k in theta}

    for read, tes in multimapped_reads.items():
        denom = sum([theta[te]/max(len_transcripts[te]-read_lens[read]+1,1) for te in tes]) + 1e-12
        for te in tes:
            e_len = max(len_transcripts[te]-read_lens[read]+1,1)
            frac[te] += (theta[te] / e_len) / denom
    return frac


def m_step(frac, len_transcripts, read_lens, multimapped_reads, avg_len):
    theta = {k: 0 for k,v in frac.items()}
    denom = sum([frac[te] /max(len_transcripts[te]-avg_len+1,1) for te in theta]) + 1e-12

    for te in theta:
        theta[te] = (frac[te]/max(len_transcripts[te]-avg_len+1,1) + 1e-12) / denom
    return theta


def em(len_transcripts, read_lens, multimapped_reads, cooling_rate, cooling_schedule, step_size, neighborhood, acceptence_prob):

    threshold = 0.01
    temp = 1
    old_abundance = init_abundance(len_transcripts, multimapped_reads)
    best = old_abundance
    avg_len = sum(read_lens.values())/len(read_lens)

#    for i in range(1,10000):
    for i in range(1,1000):
        sa_abundance = sa(old_abundance, temp, step_size, neighborhood)
        ll_old = log_likelihood(old_abundance, len_transcripts, multimapped_reads, read_lens, avg_len)
        ll_sa = log_likelihood(sa_abundance, len_transcripts, multimapped_reads, read_lens, avg_len)
#        print(ll_sa, ll_old, ll_sa + ll_old)        
        accept = accept_sa(ll_old, ll_sa, temp)
        if accept:
            old_abundance = sa_abundance
            diff = ll_sa + ll_old
    #        print(i, ll_sa - ll_old, ll_sa, ll_old, accept)
            continue
        frac = e_step(old_abundance, multimapped_reads, len_transcripts, read_lens, avg_len)
        new_abundance = m_step(frac, len_transcripts, read_lens, multimapped_reads, avg_len)
        ll_new = log_likelihood(new_abundance, len_transcripts, multimapped_reads, read_lens, avg_len)
        diff = ll_new - ll_old
        old_abundance = new_abundance
     #   print(i, diff, ll_old, ll_new)

        if abs(diff) < threshold:
            return old_abundance
        temp = reduce_temp(i, temp, cooling_rate, cooling_schedule)
    return old_abundance



## driver fxn ##
if __name__ == '__main__':

    cooling_rate = float(sys.argv[1])
    cooling_schedule = float(sys.argv[2])
    step_size = float(sys.argv[3])
    neighborhood = float(sys.argv[4])
    acceptence_prob = float(sys.argv[5])
    n_run = sys.argv[6]

    base_loc = "/home/stexocae/li_lab/saem/"
    gtf_loc = base_loc + "refs/hs1.gtf"
    unique_counts_loc = base_loc + "sim_data/star.unique.counts"
    multimapped_loc = base_loc + "sim_data/star.multi.translation"
    sam_loc = base_loc + "sim_data/alignment/segemehl.sam"

    readlen = 150

    len_transcripts = {}
    unique_counts = {}
    multi_counts = {}
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
    

    em_counts = em(len_transcripts, read_lens, multimapped_reads, cooling_rate, cooling_schedule, step_size, neighborhood, acceptence_prob)
    non_zero = {k:int(v*len(multimapped_reads)) for k,v in em_counts.items()}

    all_tes = {k:0 for k in len_transcripts}
    for te in all_tes:
        if te in unique_counts:
            all_tes[te] += unique_counts[te]
        if te in non_zero:
            all_tes[te] += non_zero[te]

    out_loc = base_loc + "manuscript/out/saem_" + n_run + ".out"
    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
