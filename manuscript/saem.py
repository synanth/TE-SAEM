import sys
import math
import random

## sa ##
def sa(old_abundance, temp, step_size, neighborhood):
    sa_abundance = old_abundance.copy()
    n_neighbors = max(1, round(len(sa_abundance)*temp*neighborhood))
    tes = list(sa_abundance.keys())
    for idx in random.sample(tes, n_neighbors):
        change = random.uniform(-step_size, step_size)
        sa_abundance[idx] = max(sa_abundance[idx] +change, 1e-100)
    sum_norm = sum(sa_abundance.values())
    sa_abundance= {k:v/sum_norm for k,v in sa_abundance.items()}
    return sa_abundance


def reduce_temp(temp, cooling_rate):
    ## testing geometric cooling ##
    new_temp = cooling_rate * temp
    return new_temp 


def accept_sa(ll_old, ll_sa, temp, n):
    accept = max(random.random(), 1e-16)
    if ll_sa > ll_old:
        return True
    deltaE = (ll_sa -ll_old)/ll_old * temp
    deltaE = math.exp((ll_sa-ll_old)/(temp*n))
    if accept < deltaE:
        print(deltaE*(temp), accept)
        return True
    return False


## em ##
def log_likelihood(theta, len_transcripts, multimapped_reads, read_lens):
    ll = 0.0
    for read, tes in multimapped_reads.items():
        xs = [
            math.log(theta[te]) -
            math.log(max(len_transcripts[te] - read_lens[read] + 1, 1))
            for te in tes
        ]
        m = max(xs)
        ll += m + math.log(sum(math.exp(x - m) for x in xs))
    return ll


def init_abundance(len_transcripts, multimapped_reads, all_tes):
    theta = {k:0 for k in all_tes}

    for te in all_tes:
        theta[te] = max(random.random(), 1e-300)
#        theta[te] = (1/max(len_transcripts[te]-avg_len+1,1)+1e-12) 
    means_sum = sum(theta.values()) 
    theta = {k:(v/means_sum) for k,v in theta.items()}
    return theta


def e_step(theta, multimapped_reads, len_transcripts, read_lens):
    frac = {k:{} for k in multimapped_reads}

    for read, tes in multimapped_reads.items():
        xs = []
        for te in tes:
            e_len = max(len_transcripts[te] - read_lens[read] + 1, 1)
            xs.append(math.log(theta[te]) - math.log(e_len))

        m = max(xs)
        Z = sum(math.exp(x - m) for x in xs)

        for te, x in zip(tes, xs):
            frac[read][te] = math.exp(x - m) / Z

    return frac


def m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes):
    theta = {k: 0 for k in all_tes}

    for read, tes in multimapped_reads.items():
        for te in tes:
            theta[te] += frac[read][te]
    theta = {k:max(v, 1e-12) for k,v in theta.items()}
    theta_sum = sum(theta.values())
    theta = {k:v/theta_sum for k,v in theta.items()}
    return theta


def em(len_transcripts, read_lens, multimapped_reads, unique_counts, cooling_rate, step_size, neighborhood):

    threshold = 0.01
    temp = 1.0
    all_tes = list(set([x for sublist in multimapped_reads.values() for x in sublist]))
    old_abundance = init_abundance(len_transcripts, multimapped_reads, all_tes)

    #for i in range(1,10000):
    for i in range(1,5000):
        sa_abundance = sa(old_abundance, temp, step_size, neighborhood)
        ll_sa = log_likelihood(sa_abundance, len_transcripts, multimapped_reads, read_lens)
        ll_old = log_likelihood(old_abundance, len_transcripts, multimapped_reads, read_lens)
        if accept_sa(ll_old, ll_sa, temp, len(multimapped_reads)):
 #           print(i, ll_old, ll_sa, temp)
            old_abundance = sa_abundance
            temp = reduce_temp(temp, cooling_rate)
            continue
        frac = e_step(old_abundance, multimapped_reads, len_transcripts, read_lens)
        new_abundance = m_step(frac, len_transcripts, read_lens, multimapped_reads, all_tes)
        ll_new = log_likelihood(new_abundance, len_transcripts, multimapped_reads, read_lens)
        diff = ll_new - ll_old
#        print(i, diff, ll_old, ll_new)
        if abs(diff) < threshold:
            return {k:v*len(multimapped_reads) for k,v in old_abundance.items()}
        temp = reduce_temp(temp, cooling_rate)
        old_abundance = new_abundance
    print(old_abundance)
    return {k:v*len(multimapped_reads) for k,v in old_abundance.items()}


## driver fxn ##
if __name__ == '__main__':

    cooling_rate = float(sys.argv[1])
    step_size = float(sys.argv[2])
    neighborhood = float(sys.argv[3])
    n_run = sys.argv[4]

    base_loc = "/home/stexocae/li_lab/saem/"
    gtf_loc = base_loc + "refs/hs1.gtf"
    unique_counts_loc = base_loc + "sim_data/star.unique.counts"
    multimapped_loc = base_loc + "sim_data/star.multi.translation"
    sam_loc = base_loc + "sim_data/alignment/star.sam"
    assembly_loc = base_loc + "sim_data/assembly/to_contigs.sam"

    len_transcripts = {}
    unique_counts = {}
    multimapped_reads = {}
    read_lens = {}

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
    

    em_counts = em(len_transcripts, read_lens, multimapped_reads, unique_counts, cooling_rate, step_size, neighborhood)
    non_zero = {k:round(v) for k,v in em_counts.items()}

    all_tes = {k:0 for k in len_transcripts}
    for te in all_tes:
        if te in unique_counts:
            all_tes[te] += unique_counts[te]
        if te in non_zero.keys():
            all_tes[te] += non_zero[te]

    out_loc = base_loc + "manuscript/out/saem_" + n_run + ".out"

    with open(out_loc, "w") as f:
        for k,v in all_tes.items():
            f.write(k + "," + str(v) + "\n")
