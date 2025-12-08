#import multiprocessing as mp
import subprocess
import random
import os


out_loc = "/home/stexocae/li_lab/saem/manuscript/out/progress/"
top_loc = out_loc + "top.txt"
seed_loc = out_loc + "seed.csv"


def init_population(n):
    pop = []
    for i in range(n):
        pop += [[random.random(), random.random(), random.random(), random.random(), random.random()]]
    return pop

def new_pop_from_top(top, mutate_rate, n):
    pop = []
    for i in range(n):
        parents = random.sample(top, 2)
        child = breed(parents[0], parents[1])
        child = mutate(child, mutate_rate)
        pop += [child]
    return pop

def breed(mother, father):
    child = []
    for x in range(5):
        if random.random() < .5:
            child += [mother[x]]
        else:
            child += [father[x]]
    return child

def mutate(params, rate):
    for x in range(5):
        if random.random() < rate:
            params[x] = random.random()
    return params


def run_saem(cooling_rate, cooling_schedule, step_size, neighborhood, acceptence_prob, n_run):
    saem_call = "python3 " + saem_loc + " " + str(cooling_rate) + " " + str(cooling_schedule) + " " + str(step_size) + " " + str(neighborhood) + " " + str(acceptence_prob) + " " + str(n_run)
    subprocess.run(saem_call, shell=True)

    eval_call = "python3 " + eval_loc + " " + str(n_run)
    output = subprocess.run(eval_call, shell=True, check=True, stdout=subprocess.PIPE).stdout.decode("utf-8").strip().split("\n")
    output = [x.split()[1] for x in output]
    results = [cooling_rate, cooling_schedule, step_size, neighborhood, acceptence_prob] + [float(x) for x in output]
    return results

def load_data():
    if not os.path.exists(top_loc):
        return (0, init_population(n_pop))
    with open(top_loc, "r") as f:
        run = int(f.readlines()[-2].strip().split(",")[0])
    with open(seed_loc, "r") as f:
        lines = [[float(y) for y in x.strip().split(",")] for x in f.readlines()]
    return (run, lines)



if __name__ == '__main__':
    base_loc = "/home/stexocae/li_lab/saem/scripts/"
    saem_loc = base_loc  + "saem.py"
    eval_loc = base_loc + "eval.py"

    n_generations = 1000
    n_best = 50
    n_pop = 100
    mutation_rate = 0.1

    results = [["cooling_rate", "cooling_schedule", "step_size", "neighborhood", "acceptence_prob", 
            "tp", "fn", "fp", "sensitivity", "precision"]]
    progress = []

    i, population = load_data()

    while i < n_generations:
        res = []
        for e, x in enumerate(population):
            res += [run_saem(x[0], x[1], x[2], x[3],x[4], e)]
        res = [r + [2*(r[-1] * r[-2])/(r[-1] + r[-2])] for r in res]
        top = sorted(res, reverse=True, key=lambda x:x[-3])[:n_best]
        best = top[0][:5]
        with open(top_loc, "a") as f:
            f.write(str(i) + ',{:.6}'.format(top[0][-3]) + ',{:.6}'.format(top[0][-2]) + ',{:.6}'.format(top[0][-1]) + "\n")
            f.write("best: " + ",".join([str(x) for x in best]) + "\n")

        with open(seed_loc, "w") as f:
            for t in top:
                f.write(",".join([str(x) for x in t[:5]]) + "\n")
        
        progress += [top[0][-1]]
        population = new_pop_from_top(top, mutation_rate, n_pop) + [best]
        i += 1
    print()
    print("Fin, best:")
    print(best)
