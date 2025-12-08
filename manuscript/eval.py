import sys

n_run = sys.argv[1]

base_loc = "/home/stexocae/li_lab/saem/"
sim_loc = base_loc + "sim_data/true_counts.csv"
saem_loc = base_loc + "manuscript/out/saem_" + n_run + ".out"

true_counts = {}
saem_counts = {}


with open(sim_loc, "r") as f:
    lines = f.readlines()
    for line in lines[1:]:
        buff = line.strip().split(",")
        true_counts[buff[0]] = int(buff[1])

unrepresented = 0
with open(saem_loc, "r") as f:
    lines = f.readlines()
    for line in lines:
        buff = line.strip().split(",")
        
        if buff[0] not in true_counts.keys():
            unrepresented += int(buff[1])
        else:
            saem_counts[buff[0]] = int(buff[1])
    saem_counts["unrepresented"] = unrepresented
all_tes = set(list(saem_counts.keys()) + list(true_counts.keys()))


tp, fn, fp = 0, 0, 0
for te in all_tes:
    if te not in saem_counts.keys():
        fn += true_counts[te]
        continue
    if te not in true_counts.keys():
        fp += saem_counts[te]
        continue
    if true_counts[te] - saem_counts[te] >= 0:
        tp += saem_counts[te]
        fn += true_counts[te] - saem_counts[te]
    else:
        tp += true_counts[te]
        fp += saem_counts[te] - true_counts[te]


print("tp: " + str(tp))
print("fn: " +  str(fn))
print("fp: " + str(fp))
print("sensitivity: " + str(tp/(tp+fn)))
print("precision: " + str(tp/(tp+fp)))
