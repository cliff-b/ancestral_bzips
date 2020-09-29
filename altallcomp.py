
mlseqs = {}
with open("/Users/Cliff/Documents/Kosuri Lab/Design docs/181114 Ancestral tree/ref_sequences/Ancestors.txt") as f:
    for line in f:
        seq, name = line.rstrip().split(sep = '\t')
        if ',' in name:
            for i in name.split(sep = ','):
                mlseqs[i] = seq
        else:
            mlseqs[name] = seq

altallseqs =  {}
with open("/Users/Cliff/Documents/Kosuri Lab/Design docs/181114 Ancestral tree/ref_sequences/Altalls.txt") as f:
    for line in f:
        seq, name = line.rstrip().split(sep = '\t')
        if ',' in name:
            for i in name.split(sep = ','):
                altallseqs[i] = seq
        else:
            altallseqs[name] = seq

difflist = []
diffdic = {}
for i in altallseqs.keys():
    diff = 0
    if i in mlseqs.keys():
        for c, j in enumerate(altallseqs[i]):
            if j != mlseqs[i][c]:
                diff += 1
        difflist.append(diff)
        diffdic[i] = [altallseqs[i], mlseqs[i], diff]

print(difflist)
print("number of different sequences", len(difflist))
print("average number of differences", sum(difflist)/len(difflist))
print("holdhere for diffdic")