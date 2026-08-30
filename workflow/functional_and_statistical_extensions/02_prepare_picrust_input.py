import sys
counts,fasta,out=sys.argv[1:]
ids=set()
with open(counts,encoding='utf-8') as f:
    next(f)
    for line in f:
        ids.add(line.split('\t',1)[0])
keep=False
with open(fasta,encoding='utf-8') as fi, open(out,'w',encoding='utf-8') as fo:
    for line in fi:
        if line.startswith('>'):
            keep=line[1:].split()[0] in ids
        if keep: fo.write(line)

