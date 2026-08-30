#!/usr/bin/env python3
import argparse, csv, gzip, math, re
from collections import Counter
from pathlib import Path

ROOT=Path('/path/to/morchella_microbiome/public_validation')
PRIMERS={
 'PRJNA935967': [('ITS1F','CTTGGTCATTTAGAGGAAGTAA'),('ITS2R','GCTGCGTTCTTCATCGATGC')],
 'PRJNA993383': [('515F','GTGCCAGCMGCCGCGGTAA'),('806R','GGACTACHVGGGTWTCTAAT')]
}
IUPAC={'A':'A','C':'C','G':'G','T':'T','M':'AC','R':'AG','W':'AT','S':'CG','Y':'CT','K':'GT','V':'ACG','H':'ACT','D':'AGT','B':'CGT','N':'ACGT'}
COMP=str.maketrans('ACGTRYMKBDHVN','TGCAYRKMVHDBN')
def rc(s): return s.translate(COMP)[::-1]
def regex(s): return re.compile(''.join('['+IUPAC.get(x,x)+']' for x in s))
def quantile(hist, q):
    n=sum(hist.values()); target=max(1, math.ceil(n*q)); seen=0
    for k in sorted(hist):
        seen += hist[k]
        if seen >= target: return k
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project',required=True); a=ap.parse_args()
    indir=ROOT/'projects'/a.project/'raw'/'fastq'; outdir=ROOT/'projects'/a.project/'qc'; outdir.mkdir(exist_ok=True)
    patterns=[]
    for name,p in PRIMERS[a.project]: patterns += [(name,regex(p)),(name+'_RC',regex(rc(p)))]
    rows=[]; pos_rows=[]
    for fq in sorted(indir.glob('*.fastq.gz')):
        lh=Counter(); reads=nbases=q20=q30=nreads=0; qsum=[]; qcount=[]; hits=Counter(); malformed=0
        with gzip.open(fq,'rt',encoding='ascii',errors='replace') as h:
            while True:
                rec=[h.readline() for _ in range(4)]
                if not rec[0]: break
                if any(x=='' for x in rec) or not rec[0].startswith('@') or not rec[2].startswith('+'):
                    malformed+=1; continue
                seq=rec[1].strip().upper(); qual=rec[3].strip()
                if len(seq)!=len(qual): malformed+=1; continue
                reads+=1; lh[len(seq)]+=1; nbases+=len(seq); nreads += ('N' in seq)
                if len(qsum)<len(qual): qsum.extend([0]*(len(qual)-len(qsum))); qcount.extend([0]*(len(qual)-len(qcount)))
                for i,c in enumerate(qual):
                    q=ord(c)-33; qsum[i]+=q; qcount[i]+=1; q20+=(q>=20); q30+=(q>=30)
                edge=seq[:40]
                for name,pat in patterns: hits[name]+=bool(pat.search(edge))
        run=fq.name.split('.')[0]
        row={'project':a.project,'run':run,'file':str(fq),'reads':reads,'bases':nbases,'malformed':malformed,
             'min_len':min(lh) if lh else '', 'q25_len':quantile(lh,.25) if lh else '', 'median_len':quantile(lh,.5) if lh else '',
             'q75_len':quantile(lh,.75) if lh else '', 'max_len':max(lh) if lh else '', 'reads_with_N':nreads,
             'base_Q20_fraction':round(q20/nbases,6) if nbases else '', 'base_Q30_fraction':round(q30/nbases,6) if nbases else ''}
        row.update({'start_hit_'+k:v for k,v in hits.items()}); rows.append(row)
        for i,(s,c) in enumerate(zip(qsum,qcount),1): pos_rows.append({'project':a.project,'run':run,'position':i,'mean_q':round(s/c,4),'bases_observed':c})
        with (outdir/f'{run}_length_hist.tsv').open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(['length','reads']); w.writerows(sorted(lh.items()))
    with (outdir/'full_fastq_qc_summary.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    with (outdir/'quality_by_position.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(pos_rows[0]),delimiter='\t'); w.writeheader(); w.writerows(pos_rows)
    with (outdir/'full_fastq_qc_versions.txt').open('w',encoding='utf-8') as f:
        import sys,platform; f.write(sys.version+'\n'+platform.platform()+'\n')
if __name__=='__main__': main()


