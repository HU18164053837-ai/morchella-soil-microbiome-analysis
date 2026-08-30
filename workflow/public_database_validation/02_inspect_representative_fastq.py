#!/usr/bin/env python3
import csv, gzip, hashlib, statistics
from pathlib import Path

ROOT = Path('/path/to/morchella_microbiome/public_validation')
OUT = ROOT / 'qc' / 'representative_structure_qc.tsv'
REPORT = ROOT / 'qc' / 'representative_structure_report_zh.md'
OUT.parent.mkdir(parents=True, exist_ok=True)

PRIMERS = {
    'ITS1_F': 'CTTGGTCATTTAGAGGAAGTAA', 'ITS2_R': 'GCTGCGTTCTTCATCGATGC',
    '515F': 'GTGCCAGCMGCCGCGGTAA', '806R_V4': 'GGACTACHVGGGTWTCTAAT',
    '341F': 'CCTAYGGGRBGCASCAG', '806R_V34': 'GGACTACNNGGGTATCTAAT',
    'ITS3': 'GCATCGATGAAGAACGCAGC', 'ITS4': 'TCCTCCGCTTATTGATATGC'
}

IUPAC = {'A':'A','C':'C','G':'G','T':'T','M':'AC','R':'AG','W':'AT','S':'CG','Y':'CT','K':'GT','V':'ACG','H':'ACT','D':'AGT','B':'CGT','N':'ACGT'}
COMP = str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN')
def rc(s): return s.translate(COMP)[::-1]
def hit(seq, primer):
    if len(seq) < len(primer): return False
    for offset in range(0, min(12, max(1, len(seq)-len(primer)+1))):
        sub=seq[offset:offset+len(primer)]
        if sum(b in IUPAC.get(p,p) for b,p in zip(sub,primer)) >= len(primer)-2: return True
    return False

rows=[]
for fq in sorted(ROOT.glob('projects/*/raw/representative/*.fastq.gz')):
    lens=[]; qmeans=[]; headers=[]; hits={k:0 for k in PRIMERS}; hits.update({k+'_RC':0 for k in PRIMERS})
    malformed=0; n=0
    with gzip.open(fq, 'rt', encoding='ascii', errors='replace') as h:
        while n < 2000:
            rec=[h.readline() for _ in range(4)]
            if not rec[0]: break
            if any(x == '' for x in rec) or not rec[0].startswith('@') or not rec[2].startswith('+') or len(rec[1].strip()) != len(rec[3].strip()):
                malformed += 1; continue
            seq=rec[1].strip().upper(); qual=rec[3].strip()
            n += 1; lens.append(len(seq)); qmeans.append(sum(ord(c)-33 for c in qual)/len(qual)); headers.append(rec[0].strip())
            for k,p in PRIMERS.items():
                hits[k] += hit(seq,p); hits[k+'_RC'] += hit(seq,rc(p))
    name=fq.name; run=name.split('_')[0].split('.')[0]
    mate='R1' if '_1.fastq' in name else ('R2' if '_2.fastq' in name else 'single')
    row={'project':fq.parts[-4], 'run':run, 'file':str(fq), 'mate':mate, 'reads_inspected':n,
         'malformed_records':malformed, 'min_len':min(lens) if lens else '', 'median_len':statistics.median(lens) if lens else '',
         'max_len':max(lens) if lens else '', 'mean_q':round(statistics.mean(qmeans),2) if qmeans else '',
         'header_example':headers[0] if headers else ''}
    row.update({f'hit_{k}':v for k,v in hits.items()}); rows.append(row)

fields=list(rows[0]) if rows else []
with OUT.open('w', encoding='utf-8', newline='') as f:
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(rows)

with REPORT.open('w', encoding='utf-8') as f:
    f.write('# 羊肚菌公共数据库外部验证：代表性 FASTQ 结构检查\n\n')
    f.write('每个 FASTQ 最多读取前 2,000 条，仅用于结构、长度、质量和引物签名诊断，不属于正式去噪。\n\n')
    for r in rows:
        ranked=sorted(((k.replace('hit_',''),v) for k,v in r.items() if k.startswith('hit_')), key=lambda x:x[1], reverse=True)[:4]
        f.write(f"- {r['project']} {r['run']} {r['mate']}: n={r['reads_inspected']}, length={r['min_len']}/{r['median_len']}/{r['max_len']}, meanQ={r['mean_q']}, malformed={r['malformed_records']}; top primer signatures: {', '.join(f'{k}={v}' for k,v in ranked)}.\n")

