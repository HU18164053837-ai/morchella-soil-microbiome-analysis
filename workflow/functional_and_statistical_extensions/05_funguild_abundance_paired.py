import csv,itertools,math,os,sys
counts_file,guild_file,outdir=sys.argv[1:]
os.makedirs(outdir,exist_ok=True)
def norm_sample(x):
    for suf in ('_ITS_trimmed.fastq.gz','_16S_trimmed.fastq.gz'):
        if x.endswith(suf): return x[:-len(suf)]
    return x
with open(guild_file,encoding='utf-8') as f:
    ann={r['ASV_ID']:r for r in csv.DictReader(f,delimiter='\t')}
with open(counts_file,encoding='utf-8') as f:
    r=csv.reader(f,delimiter='\t'); hdr=next(r); samples=[norm_sample(x) for x in hdr[1:]]; rows=list(r)
ma=[s for s in samples if s[:1] in ('M','A') and len(s)>=4 and s[1] in '234']; idx=[samples.index(s) for s in ma]
levels={'strict':{'Probable','Highly Probable'},'inclusive':{'Possible','Probable','Highly Probable'}}
def split_modes(x): return [z.strip(' |') for z in (x or '').split('-') if z.strip(' |')]
def split_guilds(x): return [z.strip(' |') for z in (x or '').replace('|','').split('-') if z.strip()]
def bh(p):
    n=len(p); order=sorted(range(n),key=lambda i:p[i]); q=[1.0]*n; cur=1.0
    for rank,i in reversed(list(enumerate(order,1))): cur=min(cur,p[i]*n/rank); q[i]=cur
    return q
def exact_paired(vals):
    pairs=[f'{y}-{p}' for y in (2,3,4) for p in (1,2,3)]; d=[]
    for pair in pairs: d.append(vals.get('A'+pair,0.0)-vals.get('M'+pair,0.0))
    obs=sum(d)/len(d); perm=[]
    for signs in itertools.product((-1,1),repeat=9): perm.append(sum(x*s for x,s in zip(d,signs))/9)
    p=(1+sum(abs(x)>=abs(obs)-1e-15 for x in perm))/(1+len(perm))
    return obs,p,d
all_outputs=[]
for confidence,allowed in levels.items():
  for ontology,splitter,field in [('trophic_mode',split_modes,'trophic_mode'),('guild',split_guilds,'guild')]:
    raw={s:{} for s in ma}; total={s:0.0 for s in ma}; annotated={s:0.0 for s in ma}
    for row in rows:
        a=ann.get(row[0]); vals=[float(row[1+i]) for i in idx]
        for s,v in zip(ma,vals): total[s]+=v
        if not a or a.get('confidence_ranking','') not in allowed: continue
        cats=splitter(a.get(field,''))
        if not cats: continue
        for s,v in zip(ma,vals):
            annotated[s]+=v
            for cat in cats: raw[s][cat]=raw[s].get(cat,0.0)+v/len(cats)
    cats=sorted({c for s in ma for c in raw[s]})
    long=[]; stats=[]
    for cat in cats:
        va={s:(raw[s].get(cat,0.0)/total[s] if total[s] else 0) for s in ma}
        vi={s:(raw[s].get(cat,0.0)/annotated[s] if annotated[s] else 0) for s in ma}
        for denom,vals in [('all_reads',va),('annotated_reads',vi)]:
            obs,p,d=exact_paired(vals)
            means={st:sum(vals[s] for s in ma if s.startswith(st))/9 for st in ('M','A')}
            loo=[]
            for gh in ('2','3','4'):
                keep=[x for i,x in enumerate(d) if str((i//3)+2)!=gh]; loo.append(sum(keep)/len(keep))
            stats.append({'confidence_set':confidence,'ontology':ontology,'denominator':denom,'function':cat,'mean_M':means['M'],'mean_A':means['A'],'mean_delta':obs,'pairs_positive':sum(x>0 for x in d),'exact_p':p,'loo_GH2':loo[0],'loo_GH3':loo[1],'loo_GH4':loo[2]})
            for s in ma: long.append({'confidence_set':confidence,'ontology':ontology,'denominator':denom,'sample_id':s,'stage':s[0],'greenhouse_id':'GH'+s[1],'pair_id':s[1:],'function':cat,'relative_abundance':vals[s]})
    # Control FDR separately for each denominator, as documented.  Pooling
    # all_reads and annotated_reads would duplicate each biological function
    # in the same correction family and make q values unnecessarily conservative.
    for denom in ('all_reads','annotated_reads'):
        subset=[x for x in stats if x['denominator']==denom]
        qs=bh([x['exact_p'] for x in subset])
        for x,q in zip(subset,qs): x['BH_q']=q
    all_outputs.extend(stats)
    prefix=f'{confidence}_{ontology}'
    with open(os.path.join(outdir,prefix+'_sample_abundance.tsv'),'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=long[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(long)
with open(os.path.join(outdir,'FUNGuild_paired_function_statistics.tsv'),'w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=all_outputs[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(all_outputs)
with open(os.path.join(outdir,'FUNGuild_analysis_notes.txt'),'w',encoding='utf-8') as f:
    f.write('Counts are fractionally allocated when an ASV has multiple trophic modes or guilds.\n')
    f.write('strict includes Probable and Highly Probable; inclusive additionally includes Possible.\n')
    f.write('Exact P values use all 2^9 paired sign flips; BH correction is within each confidence/ontology/denominator table.\n')
    f.write('Results are predicted ecological guilds, not measured functions or activities.\n')

