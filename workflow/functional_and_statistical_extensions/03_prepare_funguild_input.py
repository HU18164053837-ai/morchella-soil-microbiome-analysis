import csv,sys
inp,out=sys.argv[1:]
with open(inp,encoding='utf-8') as f, open(out,'w',encoding='utf-8',newline='') as g:
    r=csv.DictReader(f,delimiter='\t'); w=csv.writer(g,delimiter='\t'); w.writerow(['OTU','Taxonomy'])
    for x in r:
        oid=x.get('ASV_ID') or x.get(r.fieldnames[0]); parts=[]
        for k,p in [('Kingdom','k'),('Phylum','p'),('Class','c'),('Order','o'),('Family','f'),('Genus','g'),('Species','s')]:
            v=x.get(k,'') or ''; parts.append(f'{p}__{v}')
        w.writerow([oid,';'.join(parts)])

