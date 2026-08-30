import csv,json,sys
taxfile,dbfile,outfile=sys.argv[1:]
lines=open(dbfile,encoding='utf-8',errors='replace').read().splitlines()
payload=next(x.strip() for x in lines if x.lstrip().startswith('[{"taxon"'))
db=json.JSONDecoder().raw_decode(payload)[0]
lookup={str(x.get('taxon','')).replace(' ','_').lower():x for x in db}
ranks=['Species','Genus','Family','Order','Class','Phylum','Kingdom']
with open(taxfile,encoding='utf-8') as f, open(outfile,'w',encoding='utf-8',newline='') as g:
 r=csv.DictReader(f,delimiter='\t'); fields=r.fieldnames+['matched_taxon','matched_rank','trophic_mode','guild','growth_form','trait','confidence_ranking','guild_notes','citation_source']; w=csv.DictWriter(g,fields,delimiter='\t',extrasaction='ignore'); w.writeheader()
 for row in r:
  hit=None; rank=''
  for k in ranks:
   q=(row.get(k) or '').strip().replace(' ','_').lower()
   if q and q in lookup: hit=lookup[q]; rank=k; break
  row.update({'matched_taxon':hit.get('taxon','') if hit else '', 'matched_rank':rank,
   'trophic_mode':(hit.get('trophicMode') or hit.get('TrophicMode') or '') if hit else '',
   'guild':hit.get('guild','') if hit else '', 'growth_form':(hit.get('growthForm') or hit.get('growthMorphology') or '') if hit else '',
   'trait':hit.get('trait','') if hit else '', 'confidence_ranking':hit.get('confidenceRanking','') if hit else '',
   'guild_notes':hit.get('notes','') if hit else '', 'citation_source':hit.get('citationSource','') if hit else ''})
  w.writerow(row)

