#!/usr/bin/env python3
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 4:
    raise SystemExit("Usage: 15_parse_unite_blast_v2.py <hits.tsv> <query.fasta> <outdir>")

hits_path, query_path, outdir = map(Path, sys.argv[1:])
outdir.mkdir(parents=True, exist_ok=True)
ranks = [("Kingdom", "k__"), ("Phylum", "p__"), ("Class", "c__"),
         ("Order", "o__"), ("Family", "f__"), ("Genus", "g__"),
         ("Species", "s__")]


def parse_taxonomy(title):
    tax = {rank: "" for rank, _ in ranks}
    for rank, prefix in ranks:
        # UNITE titles introduce the taxonomy with a pipe, then delimit ranks
        # with semicolons: ...|k__Fungi;p__... . Accept either delimiter.
        match = re.search(r"(?:^|[;|])" + re.escape(prefix) + r"([^;|]+)", title)
        if match:
            value = match.group(1).strip()
            if value.lower() not in {"unidentified", "uncultured", ""}:
                tax[rank] = value
    sh = re.search(r"SH\d+\.\d+FU", title)
    tax["SH"] = sh.group(0) if sh else ""
    return tax


hits = defaultdict(list)
with hits_path.open(newline="", encoding="utf-8") as handle:
    for row in csv.reader(handle, delimiter="\t"):
        if len(row) != 15:
            raise ValueError(f"Unexpected BLAST column count: {len(row)}")
        hit = {
            "pident": float(row[2]), "length": int(row[3]),
            "evalue": float(row[10]), "bitscore": float(row[11]),
            "qlen": int(row[12]), "slen": int(row[13]),
            "qcov": float(row[14]), "title": row[1],
        }
        hit["taxonomy"] = parse_taxonomy(row[1])
        hits[row[0]].append(hit)


def fasta_records(path):
    name, seq = None, []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)


queries = list(fasta_records(query_path))
rows, unclassified = [], []
for qid, seq in queries:
    qhits = sorted(hits.get(qid, []), key=lambda x: (-x["bitscore"], -x["pident"]))
    if not qhits:
        rows.append({"ASV_ID": qid, **{r: "" for r, _ in ranks}, "SH": "",
                     "best_identity": "", "best_query_coverage": "", "best_evalue": "",
                     "best_bitscore": "", "top_hit_count": 0, "confidence": "unclassified"})
        unclassified.append((qid, seq))
        continue

    best = qhits[0]
    top = [h for h in qhits if h["bitscore"] >= best["bitscore"] * 0.99
           and h["pident"] >= best["pident"] - 0.5]
    consensus, blocked = {}, False
    for rank, _ in ranks:
        vals = {h["taxonomy"][rank] for h in top if h["taxonomy"][rank]}
        if blocked or len(vals) != 1:
            consensus[rank], blocked = "", True
        else:
            consensus[rank] = next(iter(vals))

    sh_vals = {h["taxonomy"]["SH"] for h in top if h["taxonomy"]["SH"]}
    sh = next(iter(sh_vals)) if len(sh_vals) == 1 else ""
    if best["pident"] >= 99 and best["qcov"] >= 90 and consensus["Species"] and sh:
        confidence = "high_species_SH_candidate"
    elif best["pident"] >= 97 and best["qcov"] >= 80 and consensus["Genus"]:
        confidence = "genus_candidate"
    elif consensus["Kingdom"]:
        confidence = "higher_rank_only"
    else:
        confidence = "ambiguous"
    rows.append({"ASV_ID": qid, **consensus, "SH": sh,
                 "best_identity": best["pident"], "best_query_coverage": best["qcov"],
                 "best_evalue": best["evalue"], "best_bitscore": best["bitscore"],
                 "top_hit_count": len(top), "confidence": confidence})

fields = ["ASV_ID"] + [r for r, _ in ranks] + ["SH", "best_identity",
          "best_query_coverage", "best_evalue", "best_bitscore", "top_hit_count", "confidence"]
with (outdir / "its_asv_taxonomy.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

with (outdir / "unclassified_asvs.fasta").open("w", encoding="ascii") as handle:
    for qid, seq in unclassified:
        handle.write(f">{qid}\n{seq}\n")

with (outdir / "taxonomy_annotation_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "count", "fraction"])
    total = len(rows)
    for rank, _ in ranks:
        count = sum(bool(r[rank]) for r in rows)
        writer.writerow([f"annotated_{rank.lower()}", count, count / total])
    for label in ("high_species_SH_candidate", "genus_candidate", "higher_rank_only",
                  "ambiguous", "unclassified"):
        count = sum(r["confidence"] == label for r in rows)
        writer.writerow([f"confidence_{label}", count, count / total])

print(f"PARSE_COMPLETE queries={len(rows)} unclassified={len(unclassified)}")

