#!/usr/bin/env python3
"""Site-normalized LibInvent R-group analysis, preserving ambiguous assignments.

Requires: rdkit, pandas, matplotlib, numpy, and analyze_reinvent_rgroups_v3.py
in the same scripts directory. Never modifies generated.csv.
"""
import argparse
from collections import Counter
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from analyze_reinvent_rgroups_v3 import assemble, canonical

def grid_gallery(items, path, n=24):
    mols, legends = [], []
    for smi, legend in items:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        mols.append(m)
        legends.append(str(legend)[:55])
        if len(mols) >= n:
            break
    if mols:
        Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(340, 270),
                             legends=legends).save(str(path))

def heatmap(matrix, rows, cols, path, title, xlabel, ylabel):
    if not rows or not cols:
        return
    fig, ax = plt.subplots(figsize=(max(8, 0.9*len(cols)+3),
                                    max(5, 0.45*len(rows)+2)))
    im = ax.imshow(matrix, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(rows)))
    ax.set_xticklabels([str(x)[:24] for x in cols], rotation=65, ha="right", fontsize=8)
    ax.set_yticklabels([str(x)[:34] for x in rows], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Molecule count")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)

def barplot(counter, path, title, limit=20):
    top=counter.most_common(limit)
    if not top:
        return
    fig,ax=plt.subplots(figsize=(11,max(4,0.36*len(top)+1)))
    ax.barh(range(len(top)),[n for _,n in top][::-1])
    ax.set_yticks(range(len(top)),[s[:48] for s,_ in top][::-1],fontsize=8)
    ax.set_xlabel("Number of molecules")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path,dpi=200)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",required=True)
    parser.add_argument("--output",required=True)
    parser.add_argument("--top-n",type=int,default=24)
    parser.add_argument("--heatmap-top",type=int,default=15)
    args=parser.parse_args()
    root=Path(args.input)
    out=Path(args.output)
    tables,figs,galleries=(out/"tables",out/"figures",out/"galleries")
    for d in (tables,figs,galleries):
        d.mkdir(parents=True,exist_ok=True)
    paths=sorted(root.glob("Hit_*/generated.csv"))
    if not paths:
        raise SystemExit(f"No Hit_*/generated.csv found under {root}")
    assignments=[]
    summaries=[]
    frequencies=[]
    combinations=[]
    for path in paths:
        hit=path.parent.name
        df=pd.read_csv(path,keep_default_na=False)
        required={"SMILES","design_id","mapped_scaffold","R-groups"}
        if not required.issubset(df.columns):
            raise ValueError(f"{path} missing {required-set(df.columns)}")
        for design,g in df.groupby("design_id",sort=True):
            cores=g["mapped_scaffold"].unique()
            if len(cores)!=1:
                raise ValueError(f"{hit}/{design}: inconsistent mapped scaffolds")
            core=str(cores[0])
            site0,site1=Counter(),Counter()
            pair_counts=Counter()
            pair_example={}
            ambiguous_pairs=Counter()
            n_direct=n_swapped=n_both=n_failed=n_invalid=0
            unique_molecules=set()
            for index,row in g.iterrows():
                parts=str(row["R-groups"]).split("|")
                smi=canonical(row["SMILES"])
                if len(parts)!=2 or not smi:
                    n_invalid+=1
                    continue
                frags=[canonical(x) for x in parts]
                if not all(frags):
                    n_invalid+=1
                    continue
                direct=(assemble(core,parts)==smi)
                swapped=(assemble(core,parts[::-1])==smi)
                if direct:
                    n_direct+=1
                if swapped:
                    n_swapped+=1
                unique_molecules.add(smi)
                if direct and swapped:
                    status="ambiguous"
                    n_both+=1
                    r0=r1=""
                    ambiguous_pairs[tuple(sorted(frags))]+=1
                elif direct:
                    status="direct"
                    r0,r1=frags
                elif swapped:
                    status="swapped"
                    r0,r1=frags[::-1]
                else:
                    status="unverified"
                    n_failed+=1
                    r0=r1=""
                if status in ("direct","swapped"):
                    site0[r0]+=1
                    site1[r1]+=1
                    pair_counts[(r0,r1)]+=1
                    pair_example.setdefault((r0,r1),smi)
                assignments.append({
                    "Hit":hit,"DesignID":design,"SourceRow":int(index)+2,
                    "SMILES":smi,"MappedScaffold":core,
                    "OriginalFragment1":frags[0],"OriginalFragment2":frags[1],
                    "Site0Fragment":r0,"Site1Fragment":r1,
                    "AssignmentStatus":status,
                    "DirectReconstruction":direct,"SwappedReconstruction":swapped,
                    "QED":row.get("QED",""),"SA":row.get("SA",""),
                    "SCScore":row.get("SCScore",""),
                    "ParentSimilarity":row.get("ParentSimilarity","")
                })
            n=len(g)
            assigned=sum(pair_counts.values())
            summary={
                "Hit":hit,"DesignID":design,"Molecules":n,
                "UniqueMolecules":len(unique_molecules),
                "DirectMatches":n_direct,"SwappedMatches":n_swapped,
                "UnambiguousAssignments":assigned,"AmbiguousAssignments":n_both,
                "Unverified":n_failed,"InvalidRows":n_invalid,
                "AssignedPercent":round(100*assigned/n,2) if n else 0,
                "UniqueSite0":len(site0),"UniqueSite1":len(site1),
                "UniqueAssignedPairs":len(pair_counts),
                "MappedScaffold":core,
            }
            summaries.append(summary)
            for site,counter in (("site_0",site0),("site_1",site1)):
                for fragment,count in counter.most_common():
                    frequencies.append({
                        "Hit":hit,"DesignID":design,"Site":site,
                        "FragmentSMILES":fragment,"Count":count,
                        "PercentOfAssigned":round(100*count/assigned,3) if assigned else 0
                    })
                barplot(counter,figs/f"{hit}_{design}_{site}_frequency.png",
                        f"{hit} {design} {site}: unambiguous assignments")
                grid_gallery([(s,f"{site} | n={v}") for s,v in counter.most_common()],
                             galleries/f"{hit}_{design}_{site}_fragments.png",args.top_n)
            for (r0,r1),count in pair_counts.most_common():
                combinations.append({
                    "Hit":hit,"DesignID":design,"Site0Fragment":r0,
                    "Site1Fragment":r1,"Count":count,
                    "ExampleSMILES":pair_example[(r0,r1)]
                })
            for (a,b),count in ambiguous_pairs.most_common():
                combinations.append({
                    "Hit":hit,"DesignID":design,"Site0Fragment":"",
                    "Site1Fragment":"","Count":count,
                    "ExampleSMILES":"",
                    "AmbiguousUnorderedPair":a+" | "+b
                })
            top0=[s for s,_ in site0.most_common(args.heatmap_top)]
            top1=[s for s,_ in site1.most_common(args.heatmap_top)]
            matrix=np.array([[pair_counts[(r0,r1)] for r1 in top1] for r0 in top0])
            heatmap(matrix,top0,top1,figs/f"{hit}_{design}_site0_site1_heatmap.png",
                    f"{hit} {design}: site-normalized R-group combinations",
                    "Site 1 fragment","Site 0 fragment")
            grid_gallery([(pair_example[k],f"Pair n={v}") for k,v in pair_counts.most_common()],
                         galleries/f"{hit}_{design}_common_assigned_pairs.png",args.top_n)
            print(f"{hit} {design}: n={n}, assigned={assigned}, ambiguous={n_both}, "
                  f"unverified={n_failed}, invalid={n_invalid}",flush=True)
    pd.DataFrame(assignments).to_csv(tables/"site_normalized_assignments.csv",index=False)
    pd.DataFrame(summaries).to_csv(tables/"design_summary.csv",index=False)
    pd.DataFrame(frequencies).to_csv(tables/"site_frequencies.csv",index=False)
    pd.DataFrame(combinations).to_csv(tables/"site_pair_frequencies.csv",index=False)
    total=sum(x["Molecules"] for x in summaries)
    assigned=sum(x["UnambiguousAssignments"] for x in summaries)
    ambiguous=sum(x["AmbiguousAssignments"] for x in summaries)
    failed=sum(x["Unverified"] for x in summaries)
    invalid=sum(x["InvalidRows"] for x in summaries)
    readme=f"""# REINVENT4 V4 site-normalized R-group analysis

Total molecules: {total}
Unambiguous site assignments: {assigned}
Ambiguous (both orientations reconstruct): {ambiguous}
Unverified reconstruction: {failed}
Invalid input rows: {invalid}

## Definitions
- Site 0 = the generation core attachment with dummy atom map number 0 (unmapped `*`).
- Site 1 = the generation core attachment with dummy atom map number 1 (`[*:1]`).
- Direct: original fragment order reconstructs the generated SMILES.
- Swapped: reversed fragment order reconstructs the generated SMILES.
- Ambiguous: both orders reconstruct the same SMILES. These rows are **excluded**
  from site-specific frequency and pair matrices, but retained in assignments.
- Unverified: neither orientation reconstructs the generated SMILES.
- PercentOfAssigned uses the number of unambiguous assignments as denominator.
- No binding affinity, biological SAR, or synthetic feasibility is inferred.
- Site labels are local to a particular design, not comparable across designs
  without an atom-mapping alignment.
- Fragments are canonicalized as single-dummy RDKit SMILES. Distinct strings
  can still represent chemically related groups.

## Tables
- tables/site_normalized_assignments.csv
- tables/design_summary.csv
- tables/site_frequencies.csv
- tables/site_pair_frequencies.csv

## Figures
- figures/*site0_site1_heatmap.png and PDF
- figures/*site_0_frequency.png and PDF
- figures/*site_1_frequency.png and PDF
- galleries/*_fragments.png and *_common_assigned_pairs.png
"""
    (out/"README.md").write_text(readme)
    print(f"TOTAL: {total}, unambiguous={assigned}, ambiguous={ambiguous}, "
          f"unverified={failed}, invalid={invalid}")
    print(f"Outputs: {out}")

if __name__=="__main__":
    main()
