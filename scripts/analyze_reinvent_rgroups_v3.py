#!/usr/bin/env python3
"""LibInvent design-aware R-group analysis (read-only)."""
import argparse
from collections import Counter
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def canonical(s):
    m=Chem.MolFromSmiles(str(s))
    return Chem.MolToSmiles(m, isomericSmiles=True) if m else None

def assemble(mapped, fragments):
    """Join single-dummy fragments to two mapped-core dummy atoms; return canonical product."""
    core=Chem.MolFromSmiles(mapped)
    if core is None: return None
    ds=[a for a in core.GetAtoms() if a.GetAtomicNum()==0]
    if len(ds)!=2 or any(a.GetDegree()!=1 for a in ds): return None
    # Atom map 1 marks the second generation site; atom map 0 is unmarked in RDKit.
    ds.sort(key=lambda a:(a.GetAtomMapNum()==1,a.GetIdx()))
    mols=[Chem.MolFromSmiles(x) for x in fragments]
    if any(m is None for m in mols): return None
    fd=[]
    for m in mols:
        d=[a for a in m.GetAtoms() if a.GetAtomicNum()==0]
        if len(d)!=1 or d[0].GetDegree()!=1:return None
        fd.append((d[0].GetIdx(),d[0].GetNeighbors()[0].GetIdx()))
    try:
        n=core.GetNumAtoms(); n1=mols[0].GetNumAtoms()
        rw=Chem.RWMol(Chem.CombineMols(Chem.CombineMols(core,mols[0]),mols[1]))
        remove=[]
        for j,d in enumerate(ds):
            offset=n+(n1 if j else 0)
            rw.AddBond(d.GetNeighbors()[0].GetIdx(),offset+fd[j][1],Chem.BondType.SINGLE)
            remove.extend([d.GetIdx(),offset+fd[j][0]])
        for i in sorted(remove,reverse=True):rw.RemoveAtom(i)
        product=rw.GetMol();Chem.SanitizeMol(product)
        return Chem.MolToSmiles(product,isomericSmiles=True)
    except Exception:return None

def gallery(records,path,n):
    mols=[];labels=[]
    for smi,label in records:
        m=Chem.MolFromSmiles(smi)
        if m:
            mols.append(m);labels.append(str(label)[:45])
        if len(mols)>=n:break
    if mols:
        Draw.MolsToGridImage(mols,molsPerRow=4,subImgSize=(330,250),legends=labels).save(str(path))

def plot_counts(counts,path,title):
    top=counts.most_common(15)
    if not top:return
    fig,ax=plt.subplots(figsize=(11,6))
    ax.barh(range(len(top)),[v for _,v in top][::-1])
    ax.set_yticks(range(len(top)),[(s[:50]+'…' if len(s)>50 else s) for s,_ in top][::-1],fontsize=8)
    ax.set_xlabel('Molecules');ax.set_title(title)
    fig.tight_layout();fig.savefig(path,dpi=200);plt.close(fig)

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--input',required=True);a.add_argument('--output',required=True)
    a.add_argument('--top-n',type=int,default=24)
    args=a.parse_args()
    root=Path(args.input);out=Path(args.output)
    tables=out/'tables';figs=out/'figures';grids=out/'galleries'
    for p in (tables,figs,grids):p.mkdir(parents=True,exist_ok=True)
    summaries=[];assignments=[];freq=[];pairfreq=[];warnings=[]
    files=sorted(root.glob('Hit_*/generated.csv'))
    if not files:raise SystemExit(f'No Hit_*/generated.csv found under {root}')
    for file in files:
        hit=file.parent.name
        df=pd.read_csv(file,keep_default_na=False)
        for col in ('SMILES','design_id','mapped_scaffold','R-groups'):
            if col not in df:raise ValueError(f'{file}: missing {col}')
        for design,g in df.groupby('design_id',sort=True):
            cores=g.mapped_scaffold.unique()
            if len(cores)!=1:
                warnings.append(dict(Hit=hit,Design=design,Issue='multiple_mapped_scaffolds'))
            core=str(cores[0]);r1=Counter();r2=Counter();pairs=Counter();examples={};ranked=[]
            direct=swapped=both=failed=invalid=0
            for idx,row in g.iterrows():
                parts=str(row['R-groups']).split('|')
                if len(parts)!=2:
                    invalid+=1;continue
                x,y=map(canonical,parts)
                if not x or not y:invalid+=1;continue
                smi=canonical(row.SMILES)
                if not smi:invalid+=1;continue
                r1[x]+=1;r2[y]+=1;pairs[(x,y)]+=1
                examples.setdefault((x,y),smi)
                fwd=assemble(core,parts)==smi
                rev=assemble(core,parts[::-1])==smi
                direct+=int(fwd);swapped+=int(rev);both+=int(fwd and rev)
                failed+=int(not fwd and not rev)
                assignments.append(dict(Hit=hit,DesignID=design,SourceRow=int(idx)+2,SMILES=smi,
                    R1=x,R2=y,DirectReconstruction=fwd,SwappedReconstruction=rev,
                    QED=row.get('QED',''),SA=row.get('SA',''),SCScore=row.get('SCScore',''),
                    ParentSimilarity=row.get('ParentSimilarity','')))
                def num(k,default):
                    try:return float(row[k])
                    except (ValueError,TypeError,KeyError):return default
                ranked.append((-num('QED',0),num('SA',99),num('SCScore',99),smi))
            n=sum(pairs.values())
            if not n:continue
            summaries.append(dict(Hit=hit,DesignID=design,Molecules=n,UniqueR1=len(r1),UniqueR2=len(r2),
                UniquePairs=len(pairs),UniquePairPercent=round(100*len(pairs)/n,2),
                DirectMatches=direct,SwappedMatches=swapped,BothMappingsMatch=both,
                ReconstructionUnverified=failed,InvalidRows=invalid,MappedScaffold=core))
            for site,counter in [('R1',r1),('R2',r2)]:
                for fragment,count in counter.most_common():
                    freq.append(dict(Hit=hit,DesignID=design,Site=site,FragmentSMILES=fragment,
                                     Count=count,Percent=round(100*count/n,3)))
                plot_counts(counter,figs/f'{hit}_{design}_{site}_frequency.png',f'{hit} {design} {site}')
                gallery([(s,f'{site} n={v}') for s,v in counter.most_common()],
                        grids/f'{hit}_{design}_{site}_fragments.png',args.top_n)
            for (x,y),count in pairs.most_common():
                pairfreq.append(dict(Hit=hit,DesignID=design,R1=x,R2=y,Count=count,
                                     ExampleSMILES=examples[(x,y)]))
            ranked.sort()
            gallery([(v[3],f'Rank {i+1}') for i,v in enumerate(ranked)],
                    grids/f'{hit}_{design}_top_molecules.png',args.top_n)
            gallery([(examples[p],f'Pair count={v}') for p,v in pairs.most_common()],
                    grids/f'{hit}_{design}_frequent_pairs.png',args.top_n)
            print(f'{hit} {design}: n={n}, unique R1={len(r1)}, R2={len(r2)}, pairs={len(pairs)}, '
                  f'direct={direct}, swapped={swapped}, unverified={failed}',flush=True)
    for name,data in [('design_summary',summaries),('rgroup_assignments',assignments),
                      ('rgroup_frequencies',freq),('rgroup_pair_frequencies',pairfreq),('warnings',warnings)]:
        pd.DataFrame(data).to_csv(tables/f'{name}.csv',index=False)
    print(f'Done: {out}')
    print('IMPORTANT: R1/R2 are sampler fragment order, not validated site labels.')
    print('Inspect DirectMatches/SwappedMatches before mapping fragments to scaffold positions.')

if __name__=='__main__':main()
