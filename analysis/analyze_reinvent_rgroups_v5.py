#!/usr/bin/env python3
"""REINVENT4 V5: symmetry-aware R-groups, parent features, clustering, prioritization.
Read-only inputs. RDKit, pandas, numpy, matplotlib required.
Usage: python analysis/analyze_reinvent_rgroups_v5.py --input results/v5/REINVENT_V5_All_8_Hits_Relaxed --output results/v5/REINVENT_V5_RGroup_Analysis_V5
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import math
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdMolDescriptors, Descriptors
from rdkit.ML.Cluster import Butina
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PARENTS = [
'NCC(CN)NCC(O)COc1cccc2ccccc12',
'Cc1nn(-c2ccc(F)cc2)c2sc(C(=O)N3C[C@@H](N)[C@H](N)C3)cc12',
'NCCCNCCOc1ccc(-c2ccccc2)cc1',
'NCCCC[C@H](N)C(=O)Nc1ccc(-c2nc3ccccc3o2)cc1',
'NCCCC[C@H](N)C(=O)N1CCc2onc(-c3cccc4ccccc34)c2C1',
'Cc1cc(N2CCC(CN3CCN(C)CC3)CC2)nc2cc3c(cc12)CCC3',
'COc1cc2c(cc1OC)-c1cc(-c3ccccc3)c(C#N)c(=N)n1CC2',
'CCn1cc(C(=O)NCCN(C)C)c(=O)c2cc(F)c(N3CCN(C)CC3)cc21',
]
# SMARTS are functional-group *presence*, not conservation of atom mapping or pharmacophore geometry.
FEATURES = {
    'primary_amine':'[NX3;H2;!$(NC=O)]',
    'secondary_amine':'[NX3;H1;!$(NC=O)]',
    'tertiary_amine':'[NX3;H0;!$(NC=O);!$(N=[*])]',
    'amide':'[NX3][CX3](=[OX1])',
    'hydroxyl':'[OX2H]',
    'ether':'[OD2]([#6])[#6]',
    'nitrile':'[CX2]#N',
    'fluorine':'[F]',
    'bromine':'[Br]',
    'piperazine':'N1CCNCC1',
    'aromatic_nitrogen':'[n]',
    'carbonyl':'[CX3]=[OX1]',
}
PATTERNS={k:Chem.MolFromSmarts(v) for k,v in FEATURES.items()}

def mol(s):
    try:return Chem.MolFromSmiles(str(s))
    except Exception:return None

def canon(s):
    m=mol(s)
    return Chem.MolToSmiles(m,canonical=True,isomericSmiles=True) if m else ''

def fp(m):
    return AllChem.GetMorganGenerator(radius=2,fpSize=2048).GetFingerprint(m)

def assemble(core_smi,parts):
    """Join two single-dummy fragments to mapped core (map0/unmapped first, map1 second)."""
    core=mol(core_smi)
    frags=[mol(x) for x in parts]
    if core is None or any(x is None for x in frags):return None
    dummies=[a for a in core.GetAtoms() if a.GetAtomicNum()==0]
    if len(dummies)!=2 or sum(a.GetAtomMapNum()==1 for a in dummies)!=1:return None
    dummies.sort(key=lambda a:(a.GetAtomMapNum()==1,a.GetIdx()))
    if any(a.GetDegree()!=1 for a in dummies):return None
    fds=[]
    for f in frags:
        ds=[a for a in f.GetAtoms() if a.GetAtomicNum()==0]
        if len(ds)!=1 or ds[0].GetDegree()!=1:return None
        fds.append(ds[0])
    try:
        offsets=[core.GetNumAtoms(),core.GetNumAtoms()+frags[0].GetNumAtoms()]
        combined=Chem.CombineMols(Chem.CombineMols(core,frags[0]),frags[1])
        rw=Chem.RWMol(combined)
        remove=[]
        for i in range(2):
            rw.AddBond(dummies[i].GetNeighbors()[0].GetIdx(),
                       offsets[i]+fds[i].GetNeighbors()[0].GetIdx(),Chem.BondType.SINGLE)
            remove.extend([dummies[i].GetIdx(),offsets[i]+fds[i].GetIdx()])
        for j in sorted(remove,reverse=True):rw.RemoveAtom(j)
        m=rw.GetMol();Chem.SanitizeMol(m)
        return Chem.MolToSmiles(m,canonical=True,isomericSmiles=True)
    except Exception:return None

def features(m):
    return {k:int(m.HasSubstructMatch(p)) for k,p in PATTERNS.items()}

def gallery(records,path,n=24):
    ms=[];leg=[]
    for smi,label in records:
        m=mol(smi)
        if m is None:continue
        ms.append(m);leg.append(str(label)[:65])
        if len(ms)>=n:break
    if ms:
        Draw.MolsToGridImage(ms,molsPerRow=4,subImgSize=(330,260),legends=leg).save(str(path))

def plot_bar(counter,path,title,n=18):
    pairs=counter.most_common(n)
    if not pairs:return
    fig,ax=plt.subplots(figsize=(11,max(4,len(pairs)*0.4+1)))
    ax.barh(range(len(pairs)),[v for _,v in pairs][::-1])
    ax.set_yticks(range(len(pairs)),[s[:48] for s,_ in pairs][::-1],fontsize=8)
    ax.set_title(title);ax.set_xlabel('Count');fig.tight_layout()
    fig.savefig(path,dpi=200);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)

def heatmap(counter0,counter1,pairs,path,title,top=15):
    xs=[k for k,_ in counter0.most_common(top)]
    ys=[k for k,_ in counter1.most_common(top)]
    if not xs or not ys:return
    a=np.array([[pairs.get((x,y),0) for y in ys] for x in xs])
    fig,ax=plt.subplots(figsize=(max(8,len(ys)*0.8+3),max(5,len(xs)*0.5+2)))
    im=ax.imshow(a,cmap='Blues',aspect='auto')
    ax.set_xticks(range(len(ys)),[x[:25] for x in ys],rotation=65,ha='right',fontsize=8)
    ax.set_yticks(range(len(xs)),[x[:35] for x in xs],fontsize=8)
    ax.set_xlabel('Site 1');ax.set_ylabel('Site 0');ax.set_title(title)
    fig.colorbar(im,ax=ax,label='Molecule count');fig.tight_layout()
    fig.savefig(path,dpi=200);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)

def clusters(fragments,cutoff):
    """Butina clustering by Morgan Tanimoto distance. IDs local to one design and site."""
    smis=sorted(set(fragments))
    valid=[(s,mol(s)) for s in smis]
    valid=[(s,m) for s,m in valid if m is not None]
    if not valid:return {}
    fps=[fp(m) for _,m in valid]
    dist=[]
    for i in range(1,len(fps)):
        sims=DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dist.extend(1-x for x in sims)
    groups=Butina.ClusterData(dist,len(fps),cutoff,isDistData=True)
    mapping={}
    for ci,group in enumerate(groups,1):
        for idx in group:mapping[valid[idx][0]]=ci
    return mapping

def num(v,default):
    try:
        x=float(v)
        return x if math.isfinite(x) else default
    except (TypeError,ValueError):return default

def prioritize(g,keep):
    """Heuristic property rank + greedy fingerprint diversity, not binding prediction."""
    candidates=[]
    for row in g:
        m=mol(row['SMILES'])
        if m is None:continue
        qed=num(row.get('QED'),0)
        sa=num(row.get('SA'),6)
        sc=num(row.get('SCScore'),6)
        sim=num(row.get('ParentSimilarity'),0)
        # normalized, transparent heuristic; avoid treating it as biological activity
        score=0.4*qed+0.25*max(0,min(1,(5-sa)/4))+0.15*max(0,min(1,(5-sc)/4))+0.2*sim
        candidates.append((score,row,fp(m)))
    candidates.sort(key=lambda x:-x[0])
    pool=candidates[:min(len(candidates),max(keep*10,100))]
    chosen=[];chosen_fp=[]
    while pool and len(chosen)<keep:
        if not chosen:
            ix=0;distance=1
        else:
            options=[]
            for i,(score,row,fingerprint) in enumerate(pool):
                max_sim=max(DataStructs.BulkTanimotoSimilarity(fingerprint,chosen_fp))
                options.append((0.65*score+0.35*(1-max_sim),i,1-max_sim))
            _,ix,distance=max(options)
        score,row,fingerprint=pool.pop(ix)
        chosen_fp.append(fingerprint)
        chosen.append({**row,'HeuristicPropertyScore':round(score,4),
                       'MinDistanceToPreviouslySelected':round(distance,4),
                       'SelectionRank':len(chosen)+1})
    return chosen

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--parents',help='Optional CSV with columns Hit,SMILES; overrides built-in eight parents')
    ap.add_argument('--top-n',type=int,default=24)
    ap.add_argument('--cluster-distance',type=float,default=0.4)
    ap.add_argument('--select-per-hit',type=int,default=30)
    args=ap.parse_args()
    if not 0<args.cluster_distance<1:ap.error('--cluster-distance must be between 0 and 1')
    out=Path(args.output);tables=out/'tables';figs=out/'figures';gals=out/'galleries'
    for d in (tables,figs,gals):d.mkdir(parents=True,exist_ok=True)
    parent_map={f'Hit_{i:04d}':s for i,s in enumerate(PARENTS,1)}
    if args.parents:
        parent_df=pd.read_csv(args.parents)
        if not {'Hit','SMILES'}.issubset(parent_df):raise ValueError('Parents CSV needs Hit,SMILES columns')
        parent_map.update(dict(zip(parent_df.Hit.astype(str),parent_df.SMILES.astype(str))))
    records=[];summaries=[];freq_rows=[];pair_rows=[];cluster_rows=[];feature_rows=[];selected=[]
    files=sorted(Path(args.input).glob('Hit_*/generated.csv'))
    if not files:raise SystemExit('No Hit_*/generated.csv found')
    for file in files:
        hit=file.parent.name;df=pd.read_csv(file,keep_default_na=False)
        if not {'SMILES','mapped_scaffold','design_id','R-groups'}.issubset(df):
            raise ValueError(f'Missing required columns in {file}')
        pm=mol(parent_map.get(hit,''))
        parent_feat=features(pm) if pm else {}
        parent_fp=fp(pm) if pm else None
        hit_records=[]
        for design,g in df.groupby('design_id',sort=True):
            if g.mapped_scaffold.nunique()!=1:raise ValueError(f'{hit}/{design} inconsistent mapped scaffold')
            core=str(g.mapped_scaffold.iloc[0]);cnt0=Counter();cnt1=Counter();cntpairs=Counter();unordered=Counter()
            drows=[];amb=failed=invalid=0
            for i,row in g.iterrows():
                smi=canon(row['SMILES']);parts=str(row['R-groups']).split('|')
                if not smi or len(parts)!=2 or not all(canon(x) for x in parts):
                    invalid+=1;continue
                parts=[canon(x) for x in parts]
                direct=assemble(core,parts)==smi
                swapped=assemble(core,parts[::-1])==smi
                if direct and swapped:
                    status='symmetry_ambiguous';amb+=1
                    r0=r1='';unordered[tuple(sorted(parts))]+=1
                elif direct:
                    status='direct';r0,r1=parts
                elif swapped:
                    status='swapped';r0,r1=parts[::-1]
                else:
                    status='unverified';failed+=1;r0=r1=''
                if status in ('direct','swapped'):
                    cnt0[r0]+=1;cnt1[r1]+=1;cntpairs[(r0,r1)]+=1
                m=mol(smi);ft=features(m)
                retained=sum(ft[k] for k in parent_feat if parent_feat[k])
                present=sum(parent_feat.values())
                rec={'Hit':hit,'DesignID':design,'SourceRow':int(i)+2,'SMILES':smi,
                     'MappedScaffold':core,'FragmentA':parts[0],'FragmentB':parts[1],
                     'Site0':r0,'Site1':r1,'UnorderedPair':' | '.join(sorted(parts)) if status=='symmetry_ambiguous' else '',
                     'AssignmentStatus':status,'QED':row.get('QED',''),'SA':row.get('SA',''),
                     'SCScore':row.get('SCScore',''),'ParentSimilarity':row.get('ParentSimilarity',''),
                     'ParentFeaturePresenceRetained':retained,'ParentFeaturePresenceTotal':present,
                     'ParentFeaturePresenceFraction':round(retained/present,4) if present else '',
                     'IndependentParentTanimoto':round(DataStructs.TanimotoSimilarity(fp(m),parent_fp),4) if parent_fp else ''}
                drows.append(rec);records.append(rec);hit_records.append(rec)
                feature_rows.append({'Hit':hit,'DesignID':design,'SMILES':smi,**ft,
                                     **{f'Parent_{k}':v for k,v in parent_feat.items()}})
            n=len(g);assigned=sum(cntpairs.values())
            summaries.append({'Hit':hit,'DesignID':design,'Molecules':n,'Assigned':assigned,
                              'SymmetryAmbiguous':amb,'Unverified':failed,'Invalid':invalid,
                              'UniqueSite0':len(cnt0),'UniqueSite1':len(cnt1),
                              'UniqueAssignedPairs':len(cntpairs),'UniqueUnorderedAmbiguousPairs':len(unordered),
                              'MappedScaffold':core})
            for site,cnt in [('site_0',cnt0),('site_1',cnt1)]:
                mapping=clusters(cnt.keys(),args.cluster_distance)
                for smi,count in cnt.most_common():
                    freq_rows.append({'Hit':hit,'DesignID':design,'Site':site,'FragmentSMILES':smi,
                                      'Count':count,'PercentOfAssigned':round(100*count/assigned,3) if assigned else 0,
                                      'ClusterID':mapping.get(smi,'')})
                    cluster_rows.append({'Hit':hit,'DesignID':design,'Site':site,
                                         'FragmentSMILES':smi,'ClusterID':mapping.get(smi,''),'Count':count})
                plot_bar(cnt,figs/f'{hit}_{design}_{site}_frequency.png',f'{hit} {design}: {site}')
                gallery([(s,f'n={c} | cluster={mapping.get(s)}') for s,c in cnt.most_common()],
                        gals/f'{hit}_{design}_{site}_fragments.png',args.top_n)
            for (a,b),count in cntpairs.items():
                pair_rows.append({'Hit':hit,'DesignID':design,'Type':'ordered','Site0':a,'Site1':b,
                                  'UnorderedPair':'','Count':count})
            for (a,b),count in unordered.items():
                pair_rows.append({'Hit':hit,'DesignID':design,'Type':'symmetry_unordered','Site0':'','Site1':'',
                                  'UnorderedPair':a+' | '+b,'Count':count})
            heatmap(cnt0,cnt1,cntpairs,figs/f'{hit}_{design}_pair_heatmap.png',
                    f'{hit} {design}: site-resolved pairs')
            # Gallery of molecules from this design, including ambiguous molecules
            sorted_rows=sorted(drows,key=lambda r:(-num(r['QED'],0),num(r['SA'],99)))
            gallery([(r['SMILES'],f"QED={r['QED']} SA={r['SA']}") for r in sorted_rows],
                    gals/f'{hit}_{design}_top_molecules.png',args.top_n)
            print(f'{hit}/{design}: {n} molecules, assigned={assigned}, symmetry={amb}, unverified={failed}',flush=True)
        if pm:
            gallery([(Chem.MolToSmiles(pm),'Original parent')]+[(r['SMILES'],f"QED={r['QED']}") for r in
                    sorted(hit_records,key=lambda r:-num(r['QED'],0))],gals/f'{hit}_parent_vs_generated.png',args.top_n)
        selected.extend(prioritize(hit_records,args.select_per_hit))
    for name,rows in [('assignments',records),('design_summary',summaries),
                      ('rgroup_frequencies',freq_rows),('rgroup_pairs',pair_rows),
                      ('rgroup_clusters',cluster_rows),('functional_groups',feature_rows),
                      ('prioritized_compounds',selected)]:
        pd.DataFrame(rows).to_csv(tables/f'{name}.csv',index=False)
    for hit,g in pd.DataFrame(selected).groupby('Hit'):
        gallery([(r.SMILES,f'Rank {r.SelectionRank} | score {r.HeuristicPropertyScore}')
                 for r in g.itertuples()],gals/f'{hit}_prioritized.png',args.top_n)
    s=pd.DataFrame(summaries)
    total=int(s.Molecules.sum());amb=int(s.SymmetryAmbiguous.sum());assigned=int(s.Assigned.sum())
    (out/'README.md').write_text(f'''# REINVENT4 V5 analysis

Total: {total}; unambiguous: {assigned}; symmetry-ambiguous: {amb}; unverified: {int(s.Unverified.sum())}; invalid: {int(s.Invalid.sum())}.

All generated molecules remain in assignments.csv. Symmetry-equivalent pairs are counted as **unordered pairs** in rgroup_pairs.csv, not site-specific R1/R2 counts. Site labels are **local to a design**. The assignment reflects structural reconstruction and does not prove which site the sampler used first.

Functional group retention is **presence/absence of SMARTS motifs** in the parent and generated molecule, not atom-level preservation of original groups or activity. Parent SMILES are built-in from the supplied eight hits, or can be overridden with --parents CSV (Hit,SMILES). Compare input normalization settings if parents were changed during generation.

Butina clustering uses Morgan radius 2, 2048-bit fingerprints on dummy-containing fragments, with distance cutoff {args.cluster_distance}; clusters are independent per design/site. This is chemical similarity clustering, not SAR activity clustering.

Prioritization: score=0.40*QED + 0.25*clamp((5-SA)/4) + 0.15*clamp((5-SCScore)/4) + 0.20*ParentSimilarity. Select from top max(100,10*N) candidates per hit by maximizing 0.65*score+0.35*(1-max Tanimoto to previously selected). No docking, affinity, selectivity, or synthesizability validation is implied. ParentSimilarity is read from generated.csv, not the independent fingerprint column.

Outputs: tables/assignments.csv, design_summary.csv, rgroup_frequencies.csv, rgroup_pairs.csv, rgroup_clusters.csv, functional_groups.csv, prioritized_compounds.csv; figures/ and galleries/.
''')
    print(f'COMPLETE: total={total} assigned={assigned} symmetry_ambiguous={amb}; output={out}')

if __name__=='__main__':main()
