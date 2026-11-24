#!/usr/bin/env python3
"""Comprehensive REINVENT4 post-generation analysis. Requires rdkit, numpy, pandas, matplotlib, scikit-learn."""
import argparse, html, json, math, warnings
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.decomposition import PCA
RDLogger.DisableLog('rdApp.warning')

PROPS = ['MW','cLogP','TPSA','HBD','HBA','QED','SA','SCScore','RotB','FormalCharge','FractionCSP3','LipinskiViolations','ParentSimilarity','NLL']
NUMERIC = PROPS + ['RAscore']
PALETTE = plt.get_cmap('tab10')

def savefig(fig, base, dpi=300):
    fig.tight_layout()
    fig.savefig(str(base)+'.png', dpi=dpi, bbox_inches='tight')
    fig.savefig(str(base)+'.pdf', bbox_inches='tight')
    plt.close(fig)

def boolean_series(s):
    return s.astype('string').str.strip().str.lower().map({'true':1.,'false':0.,'1':1.,'0':0.,'yes':1.,'no':0.})

def plot_box(df, prop, directory):
    groups = [(k, g[prop].dropna().to_numpy(dtype=float)) for k,g in df.groupby('Hit',sort=True)]
    groups = [(k,v) for k,v in groups if len(v)]
    if not groups: return
    fig,ax=plt.subplots(figsize=(10,5.5))
    ax.boxplot([v for _,v in groups], tick_labels=[k.replace('Hit_','Hit ') for k,_ in groups],showfliers=False)
    ax.set(xlabel='Parent hit',ylabel=prop,title=f'{prop} by parent hit')
    ax.grid(axis='y',alpha=.2)
    savefig(fig,directory/f'{prop}_boxplot')

def plot_hist(df, prop, directory):
    if not df[prop].notna().any(): return
    fig,ax=plt.subplots(figsize=(9,5.5))
    for i,(k,g) in enumerate(df.groupby('Hit',sort=True)):
        v=g[prop].dropna().to_numpy(dtype=float)
        if len(v): ax.hist(v,bins=35,density=True,histtype='step',linewidth=1.5,label=k.replace('Hit_','Hit '),color=PALETTE(i))
    ax.set(xlabel=prop,ylabel='Density',title=f'{prop} distributions')
    ax.legend(ncol=2,fontsize=8); ax.grid(alpha=.15)
    savefig(fig,directory/f'{prop}_histogram')

def plot_bar(series, ylabel, title, base, percent=False):
    series=series.dropna()
    if series.empty:return
    fig,ax=plt.subplots(figsize=(9,5))
    bars=ax.bar([x.replace('Hit_','Hit ') for x in series.index],series.values,color=[PALETTE(i) for i in range(len(series))])
    ax.set(ylabel=ylabel,title=title)
    if percent:ax.set_ylim(0,105)
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2,b.get_height(),f'{b.get_height():.1f}' if percent else f'{b.get_height():,.0f}',ha='center',va='bottom',fontsize=8)
    ax.grid(axis='y',alpha=.2);savefig(fig,base)

def scatter(df,x,y,base):
    sub=df[['Hit',x,y]].dropna()
    if sub.empty:return
    fig,ax=plt.subplots(figsize=(8,6))
    for i,(hit,g) in enumerate(sub.groupby('Hit')):
        ax.scatter(g[x],g[y],s=9,alpha=.35,label=hit.replace('Hit_','Hit '),rasterized=True,color=PALETTE(i))
    ax.set(xlabel=x,ylabel=y,title=f'{y} vs {x}')
    ax.legend(markerscale=2,ncol=2,fontsize=8)
    ax.grid(alpha=.15);savefig(fig,base)

def heatmap(matrix, labels, title, base, vmin=None, vmax=None):
    fig,ax=plt.subplots(figsize=(max(7,len(labels)*.75),max(6,len(labels)*.65)))
    im=ax.imshow(matrix,cmap='viridis',vmin=vmin,vmax=vmax,aspect='auto')
    ax.set_xticks(range(len(labels)),labels,rotation=45,ha='right')
    ax.set_yticks(range(len(labels)),labels)
    fig.colorbar(im,ax=ax,shrink=.8)
    ax.set_title(title);savefig(fig,base)

def pareto_mask(values):
    # All columns are objectives to maximize. O(N^2) but only applied to shortlisted candidates.
    n=len(values);keep=np.ones(n,dtype=bool)
    for i in range(n):
        if not keep[i]:continue
        dominated=np.all(values>=values[i],axis=1)&np.any(values>values[i],axis=1)
        if dominated.any():keep[i]=False
    return keep

def run(args):
    root=Path(args.input).resolve();out=Path(args.output).resolve()
    if not root.is_dir():raise SystemExit(f'Input directory not found: {root}')
    if out==root:raise SystemExit('Output must differ from input')
    figs=out/'figures';tables=out/'tables';figs.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=args.fp_size)
    hits=sorted([p for p in root.glob('Hit_*') if p.is_dir() and (p/'generated.csv').exists()])
    if not hits:raise SystemExit(f'No Hit_*/generated.csv found under {root}')
    frames=[];fps_by_hit={};valid_counts={};invalid_rows=[];parent_map={};scaffold_match_info=[]
    for p in hits:
        d=pd.read_csv(p/'generated.csv',low_memory=False)
        if 'SMILES' not in d:raise SystemExit(f'Missing SMILES in {p}')
        d['Hit']=p.name;valid=[];mols=[]
        for idx,s in enumerate(d.SMILES):
            m=Chem.MolFromSmiles(s) if isinstance(s,str) and s.strip() else None
            if m is None:invalid_rows.append({'Hit':p.name,'row':idx+2,'SMILES':str(s)})
            else:valid.append(idx);mols.append(m)
        d=d.iloc[valid].copy().reset_index(drop=True)
        if not len(d):continue
        for prop in NUMERIC:
            if prop not in d:d[prop]=np.nan
            d[prop]=pd.to_numeric(d[prop],errors='coerce')
        for flag in ['ScaffoldRetained','PAINS','BRENK','RequiredFeaturesRetained','ForbiddenFeaturesPresent','RouteFound']:
            if flag in d:d[flag+'_numeric']=boolean_series(d[flag])
        d['CanonicalSMILES']=[Chem.MolToSmiles(m,isomericSmiles=True) for m in mols]
        d['MurckoScaffold']=[MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols]
        d['MurckoGenericScaffold']=[MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False) for m in mols]
        # Independently check literal scaffold SMILES when present, without claiming equivalence to pipeline's mapped core.
        exact=[];stereo=[]
        for mol,raw in zip(mols,d['scaffold'] if 'scaffold' in d else [None]*len(d)):
            sm=Chem.MolFromSmiles(raw) if isinstance(raw,str) and raw.strip() else None
            exact.append(float(mol.HasSubstructMatch(sm)) if sm else np.nan)
            stereo.append(float(mol.HasSubstructMatch(sm,useChirality=True)) if sm else np.nan)
        d['IndependentScaffoldSubstructureMatch']=exact
        d['IndependentScaffoldStereoMatch']=stereo
        fps=[gen.GetFingerprint(m) for m in mols]
        fps_by_hit[p.name]=fps;valid_counts[p.name]=len(d)
        frames.append(d)
        print(f'{p.name}: {len(d):,} valid molecules, {d.MurckoScaffold.nunique():,} Murcko scaffolds',flush=True)
    if not frames:raise SystemExit('No valid molecules')
    df=pd.concat(frames,ignore_index=True)
    df.to_csv(tables/'all_molecules_annotated.csv',index=False)
    pd.DataFrame(invalid_rows,columns=['Hit','row','SMILES']).to_csv(tables/'invalid_smiles.csv',index=False)
    counts=df.groupby('Hit').size().rename('Molecules')
    summary=pd.DataFrame(index=counts.index);summary['Molecules']=counts
    summary['UniqueMolecules']=df.groupby('Hit').CanonicalSMILES.nunique()
    summary['UniqueMurckoScaffolds']=df.groupby('Hit').MurckoScaffold.nunique()
    summary['ScaffoldToMoleculeRatio']=summary.UniqueMurckoScaffolds/summary.Molecules
    for col in PROPS:
        summary['Mean_'+col]=df.groupby('Hit')[col].mean()
        summary['Median_'+col]=df.groupby('Hit')[col].median()
    for flag in ['ScaffoldRetained','PAINS','BRENK','RequiredFeaturesRetained','ForbiddenFeaturesPresent','RouteFound']:
        name=flag+'_numeric'
        if name in df:summary[flag+'_Percent']=df.groupby('Hit')[name].mean()*100
    for flag in ['IndependentScaffoldSubstructureMatch','IndependentScaffoldStereoMatch']:
        summary[flag+'_Percent']=df.groupby('Hit')[flag].mean()*100
    # Distribution of exact pairwise similarities within each hit. Optional cap controls cost and memory.
    rng=np.random.default_rng(args.seed);diversity=[];hist_rows=[]
    for hit,fps in fps_by_hit.items():
        n=len(fps);total=n*(n-1)//2
        if n<2:continue
        # Exact all-pairs for <= max_exact; otherwise random pairs with replacement.
        if n<=args.max_exact:
            chunks=[np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i]),dtype=float) for i in range(1,n)]
            sims=np.concatenate(chunks)
            method='exact'
        else:
            count=min(args.pair_samples,total)
            ii=rng.integers(0,n,size=count);jj=rng.integers(0,n-1,size=count);jj+=jj>=ii
            sims=np.array([DataStructs.TanimotoSimilarity(fps[a],fps[b]) for a,b in zip(ii,jj)])
            method='sampled'
        diversity.append({'Hit':hit,'PairsTotal':total,'PairsEvaluated':len(sims),'Method':method,'MeanSimilarity':sims.mean(),'MeanDissimilarity':1-sims.mean(),'MedianSimilarity':np.median(sims),'P90Similarity':np.quantile(sims,.9)})
        hist,edges=np.histogram(sims,bins=np.linspace(0,1,51))
        hist_rows.extend({'Hit':hit,'BinStart':edges[j],'BinEnd':edges[j+1],'Count':int(h)} for j,h in enumerate(hist))
        fig,ax=plt.subplots(figsize=(8,5));ax.hist(sims,bins=np.linspace(0,1,51),color=PALETTE(list(fps_by_hit).index(hit)))
        ax.set(xlabel='Pairwise Morgan Tanimoto similarity',ylabel='Molecular pairs',title=f'{hit} within-library similarity ({method})');savefig(fig,figs/f'{hit}_pairwise_similarity')
    diversity_df=pd.DataFrame(diversity).set_index('Hit');summary=summary.join(diversity_df)
    diversity_df.to_csv(tables/'pairwise_diversity.csv',float_format='%.6f')
    pd.DataFrame(hist_rows).to_csv(tables/'pairwise_similarity_histograms.csv',index=False)
    # Cross-hit duplicates and scaffold overlaps.
    grouped=df.groupby('CanonicalSMILES').Hit.nunique()
    shared=set(grouped[grouped>1].index)
    df[df.CanonicalSMILES.isin(shared)].to_csv(tables/'cross_hit_duplicate_molecules.csv',index=False)
    labels=list(summary.index);scaffold_sets={h:set(df.loc[df.Hit==h,'MurckoScaffold']) for h in labels}
    overlap=np.zeros((len(labels),len(labels)),dtype=int)
    jaccard=np.zeros_like(overlap,dtype=float)
    for i,a in enumerate(labels):
        for j,b in enumerate(labels):
            intersect=len(scaffold_sets[a]&scaffold_sets[b]);union=len(scaffold_sets[a]|scaffold_sets[b]);overlap[i,j]=intersect;jaccard[i,j]=intersect/union if union else np.nan
    pd.DataFrame(overlap,index=labels,columns=labels).to_csv(tables/'murcko_scaffold_overlap_counts.csv')
    pd.DataFrame(jaccard,index=labels,columns=labels).to_csv(tables/'murcko_scaffold_jaccard.csv')
    heatmap(jaccard,[x.replace('Hit_','Hit ') for x in labels],'Murcko scaffold Jaccard overlap',figs/'murcko_scaffold_overlap',vmin=0,vmax=1)
    # Top scaffolds; empty Murcko scaffold denotes acyclic molecules.
    top=df.groupby(['Hit','MurckoScaffold']).size().rename('Count').reset_index().sort_values(['Hit','Count'],ascending=[True,False])
    top.to_csv(tables/'scaffold_frequencies.csv',index=False)
    top.groupby('Hit').head(15).to_csv(tables/'top15_scaffolds_per_hit.csv',index=False)
    # Statistical summaries, threshold compliance.
    stat=df.groupby('Hit')[PROPS].agg(['count','mean','std','median','min','max'])
    stat.to_csv(tables/'property_statistics.csv')
    rules={'QED_ge_0.5':df.QED>=.5,'QED_ge_0.7':df.QED>=.7,'SA_le_4':df.SA<=4,'SCScore_le_3.5':df.SCScore<=3.5,'SCScore_le_5':df.SCScore<=5,'Lipinski_zero_violations':df.LipinskiViolations==0}
    pass_table=pd.DataFrame({k:v.groupby(df.Hit).mean()*100 for k,v in rules.items()})
    pass_table.to_csv(tables/'threshold_pass_rates_percent.csv')
    # Pareto front on QED (higher), SA and SCScore (lower) for a prescreened candidate pool.
    # Rank is explicitly heuristic, not a binding prediction or synthetic-route guarantee.
    def norm(v,lo,hi,reverse=False):
        x=np.clip((v-lo)/(hi-lo),0,1)
        return 1-x if reverse else x
    q=norm(df.QED.to_numpy(dtype=float),0,1)
    sa=norm(df.SA.to_numpy(dtype=float),1,10,True)
    sc=norm(df.SCScore.to_numpy(dtype=float),1,5,True)
    parent=norm(df.ParentSimilarity.to_numpy(dtype=float),0,1)
    weights=np.array([args.weight_qed,args.weight_sa,args.weight_sc,args.weight_parent],dtype=float)
    if np.any(weights<0) or not weights.sum():raise SystemExit('Ranking weights must be nonnegative and sum > 0')
    vals=np.column_stack([q,sa,sc,parent]);missing=~np.isfinite(vals).all(axis=1)
    df['RankingScore']=np.where(missing,np.nan,vals@weights/weights.sum())
    ranked=df.sort_values('RankingScore',ascending=False,na_position='last')
    ranked[['Hit','CanonicalSMILES','RankingScore','QED','SA','SCScore','ParentSimilarity','IndependentScaffoldSubstructureMatch']].to_csv(tables/'candidate_ranking_all.csv',index=False)
    shortlist=ranked.dropna(subset=['RankingScore']).groupby('Hit',group_keys=False).head(args.top_n).copy()
    shortlist.to_csv(tables/'top_candidates_per_hit.csv',index=False)
    if len(shortlist):
        pareto=pareto_mask(np.column_stack([shortlist.QED.to_numpy(),-shortlist.SA.to_numpy(),-shortlist.SCScore.to_numpy()]))
        shortlist.loc[pareto].to_csv(tables/'pareto_candidates_shortlist.csv',index=False)
    # Core figures.
    plot_bar(summary.Molecules,'Valid molecules','Generated molecules by hit',figs/'molecule_counts')
    plot_bar(summary.UniqueMurckoScaffolds,'Unique scaffolds','Bemis–Murcko scaffold counts',figs/'murcko_scaffold_counts')
    plot_bar(summary.ScaffoldToMoleculeRatio*100,'Scaffolds / molecules (%)','Scaffold diversity ratio',figs/'scaffold_diversity_ratio',True)
    plot_bar(summary.MeanDissimilarity,'Mean dissimilarity','Mean pairwise Morgan dissimilarity',figs/'mean_pairwise_dissimilarity')
    for flag in ['ScaffoldRetained_Percent','IndependentScaffoldSubstructureMatch_Percent','IndependentScaffoldStereoMatch_Percent','PAINS_Percent','BRENK_Percent','RouteFound_Percent']:
        if flag in summary:plot_bar(summary[flag],'%','Per-hit '+flag,figs/flag,True)
    for prop in PROPS:
        if df[prop].notna().any():
            plot_box(df,prop,figs)
            plot_hist(df,prop,figs)
    for x,y in [('ParentSimilarity','QED'),('ParentSimilarity','SCScore'),('QED','SA'),('QED','SCScore'),('MW','cLogP'),('TPSA','cLogP')]:
        if df[[x,y]].notna().any().all():scatter(df,x,y,figs/f'{y}_vs_{x}')
    corr=df[PROPS].corr(numeric_only=True)
    corr.to_csv(tables/'descriptor_correlations.csv')
    heatmap(corr.to_numpy(),list(corr.columns),'Descriptor Pearson correlations',figs/'descriptor_correlations',vmin=-1,vmax=1)
    # Descriptor PCA: standardize with global mean/std, median-impute missing; exclude completely empty columns.
    features=[x for x in ['MW','cLogP','TPSA','HBD','HBA','RotB','FractionCSP3','FormalCharge'] if df[x].notna().any()]
    if len(features)>=2 and len(df)>=3:
        X=df[features].copy();X=X.fillna(X.median());std=X.std(ddof=0).replace(0,1);X=(X-X.mean())/std
        pca=PCA(n_components=2,random_state=args.seed);coords=pca.fit_transform(X)
        pca_df=df[['Hit','CanonicalSMILES']].copy();pca_df['PC1']=coords[:,0];pca_df['PC2']=coords[:,1]
        pca_df.to_csv(tables/'descriptor_pca_coordinates.csv',index=False)
        fig,ax=plt.subplots(figsize=(9,7))
        for i,(hit,g) in enumerate(pca_df.groupby('Hit')):
            ax.scatter(g.PC1,g.PC2,s=6,alpha=.3,label=hit.replace('Hit_','Hit '),color=PALETTE(i),rasterized=True)
        ax.set(xlabel=f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',ylabel=f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)',title='Descriptor-based chemical space (PCA)')
        ax.legend(markerscale=3,ncol=2,fontsize=8);savefig(fig,figs/'chemical_space_descriptor_PCA')
        pd.DataFrame({'Feature':features,'PC1':pca.components_[0],'PC2':pca.components_[1]}).to_csv(tables/'pca_loadings.csv',index=False)
    summary.to_csv(tables/'hit_summary.csv',float_format='%.6f')
    # Compact report with links and embedded relative PNGs.
    fig_files=sorted(figs.glob('*.png'))
    relative=lambda p:p.relative_to(out).as_posix()
    table_links=''.join(f'<li><a href="{html.escape(relative(p))}">{html.escape(p.name)}</a></li>' for p in sorted(tables.glob('*.csv')))
    gallery=''.join(f'<figure><a href="{html.escape(relative(p))}"><img loading="lazy" src="{html.escape(relative(p))}" alt="{html.escape(p.stem)}"></a><figcaption>{html.escape(p.stem.replace("_"," "))}</figcaption></figure>' for p in fig_files)
    key_cols=[c for c in ['Molecules','UniqueMurckoScaffolds','Mean_QED','Mean_SA','Mean_SCScore','Mean_ParentSimilarity','MeanDissimilarity','ScaffoldRetained_Percent','IndependentScaffoldSubstructureMatch_Percent'] if c in summary]
    summary_html=summary[key_cols].round(3).to_html(classes='data',border=0)
    report=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>REINVENT4 Post-generation Analysis</title><style>body{{font-family:Arial,sans-serif;margin:2rem auto;max-width:1300px;padding:0 1rem;color:#182333}}h1,h2{{color:#183a59}}table{{border-collapse:collapse;display:block;overflow-x:auto}}td,th{{padding:7px 10px;border-bottom:1px solid #ddd;text-align:right}}th:first-child,td:first-child{{text-align:left}}.gallery{{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:18px}}figure{{margin:0;border:1px solid #ddd;padding:12px;border-radius:8px}}img{{width:100%;height:auto}}figcaption{{font-size:.85rem}}.note{{background:#f4f6f8;padding:12px}}</style></head><body><h1>REINVENT4 post-generation analysis</h1><p>Input: {html.escape(str(root))}</p><p><b>{len(df):,}</b> valid rows; <b>{df.CanonicalSMILES.nunique():,}</b> distinct canonical SMILES; <b>{len(shared):,}</b> structures shared between hits.</p><p class="note">Scaffold substructure checks are based on the CSV scaffold SMILES and may differ from the generator's mapped-scaffold definition. ParentSimilarity is taken from generator output, not recomputed. Pairwise similarity uses Morgan radius 2, {args.fp_size}-bit fingerprints. Ranking is heuristic, not a predicted binding affinity. Missing alert and route values are not negative results. PCA is descriptor-based, not fingerprint-based.</p><h2>Hit summary</h2>{summary_html}<h2>Data tables</h2><ul>{table_links}</ul><h2>Figures ({len(fig_files)})</h2><div class="gallery">{gallery}</div></body></html>'''
    (out/'report.html').write_text(report,encoding='utf-8')
    (out/'run_settings.json').write_text(json.dumps(vars(args),indent=2,default=str),encoding='utf-8')
    print(f'\nDONE: {len(df):,} valid molecules; {df.CanonicalSMILES.nunique():,} distinct canonical SMILES; {len(shared):,} cross-hit shared structures')
    print(f'Figures: {len(fig_files)} PNG + PDF pairs')
    print(f'HTML report: {out/"report.html"}')
    print(f'Summary CSV: {tables/"hit_summary.csv"}')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True,help='Directory containing Hit_*/generated.csv')
    ap.add_argument('--output',required=True,help='Analysis output directory')
    ap.add_argument('--fp-size',type=int,default=2048)
    ap.add_argument('--max-exact',type=int,default=2000,help='Maximum hit size for exact all-pairs diversity')
    ap.add_argument('--pair-samples',type=int,default=200000)
    ap.add_argument('--top-n',type=int,default=50)
    ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--weight-qed',type=float,default=.35)
    ap.add_argument('--weight-sa',type=float,default=.25)
    ap.add_argument('--weight-sc',type=float,default=.25)
    ap.add_argument('--weight-parent',type=float,default=.15)
    run(ap.parse_args())
