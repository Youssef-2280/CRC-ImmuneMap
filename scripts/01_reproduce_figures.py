#!/usr/bin/env python3
"""Regenerate figures from approved primary_07e patient-level CSVs (no UCell rerun)."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

BASE=Path(__file__).resolve().parents[1]
DATA=BASE/'data'/'derived'
OUT=BASE/'figures'/'main'; OUT.mkdir(parents=True,exist_ok=True)
labels={'SPP1_TAM':'SPP1 TAM (macrophages)','FAP_CAF':'FAP CAF (fibroblasts)',
        'CD8_Dysfunction':'CD8 dysfunction','CD8_Cytotoxicity':'CD8 cytotoxicity'}
programs=['SPP1_TAM','FAP_CAF','CD8_Dysfunction','CD8_Cytotoxicity']
patient_colors={p:c for p,c in zip(['KUL01','KUL19','KUL21','KUL28','KUL30','KUL31'],
                 ['#bf5b4f','#c29a28','#24856b','#3297ad','#6684c5','#a46ba6'])}
BLUE='#345a79';INK='#263746';GRAY='#657782'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.labelcolor':INK,'text.color':INK,
                     'xtick.color':GRAY,'ytick.color':GRAY,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.edgecolor':'#c6cdd1',
                     'grid.alpha':.45,'savefig.dpi':210,'svg.fonttype':'none'})

def finish(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=220,facecolor='white',bbox_inches='tight')
    plt.close(fig)
    print('Saved',name)

D=pd.read_csv(DATA/'discovery_primary_associations.csv')
assert len(D)==3 and set(D.n_patients)=={11}
# Figure 1: exact discovery correlations, not CI forest plot.
comparison_text={
'SPP1_TAM__FAP_CAF':'SPP1 TAM vs FAP CAF',
'SPP1_TAM__CD8_Dysfunction':'SPP1 TAM vs CD8 dysfunction',
'FAP_CAF__CD8_Dysfunction':'FAP CAF vs CD8 dysfunction'}
D=D.set_index('comparison').loc[list(comparison_text)].reset_index()
fig,ax=plt.subplots(figsize=(10.4,4.5))
fig.patch.set_facecolor('white')
ys=np.arange(3)
ax.axvline(0,color='#97a5ac',linestyle=(0,(4,4)),lw=1.2,zorder=0)
for i,r in enumerate(D.itertuples(index=False)):
    ax.plot([0,r.rho],[i,i],lw=3,color=BLUE,alpha=.7)
    ax.scatter(r.rho,i,s=155,color=BLUE,edgecolor='white',lw=1.6,zorder=3)
    ax.text(.75,i,f'ρ = {r.rho:+.3f}   |   BH p = {r.p_BH_permutation:.3f}',
            ha='left',va='center',fontweight='medium',fontsize=10)
ax.set_yticks(ys,[comparison_text[x] for x in D.comparison]);ax.invert_yaxis()
ax.set_xlim(-.35,1.55);ax.set_xticks([-.25,0,.25,.5]);ax.set_xlabel('Patient-level Spearman ρ')
ax.set_title('Discovery: no statistically significant cross-program association',loc='left',fontsize=14,fontweight='bold',pad=24)
fig.text(.125,.91,'GSE132465  •  n = 11 complete-case tumor patients  •  permutation p-values, BH-adjusted',fontsize=9.5,color=GRAY)
ax.grid(axis='x',color='#e2e7ea');ax.set_axisbelow(True)
fig.text(.125,.045,'Exploratory comparison; dots are effect estimates, not confidence intervals.',fontsize=9,color=GRAY)
fig.subplots_adjust(left=.31,right=.97,top=.78,bottom=.2)
finish(fig,'01_discovery_associations')

# Pair scores calculate from locally sourced region aggregates, verify published summaries.
raw=pd.read_csv(DATA/'candidate_scores_by_patient_region.csv')
combined=pd.read_csv(DATA/'candidate_scores_core_plus_border_per_patient.csv')
N=raw[(raw.variant=='primary_07e')&(raw.tissue_region=='Normal')].copy()
T=combined[(combined.variant=='primary_07e')&(combined.tissue_region=='Tumor_Core_plus_Border')].copy()
paired=T.merge(N,on=['variant','patient_id','lineage','program'],suffixes=('_tumor','_normal'))
paired=paired[(paired.n_cells_tumor>=20)&(paired.n_cells_normal>=20)].copy()
paired['delta']=paired.mean_score_tumor-paired.mean_score_normal
pub=pd.read_csv(DATA/'validation_primary_paired_ge20.csv')
pub=pub[pub.tumor_region=='Tumor_Core_plus_Border'].set_index('program')
for key in programs:
    g=paired[paired.program==key];r=pub.loc[key]
    assert len(g)==r.n_matched_patients,(key,len(g))
    assert int((g.delta>0).sum())==r.n_tumor_higher
    assert abs(float(g.delta.median())-float(r.median_paired_delta))<1e-8

# Figure 2: matched slopes, colors per donor; correct denominators.
fig,axs=plt.subplots(2,2,figsize=(11.6,8.1));axs=axs.flatten()
for ax,pr in zip(axs,programs):
    g=paired[paired.program==pr].sort_values('patient_id')
    for r in g.itertuples(index=False):
        color=patient_colors[r.patient_id]
        ax.plot([0,1],[r.mean_score_normal,r.mean_score_tumor],lw=1.8,alpha=.78,color=color,zorder=2)
        ax.scatter([0,1],[r.mean_score_normal,r.mean_score_tumor],s=47,color=color,edgecolor='white',linewidth=.7,zorder=3)
    arr=pd.concat([g.mean_score_normal,g.mean_score_tumor]).to_numpy()
    rng=max(arr)-min(arr); pad=max(.015,.13*rng)
    ax.set_ylim(min(arr)-pad,max(arr)+pad)
    ax.set_xlim(-.22,1.22);ax.set_xticks([0,1],['Normal','Tumor Core + Border'])
    ax.set_title(f'{labels[pr]}   (n={len(g)})',fontsize=11.5,fontweight='bold',loc='left',pad=11)
    ax.text(.03,.92,f'Median matched Δ: {g.delta.median():+.3f}',transform=ax.transAxes,fontsize=9,color=GRAY)
    ax.grid(axis='y',color='#e3e9ec');ax.set_axisbelow(True)
axs[0].set_ylabel('Mean UCell score');axs[2].set_ylabel('Mean UCell score')
fig.suptitle('Validation: matched patient scores in tumor vs normal',x=.07,y=.975,ha='left',fontsize=15,fontweight='bold')
handles=[plt.Line2D([0],[0],color=color,marker='o',linestyle='-',lw=1.3,markersize=5,label=p) for p,color in patient_colors.items()]
fig.legend(handles=handles,loc='lower center',ncol=6,bbox_to_anchor=(.53,.052),frameon=False,fontsize=9)
fig.text(.5,.018,'≥20 eligible cells in each paired side, per program; patients eligible for each program differ. No p-values claimed.',ha='center',fontsize=8.8,color=GRAY)
fig.subplots_adjust(left=.10,right=.96,top=.87,bottom=.16,wspace=.22,hspace=.36)
finish(fig,'02_validation_scores')

# Figure 3: regional paired deltas from source patient/region rows, independent observations only across patients.
region_names=[('Tumor_Core','Core'),('Tumor_Border','Border'),('Tumor_Core_plus_Border','Core + Border')]
all_tumor=pd.concat([raw[raw.variant=='primary_07e'],combined[combined.variant=='primary_07e']],ignore_index=True)
all_tumor=all_tumor[all_tumor.tissue_region.isin([k for k,_ in region_names])]
per_region=all_tumor.merge(N,on=['variant','patient_id','lineage','program'],suffixes=('_tumor','_normal'))
per_region=per_region[(per_region.n_cells_tumor>=20)&(per_region.n_cells_normal>=20)].copy()
per_region['delta']=per_region.mean_score_tumor-per_region.mean_score_normal
for (region,pr),g in per_region.groupby(['tissue_region_tumor','program']):
    row=pub.loc[pr] if region=='Tumor_Core_plus_Border' else pd.read_csv(DATA/'validation_primary_paired_ge20.csv').set_index(['tumor_region','program']).loc[(region,pr)]
    assert len(g)==row.n_matched_patients
    assert abs(g.delta.median()-row.median_paired_delta)<1e-8
fig,axs=plt.subplots(2,2,figsize=(11.7,8));axs=axs.flatten()
for ax,pr in zip(axs,programs):
    g=per_region[per_region.program==pr]
    ax.axhline(0,color='#96a4ae',lw=1.25,linestyle=(0,(4,3)))
    for i,(region,short) in enumerate(region_names):
        q=g[g.tissue_region_tumor==region]
        for p in q.itertuples(index=False):
            ox=(list(patient_colors).index(p.patient_id)-2.5)*.056
            ax.scatter(i+ox,p.delta,s=56,color=patient_colors[p.patient_id],edgecolor='white',linewidth=.75,zorder=3)
        if not q.empty:
            md=q.delta.median();ax.hlines(md,i-.30,i+.30,lw=3.2,color=INK,zorder=4)
            ax.text(i,1.02,f'n={len(q)}',ha='center',va='bottom',transform=ax.get_xaxis_transform(),fontsize=9,color=GRAY)
    ax.set_xticks(range(3),[b for a,b in region_names]);ax.set_xlim(-.43,2.43)
    ax.set_title(labels[pr],loc='left',fontweight='bold',fontsize=11.5,pad=16)
    ax.grid(axis='y',color='#e4e9ec');ax.set_axisbelow(True)
    ax.set_ylabel('Tumor − Normal (mean UCell)')
fig.suptitle('Matched tumor-minus-normal differences across regions',x=.07,y=.983,ha='left',fontsize=14.8,fontweight='bold')
fig.text(.5,.025,'Each dot = one eligible patient; black segments = median. Core and Border are NOT independent cohorts.',ha='center',color=GRAY,fontsize=9)
fig.subplots_adjust(left=.10,right=.97,top=.86,bottom=.11,hspace=.39,wspace=.24)
finish(fig,'03_paired_differences')

# Figure 4: verified heatmap of median delta; same numeric scale within figure (not standardized effect).
S=pd.read_csv(DATA/'validation_selection_sensitivity_combined_ge20.csv')
S=S[(S.tumor_region=='Tumor_Core_plus_Border')&(S.min_cells_per_region==20)]
variants=['primary_07e','crosslineage_sensitivity','omit_high_risk_clusters','all_published_cd8_sensitivity']
varlabels=['Approved primary','Exclude cross-lineage flags','Omit high-risk clusters','All published CD8 sensitivity']
a=np.zeros((4,4)); counts=np.zeros((4,4),dtype=int); dirs=np.zeros((4,4),dtype=int)
for i,pr in enumerate(programs):
 for j,v in enumerate(variants):
    sub=S[(S.program==pr)&(S.variant==v)]
    assert len(sub)==1,(pr,v,len(sub))
    r=sub.iloc[0];a[i,j]=r.median_paired_delta;counts[i,j]=r.n_matched_patients;dirs[i,j]=r.n_tumor_higher
fig,ax=plt.subplots(figsize=(11.8,5.9))
norm=TwoSlopeNorm(vmin=-.10,vcenter=0,vmax=.40)
img=ax.imshow(a,cmap='PuOr',norm=norm,aspect='auto')
ax.set_yticks(range(4),['SPP1 TAM','FAP CAF','CD8 dysfunction','CD8 cytotoxicity'])
ax.set_xticks(range(4),varlabels,rotation=18,ha='right')
for i in range(4):
 for j in range(4):
    col='white' if abs(a[i,j])>.27 else INK
    ax.text(j,i,f'Δ {a[i,j]:+.3f}\n{dirs[i,j]}/{counts[i,j]} ↑',ha='center',va='center',fontsize=9.7,color=col,fontweight='medium')
for edge in np.arange(-.5,4,.99999):
 ax.axhline(edge,color='white',lw=3);ax.axvline(edge,color='white',lw=3)
cb=fig.colorbar(img,ax=ax,shrink=.8,pad=.03)
cb.set_label('Median paired tumor − normal UCell score (raw Δ)')
ax.set_title('Sensitivity: direction versus inclusion rules',loc='left',fontsize=15,fontweight='bold',pad=19)
fig.text(.13,.04,'Each cell: median Δ and number tumor-higher / number eligible. Program scores are NOT standardized cross-program effect sizes.',fontsize=9,color=GRAY)
fig.subplots_adjust(left=.20,right=.93,top=.84,bottom=.24)
finish(fig,'04_sensitivity_analysis')

# Figure 5: number of compartments >=20 cells by region and patient, 6 distinct patients.
C=pd.read_csv(DATA/'patient_three_program_coverage.csv')
C=C[(C.variant=='primary_07e')&(C.tissue_region.isin([r for r,_ in region_names]))].copy()
C['n_programs']=C[['SPP1_TAM','FAP_CAF','CD8_Dysfunction']].ge(20).sum(axis=1)
pats=sorted(C.patient_id.unique());regionIDs=[r for r,_ in region_names]
assert len(pats)==6
matrix=np.zeros((6,3),dtype=int)
for i,p in enumerate(pats):
 for j,r in enumerate(regionIDs):
    q=C[(C.patient_id==p)&(C.tissue_region==r)]
    assert len(q)==1
    matrix[i,j]=int(q.n_programs.iloc[0])
assert (matrix[:,0]==3).sum()==2 and (matrix[:,1]==3).sum()==3 and (matrix[:,2]==3).sum()==3
fig,ax=plt.subplots(figsize=(9.4,6.3))
cmap=ListedColormap(['#f4f4f4','#d8e5ec','#91b5c4','#397b91']);nm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4)
img=ax.imshow(matrix,cmap=cmap,norm=nm,aspect='auto')
ax.set_yticks(range(6),pats);ax.set_xticks(range(3),['Tumor Core\n(all three: 2/6)','Tumor Border\n(all three: 3/6)','Core + Border\n(all three: 3/6)'])
for i in range(6):
 for j in range(3):
    ax.text(j,i,f'{matrix[i,j]} / 3',fontsize=12.5,ha='center',va='center',color='white' if matrix[i,j]==3 else INK,fontweight='bold')
ax.set_title('Validation: three-program coverage is sparse',loc='left',fontsize=15,fontweight='bold',pad=20)
for edge in np.arange(-.5,3,.99999):ax.axvline(edge,color='white',lw=4)
for edge in np.arange(-.5,6,.99999):ax.axhline(edge,color='white',lw=4)
fig.text(.16,.06,'Number of programs with ≥20 candidate cells for each patient/region. Core + Border = ONE patient-level aggregate.',fontsize=9,color=GRAY)
fig.subplots_adjust(left=.17,right=.94,top=.79,bottom=.18)
finish(fig,'05_patient_coverage')

assert len(list(OUT.glob('*.png')))==5
assert len(list(OUT.glob('*.svg')))==0
print('SUCCESS: all 5 figures generated and original summary values verified.')
