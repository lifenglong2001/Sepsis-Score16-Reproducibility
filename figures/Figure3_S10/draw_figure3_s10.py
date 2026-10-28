"""Source-driven B1000 figure update. Approved plotting functions, fixed layouts; no fitting."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
ROOT=Path(__file__).resolve().parent
S=ROOT/'source_data'
def read(p):return pd.read_csv(S/p,sep='\t' if str(p).endswith('.tsv') else ',')
DEEP_TEAL = '#3F827C'
DUSTY_BLUE = '#7298B6'
SAGE = '#A3B8A9'
TERRACOTTA = '#C77B61'
OCHRE = '#C8A45E'
DUSTY_ROSE = '#D39AA3'
WARM_SAND = '#D8CDBF'
WARM_TAUPE = '#BFAF9F'
INK = '#333333'
AXIS_GRAY = '#77716C'
LIGHT = '#E5DDD2'
BLUE, TEAL, PUR, AMBER, GRAY, ROSE = DUSTY_BLUE, DEEP_TEAL, DEEP_TEAL, OCHRE, WARM_TAUPE, DUSTY_ROSE
CTS_COLORS = [DUSTY_BLUE, SAGE, DUSTY_ROSE]
DIVERGING = LinearSegmentedColormap.from_list(
    'unified_blue_ivory_terracotta', ['#7F929D', '#F5F1E9', TERRACOTTA], N=256
)
mpl.rcParams.update({
    'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':7.5,
    'axes.titleweight':'normal','axes.linewidth':.7,'xtick.major.width':.7,
    'ytick.major.width':.7,'xtick.major.size':2.5,'ytick.major.size':2.5,
    'legend.fontsize':6.3,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'savefig.facecolor':'white','text.color':INK,'axes.labelcolor':INK,
    'xtick.color':INK,'ytick.color':INK
})


LABELS = {'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY': 'ROS pathway', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING': 'WNT/β-catenin', 'HALLMARK_IL2_STAT5_SIGNALING': 'IL-2/STAT5 signaling', 'HALLMARK_TGF_BETA_SIGNALING': 'TGFβ signaling', 'HALLMARK_MYC_TARGETS_V1': 'MYC targets V1', 'HALLMARK_MYC_TARGETS_V2': 'MYC targets V2', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING': 'PI3K/AKT/mTOR signaling', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNFα/NF-κB signaling', 'HALLMARK_IL6_JAK_STAT3_SIGNALING': 'IL-6/JAK/STAT3', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION': 'Oxidative phosphorylation', 'HALLMARK_MTORC1_SIGNALING': 'mTORC1 signaling'}
def lab(k):
    return LABELS.get(k, str(k).replace('HALLMARK_', '').replace('_', ' ').capitalize()).replace('Myc', 'MYC').replace('E2f', 'E2F').replace('G2m', 'G2M').replace('Uv ', 'UV ').replace('Dna ', 'DNA ').replace('Kras', 'KRAS')

def clean(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(pad=3)

def inner(ax, box):
    b = ax.get_position()
    return ax.figure.add_axes([b.x0 + box[0] * b.width, b.y0 + box[1] * b.height, box[2] * b.width, box[3] * b.height])


HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=HERE/'source_data');ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
N=args.source_dir;SO=args.output_dir;SO.mkdir(parents=True,exist_ok=True);PANELS=[]
EXT=['GSE63042','GSE95233','VANISH']
MEM=['HALLMARK_'+p for p in ['HYPOXIA','MTORC1_SIGNALING','GLYCOLYSIS','CHOLESTEROL_HOMEOSTASIS','HEME_METABOLISM','G2M_CHECKPOINT','MITOTIC_SPINDLE','REACTIVE_OXYGEN_SPECIES_PATHWAY','ESTROGEN_RESPONSE_LATE','EPITHELIAL_MESENCHYMAL_TRANSITION','XENOBIOTIC_METABOLISM','ANGIOGENESIS','P53_PATHWAY','ADIPOGENESIS','NOTCH_SIGNALING','WNT_BETA_CATENIN_SIGNALING']]
LABELS.update({'Score16':'Score16','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION':'Epithelial–mesenchymal transition','HALLMARK_ESTROGEN_RESPONSE_LATE':'Estrogen response late','HALLMARK_P53_PATHWAY':'p53 pathway'})
def load(name,legacy=False):return pd.read_csv(N/name)
def record(panel,d,source,note=''):
 PANELS.append(dict(panel=panel,source=source,notes=note,records=d.to_dict(orient='records')))
def export(fig,stem):
 fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[]
 for text in fig.findobj(mpl.text.Text):
  if not text.get_visible() or not text.get_text():continue
  b=text.get_window_extent(renderer)
  if b.x0<-.8 or b.y0<-.8 or b.x1>fig.bbox.x1+.8 or b.y1>fig.bbox.y1+.8:bad.append(text.get_text())
 assert not bad,('clipped text',bad)
 fig.savefig(str(stem)+'.pdf',dpi=300);fig.savefig(str(stem)+'.png',dpi=300);fig.savefig(str(stem)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});fig.savefig(str(stem)+'_preview.png',dpi=160)
def row(d,**f):
 for k,v in f.items():d=d[d[k].eq(v)]
 assert len(d)==1,(f,len(d));return d.iloc[0]

def ref(ax,at=0):ax.axvline(at,color=AXIS_GRAY,lw=.7,ls=(0,(4,3)),zorder=0)

def ticks_log(ax,limits):
 ax.set_xscale('log');ax.set_xlim(*limits);ts=[v for v in [.25,.5,1,2,4,8,16,32] if limits[0]<=v<=limits[1]];ax.set_xticks(ts,[str(v) for v in ts]);ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())

def forest(ax,d,labels,diamond=None,pi=None,xlabel='OR per 1-SD Score16',xlim=None,numbers=False,color=DEEP_TEAL):
 assert len(d)==len(labels);ys=np.arange(len(d))[::-1];ref(ax,1)
 for i,(y,r) in enumerate(zip(ys,d.itertuples())):
  lo,hi,x=float(r.CI_low),float(r.CI_high),float(r.OR);assert 0<lo<=x<=hi
  if i==pi:ax.plot([lo,hi],[y,y],color=WARM_TAUPE,lw=1.1);ax.vlines([lo,hi],y-.08,y+.08,color=WARM_TAUPE,lw=.7)
  elif i==diamond:ax.add_patch(Polygon([[lo,y],[x,y+.14],[hi,y],[x,y-.14]],facecolor=color,edgecolor=color,lw=.5,zorder=3))
  else:ax.errorbar(x,y,xerr=[[x-lo],[hi-x]],fmt='o',ms=3.7,color=color,mec='white',mew=.35,lw=.9,capsize=2)
  if numbers:ax.text(1.03,y,f'{x:.2f} [{lo:.2f}, {hi:.2f}]' if i!=pi else f'{lo:.2f}–{hi:.2f}',transform=ax.get_yaxis_transform(),va='center',fontsize=6.3)
 if xlim is None:xlim=(min(.7,float(d.CI_low.min())*.85),float(d.CI_high.max())*1.12)
 assert d.CI_low.min()>=xlim[0] and d.CI_high.max()<=xlim[1],('CI_CLIPPED',d,xlim)
 ticks_log(ax,xlim);ax.set_yticks(ys,labels,fontsize=6.5);ax.tick_params(axis='y',length=0);ax.set_ylim(-.6,len(d)-.4);ax.set_xlabel(xlabel);clean(ax)

def errorpoints(ax,x,lo,hi,y,color=DEEP_TEAL,marker='o',label=None,filled=True,ms=3.6):
 x,lo,hi=map(np.asarray,[x,lo,hi]);ax.errorbar(x,y,xerr=[x-lo,hi-x],fmt=marker,mfc=color if filled else 'white',mec='white' if filled else color,mew=.4,ecolor=color,ms=ms,lw=.8,capsize=1.7,label=label)

def meta(group='external'):
 return row(load('05 Score16 Adult Meta and Sensitivity.csv'),axis='Score16',group=group,analysis=('primary_three_units' if group=='external' else 'primary_discovery_four_units'),omitted='none')

def f3a(ax):
 d=load('04 Score16 Adult Cohort Results.csv');d=d[d.axis.eq('Score16')&d.cohort.isin(EXT)].set_index('cohort').loc[EXT].reset_index();r=meta();pi=r.copy();pi['CI_low']=r.PI_low;pi['CI_high']=r.PI_high;d=pd.concat([d,pd.DataFrame([r,pi])],ignore_index=True);forest(ax,d,EXT+['Pooled · HKSJ','Prediction interval'],diamond=3,pi=4,xlim=(.65,7.2),numbers=True);ax.set_title('Three adult external study units',fontsize=7.5);record('Figure3a',d,'04 cohort +05 external meta','333 patients/93 deaths; three independent retained study units; PI separate row')

def conditional_models():
 v=load('06 Score16 VANISH CTS SRS HiDEF.csv');m=load('07 Score16 MARS.csv');d=pd.concat([v,m]);return d[d.axis.eq('Score16')&~d.clinical_adjustment].set_index('benchmark').loc[['CTS','SRS','MARS','HiDEF']].reset_index()

def f3b(ax):
 d=conditional_models();forest(ax,d,['CTS · VANISH','SRS · VANISH','MARS · GSE65682 (discovery)','Hi-DEF · VANISH'],xlim=(.65,3.85));ax.set_title('Score16 + established state',fontsize=7);record('Figure3b',d,'06 VANISH +07 MARS','MARS discovery; others reuse VANISH; no pooled effect; MARS discovery')

def f3c(ax):
 r=row(load('06 Score16 VANISH CTS SRS HiDEF.csv'),axis='Score16',analysis='CTS_unadjusted');errorpoints(ax,[r.AUC_base],[r.AUC_base_low],[r.AUC_base_high],[1],WARM_TAUPE,'^',filled=False);errorpoints(ax,[r.AUC_full],[r.AUC_full_low],[r.AUC_full_high],[0]);ax.set_yticks([1,0],['CTS only','CTS + Score16'],fontsize=6.5);ax.set(xlim=(.49,.83),ylim=(-1.1,1.55),xlabel='Within-cohort apparent AUC');ax.set_xticks([.5,.6,.7,.8]);ax.text(.98,.06,f'ΔAUC = {r.delta_AUC:.3f}\n95% CI {r.delta_AUC_low:.3f}–{r.delta_AUC_high:.3f}',transform=ax.transAxes,ha='right',va='bottom',fontsize=6.2);ax.tick_params(axis='y',length=0);clean(ax);record('Figure3c',pd.DataFrame([r]),'06 Score16 VANISH CTS SRS HiDEF.csv','paired bootstrap2000; apparent, not external incremental gain')

def f3d(ax):
 d=load('10f External Portability Models.csv');d=d[d.axis.eq('Score16')].set_index('cohort').loc[EXT].reset_index();y=np.arange(3)[::-1];ax.hlines(y,d.original_AUC,d.portable_AUC,color=LIGHT,lw=1);ax.scatter(d.original_AUC,y,s=16,color=DEEP_TEAL,edgecolor='white',lw=.35,label='ssGSEA Score16',zorder=3);ax.scatter(d.portable_AUC,y,s=18,marker='s',facecolor='white',edgecolor=DUSTY_BLUE,lw=.8,label='Portable rank',zorder=3);ax.set_yticks(y,EXT,fontsize=6.5);ax.set(xlim=(.53,.80),ylim=(-.5,3.45),xlabel='Mortality AUC');ax.set_xticks([.55,.65,.75]);ax.legend(loc='upper center',ncol=1,frameon=False,fontsize=6,bbox_to_anchor=(.55,1.08));ax.tick_params(axis='y',length=0);clean(ax);record('Figure3d',d,'10f External Portability Models.csv','no new fit; discovery-only reference parameters')


def panel(fig,box,fn,letter,pos):a=fig.add_axes(box);fn(a);fig.text(*pos,letter,fontsize=10,fontweight='bold',va='top');return a

def portable_delta(ax):
 d=load('10f External Portability Models.csv');d=d[d.axis.eq('Score16')].set_index('cohort').loc[EXT].reset_index();errorpoints(ax,d.delta_AUC,d.delta_AUC_low,d.delta_AUC_high,np.arange(3)[::-1],DUSTY_BLUE,'s',filled=False);ref(ax);ax.set(yticks=np.arange(3)[::-1],yticklabels=EXT,xlabel='Portable − ssGSEA paired ΔAUC',ylim=(-.6,2.6));ax.tick_params(axis='y',length=0);clean(ax);record('S10a',d,'10f External Portability Models.csv','same patients paired1000; not gain vs clinical model')

def portable_or(ax):
 d=load('10f External Portability Models.csv');d=d[d.axis.eq('Score16')].set_index('cohort').loc[EXT].reset_index();forest(ax,d,EXT,xlabel='OR per discovery-reference SD');ax.set_title('Portable rank score',fontsize=7.3);record('S10b',d,'10f External Portability Models.csv','unit differs from within-cohort SD; no pooling')

def portable_correlations(ax):
 d=load('10f External Portability Models.csv');d=d[d.axis.eq('Score16')].set_index('cohort').loc[EXT].reset_index();y=np.arange(3)[::-1];ax.hlines(y,0,d.rho_original,color=LIGHT,lw=1);ax.scatter(d.rho_original,y,s=16,color=DUSTY_BLUE,edgecolor='white',lw=.35)
 for yy,x in zip(y,d.rho_original):ax.annotate(f'{x:.3f}',xy=(x,yy),xytext=(5,0),textcoords='offset points',va='center',fontsize=6.3)
 ax.set(yticks=y,yticklabels=EXT,xlim=(0,1.14),ylim=(-.6,2.6),xlabel='Spearman ρ with ssGSEA Score16');ax.tick_params(axis='y',length=0);clean(ax);record('S10c',d,'10f External Portability Models.csv')

def portable_coverage(ax):
 d=load('10a Portable Gene Coverage.csv');ids=d.matrix.unique().tolist();p=d[d.pathway.isin(MEM)].pivot(index='pathway',columns='matrix',values='coverage').loc[MEM,ids];im=ax.imshow(p,aspect='auto',cmap=mpl.colors.LinearSegmentedColormap.from_list('coverage',['#F5F1E9',DUSTY_BLUE]),vmin=0,vmax=1);ax.set_yticks(range(16),[lab(x) for x in MEM],fontsize=6.0);ax.set_xticks(range(len(ids)),ids,rotation=55,ha='right',fontsize=5.8);ax.tick_params(length=0)
 for sp in ax.spines.values():sp.set_visible(False)
 cax=inner(ax,[1.035,.22,.025,.50]);cb=ax.figure.colorbar(im,cax=cax);cb.set_ticks([0,.5,1]);cb.ax.set_title('Coverage',fontsize=6,pad=6);cb.ax.tick_params(labelsize=6,length=2);record('S10d',d,'10a Portable Gene Coverage.csv','seven expression matrices: four unchanged discovery references and three external units')

def s10():
 fig=plt.figure(figsize=(170/25.4,207/25.4));panel(fig,[.23,.76,.29,.15],portable_delta,'a',(.02,.98));panel(fig,[.74,.76,.21,.15],portable_or,'b',(.60,.98));panel(fig,[.23,.52,.29,.13],portable_correlations,'c',(.02,.70));panel(fig,[.38,.13,.53,.28],portable_coverage,'d',(.02,.46));export(fig,SO/'Figure_S10');plt.close(fig)






def figure3():
 fig=plt.figure(figsize=(170/25.4,208/25.4))
 boxes=[(.20,.665,.60,.27),(.22,.38,.28,.17),(.61,.38,.29,.17),(.22,.075,.68,.17)]
 letters=[(.02,.97),(.02,.605),(.59,.605),(.02,.285)]
 for i,(fn,box,pos) in enumerate(zip([f3a,f3b,f3c,f3d],boxes,letters)):
  fn(fig.add_axes(box));fig.text(*pos,chr(97+i),fontsize=10,fontweight='bold',ha='left',va='top')
 export(fig,SO/'Figure3');plt.close(fig)
figure3()
(SO/'figure_data_ledger.json').write_text(json.dumps(PANELS,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
print('Rendered Figure3, S10 from approved source functions and actual CSVs; no models fitted.')

