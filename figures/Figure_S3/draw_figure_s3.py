"""Figure S3 from frozen results only; no model fitting or pathway selection.

Run with Python and installed numpy/pandas/matplotlib:
    python draw_figure_s3.py --output-dir /path/to/new/outputs
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.ticker import NullLocator

HERE=Path(__file__).resolve().parent
DISC=['GSE272769','GSE65682','GSE185263','GSE236713']
# Palette, typography and line widths from the approved Figure 1–5 renderer.
TEAL='#3F827C';BLUE='#7298B6';GRAY='#77716C';INK='#333333'
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,
 'axes.titlesize':7.5,'axes.titleweight':'normal','axes.linewidth':.7,
 'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':2.5,
 'ytick.major.size':2.5,'legend.fontsize':6.3,'pdf.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white','text.color':INK,
 'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})

def read(name):return pd.read_csv(HERE/'source_data'/name)
def setup(ax,ylabels):
 ax.axvline(1,color=GRAY,lw=.7,ls=(0,(4,3)),zorder=0)
 ax.set_xscale('log');ax.set_xlim(.75,3.5)
 ax.set_xticks([1,1.5,2,3],['1','1.5','2','3'])
 ax.xaxis.set_minor_locator(NullLocator())
 ax.set_yticks(np.arange(len(ylabels))[::-1],ylabels,fontsize=6.5)
 ax.set_ylim(-.55,len(ylabels)-.45)
 ax.set_xlabel('OR per 1-SD score',labelpad=5)
 ax.tick_params(axis='y',length=0,pad=4)
 ax.spines[['top','right']].set_visible(False)

def check_ci(d):
 assert all(d.CI_low>0) and all(d.CI_low<=d.OR) and all(d.OR<=d.CI_high)
 assert d.CI_low.min()>=.75 and d.CI_high.max()<=3.5,'CI outside plotting bounds'

def paired(ax,primary,other,other_label):
 primary=primary.set_index('cohort').loc[DISC].reset_index()
 other=other.set_index('cohort').loc[DISC].reset_index()
 setup(ax,DISC);y=np.arange(4)[::-1]
 for d,off,color,marker,label,filled in [(primary,.12,TEAL,'o','Primary Score16',True),
    (other,-.12,BLUE,'s',other_label,False)]:
  check_ci(d)
  ax.errorbar(d.OR,y+off,xerr=[d.OR-d.CI_low,d.CI_high-d.OR],fmt=marker,
   ms=3.5,mfc=color if filled else 'white',mec='white' if filled else color,
   mew=.4,ecolor=color,lw=.8,capsize=1.7,label=label)
 ax.legend(loc='lower right',bbox_to_anchor=(1.02,1.02),frameon=False,
           fontsize=6.1,borderaxespad=0,handlelength=1.6,labelspacing=.35)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
 out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
 meta=read('discovery fixedscore meta.csv')
 p=meta[meta.analysis.eq('primary_discovery_four_units')]
 omissions=meta[meta.analysis.eq('leave_one_study_out')].set_index('omitted').loc[DISC].reset_index()
 assert len(p)==1 and len(omissions)==4 and set(meta.axis)=={'Score16'}
 d=pd.concat([p,omissions],ignore_index=True)
 check_ci(d)
 collapsed_meta=read('cell cycle collapsed meta.csv')
 assert len(collapsed_meta)==1 and collapsed_meta['axis'].iloc[0]=='Score16_cell_cycle_collapsed'
 check_ci(collapsed_meta)  # Supplied pooled result is a traceability check; no model refit.
 primary=read('discovery cohort effects.csv')
 if 'axis' in primary:primary=primary[primary.axis.eq('Score16')]
 fig=plt.figure(figsize=(170/25.4,145/25.4))
 a=fig.add_axes([.23,.635,.52,.265])
 setup(a,['Full discovery']+['Omit '+v for v in DISC])
 y=np.arange(5)[::-1]
 for i,r in d.iterrows():
  a.errorbar(r.OR,y[i],xerr=[[r.OR-r.CI_low],[r.CI_high-r.OR]],
   fmt='D' if i==0 else 'o',ms=3.8 if i==0 else 3.5,color=TEAL,
   mec='white',mew=.35,lw=.8,capsize=1.8)
  a.text(1.035,y[i],f'{r.OR:.2f} [{r.CI_low:.2f}, {r.CI_high:.2f}]',
         transform=a.get_yaxis_transform(),va='center',fontsize=6.4)
 b=fig.add_axes([.21,.12,.28,.295]);c=fig.add_axes([.71,.12,.27,.295])
 paired(b,primary,read('cell cycle collapsed.csv'),'Cell-cycle collapsed')
 paired(c,primary,read('qc exclusion.csv'),'QC exclusions')
 for letter,x,y,title in [('a',.02,.963,'Fixed-score discovery-cohort omission'),
    ('b',.02,.535,'Cell-cycle de-redundancy'),('c',.53,.535,'QC exclusion')]:
  fig.text(x,y,letter,fontsize=10,fontweight='bold',va='top')
  fig.text(.23 if letter=='a' else .21 if letter=='b' else .71,y,title,fontsize=7.5,va='top')
 fig.canvas.draw();renderer=fig.canvas.get_renderer();clipped=[]
 for text in fig.findobj(mpl.text.Text):
  if not text.get_visible() or not text.get_text():continue
  bb=text.get_window_extent(renderer)
  if bb.x0<-1 or bb.y0<-1 or bb.x1>fig.bbox.width+1 or bb.y1>fig.bbox.height+1:clipped.append(text.get_text())
 assert not clipped,('Text outside canvas',clipped)
 for ext in ('pdf','png','svg'):fig.savefig(out/f'Figure_S3.{ext}',dpi=300)
 inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((HERE/'source_data').glob('*.csv'))}
 qc={'panels':['a: fixed-score discovery omission','b: cell-cycle de-redundancy','c: QC exclusion'],
     'fit_or_resampling_performed':False,'inputs_sha256':inputs,'text_outside_canvas':clipped,
     'omission_ORs':dict(zip(omissions.omitted,omissions.OR)),
     'font':'Arial','size_mm':[170,145],'point_size':3.5,'line_width':.8}
 (out/'Figure_S3_render_QC.json').write_text(json.dumps(qc,indent=2),encoding='utf-8')
 print(json.dumps(qc,indent=2));plt.close(fig)

if __name__=='__main__':main()
