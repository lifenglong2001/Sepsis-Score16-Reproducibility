from pathlib import Path
import numpy as np,pandas as pd,matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
ROOT=Path(__file__).resolve().parent
S=ROOT/'source data'
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

def f1a(ax):
    ax.axis('off')
    for x, name, detail, c in [(0, 'Diagnosis', 'GSE236713 · GSE63042', DUSTY_BLUE), (0.35, 'Severity', 'GSE185263 · GSE154918', SAGE), (0.7, 'Mortality', '4 discovery → 3 adult external', DEEP_TEAL)]:
        ax.plot([x, x + 0.27], [0.91, 0.91], color=c, lw=1.8)
        ax.text(x, 0.78, name, fontweight='bold', color=c, va='top')
        ax.text(x, 0.62, detail, fontsize=6.5, va='top')
    ax.text(0, 0.31, '50 Hallmarks', fontweight='bold')
    ax.annotate('', xy=(0.34, 0.34), xytext=(0.22, 0.34), arrowprops={'arrowstyle': '->', 'color': AXIS_GRAY, 'lw': 0.8})
    ax.text(0.35, 0.31, 'Cohort-specific effects')
    ax.annotate('', xy=(0.69, 0.34), xytext=(0.61, 0.34), arrowprops={'arrowstyle': '->', 'color': AXIS_GRAY, 'lw': 0.8})
    ax.text(0.7, 0.31, 'Score16', fontweight='bold', color=DEEP_TEAL)
    ax.text(0, 0.03, 'Cellular context: whole blood / neutrophils → PBMC', fontsize=6.5)
    ax.set(xlim=(-0.01, 1.01), ylim=(-0.06, 1))

def f1b(ax):
    d = fig1[fig1.panel.eq('B')]
    for cohort, marker in [('GSE236713', 'o'), ('GSE63042', 's')]:
        v = d[d.cohort.eq(cohort)].set_index('clinical_task').loc[['Diagnosis', 'Mortality']]
        ax.plot([0, 1], v.estimate, color=LIGHT, lw=0.8, zorder=1)
        for x, (_, r), c in zip([0, 1], v.iterrows(), [DUSTY_BLUE, DEEP_TEAL]):
            ax.errorbar(x, r.estimate, yerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]], fmt=marker, color=c, mfc=c, mec='white', ms=4, capsize=2, lw=0.9, zorder=3)
    ax.axhline(0.5, color=AXIS_GRAY, ls=(0, (4, 3)), lw=0.7)
    ax.set(xticks=[0, 1], xticklabels=['Diagnosis', 'Mortality'], ylim=(0.2, 1.03), xlim=(-0.15, 1.15), ylabel='Strict OOF AUC')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=INK, mec='white', label='GSE236713'), Line2D([], [], marker='s', ls='', mfc=INK, mec='white', label='GSE63042')], loc='upper right', frameon=False)
    clean(ax)

def f1c(ax):
    d = fig1[fig1.element.eq('full-adjusted SOFA forest')].set_index('pathway').loc[['HALLMARK_GLYCOLYSIS', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY', 'HALLMARK_IL6_JAK_STAT3_SIGNALING']]
    ys = np.arange(3)[::-1]
    ax.errorbar(d.estimate, ys, xerr=[d.estimate - d.ci_low, d.ci_high - d.estimate], fmt='o', mfc=SAGE, mec='white', ecolor=SAGE, ms=4, lw=0.9, capsize=2)
    ax.axvline(0, color=AXIS_GRAY, ls='--', lw=0.7)
    ax.set_yticks(ys, [lab(k).replace('Oxidative phosphorylation', 'Oxidative\nphosphorylation') for k in d.index])
    ax.set(xlabel='Adjusted SOFA points per 1-SD score', xlim=(-0.4, 1.45), ylim=(-0.6, 2.6))
    ax.tick_params(axis='y', length=0)
    clean(ax)

def f1d(parent):
    v = fig1[fig1.panel.eq('D')].set_index('element')
    KEYS = ['IL6 cross-sectional severity', 'IL6 adjusted SOFA', 'IL6 mortality discovery meta', 'IL6 external mortality']
    SAGE_D = '#93A698'

    def mini(parent, box):
        b = parent.get_position()
        return parent.figure.add_axes([b.x0 + box[0] * b.width, b.y0 + box[1] * b.height, box[2] * b.width, box[3] * b.height])

    def clean(ax):
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.spines['bottom'].set_linewidth(0.7)
        ax.set_yticks([])
        ax.tick_params(axis='x', labelsize=7, pad=2, length=2.5, width=0.7)
    parent.axis('off')
    parent.text(0.02, 1.05, 'IL-6/JAK/STAT3', fontsize=7.5, fontweight='bold', va='baseline')
    for x, label in [(0.155, 'Ordered severity'), (0.458, 'Adjusted 24-h SOFA'), (0.837, 'Mortality')]:
        parent.text(x, 0.88, label, ha='center', va='center', fontsize=7)
    a = mini(parent, [0.025, 0.12, 0.255, 0.57])
    a.axvline(0, ymax=0.66, color=AXIS_GRAY, lw=0.65, ls=(0, (3, 3)), zorder=0)
    r = v.loc[KEYS[0]]
    a.plot(r.estimate, 0, 'o', ms=4.4, color=SAGE_D, zorder=3)
    a.text(0.5, 0.88, f'ρ = {r.estimate:.3f}', transform=a.transAxes, ha='center', fontsize=7)
    a.set(xlim=(-1, 1), ylim=(-0.65, 1.35), xticks=[-1, 0, 1], xlabel='Spearman ρ')
    a.xaxis.label.set_fontsize(7)
    clean(a)
    b = mini(parent, [0.337, 0.12, 0.245, 0.57])
    b.axvline(0, ymax=0.66, color=AXIS_GRAY, lw=0.65, ls=(0, (3, 3)), zorder=0)
    r = v.loc[KEYS[1]]
    b.errorbar(r.estimate, 0, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]], fmt='o', ms=4.4, color=SAGE_D, lw=0.9, capsize=2, zorder=3)
    b.text(0.5, 0.88, f'β = {r.estimate:.3f} [{r.ci_low:.3f}, {r.ci_high:.3f}]', transform=b.transAxes, ha='center', fontsize=6.5)
    b.set(xlim=(-0.4, 0.7), ylim=(-0.65, 1.35), xticks=[-0.4, 0, 0.4], xlabel='SOFA β per 1-SD score')
    b.xaxis.label.set_fontsize(7)
    clean(b)
    c = mini(parent, [0.7, 0.12, 0.135, 0.57])
    c.axvline(1, color=AXIS_GRAY, lw=0.65, ls=(0, (3, 3)), zorder=0)
    for key, y, marker, label in [(KEYS[2], 1, 'D', 'Discovery'), (KEYS[3], 0, 'o', 'GSE63042')]:
        r = v.loc[key]
        c.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]], fmt=marker, ms=4 if marker == 'D' else 4.4, color=DEEP_TEAL, lw=0.9, capsize=2, zorder=3)
        parent.text(0.61, 0.12 + 0.57 * (y + 0.45) / 2.05, label, transform=parent.transAxes, ha='left', va='center', fontsize=6.5)
        parent.text(1.015, 0.12 + 0.57 * (y + 0.45) / 2.05, f'{r.estimate:.3f} [{r.ci_low:.3f}, {r.ci_high:.3f}]', transform=parent.transAxes, ha='right', va='center', fontsize=6.5, zorder=5)
    c.set(xscale='log', xlim=(0.4, 1.55), ylim=(-0.45, 1.6), xlabel='OR per 1-SD score')
    c.set_xticks([0.5, 1, 1.5], ['0.5', '1', '1.5'])
    c.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    c.xaxis.label.set_fontsize(7)
    clean(c)
    return [a, b, c]

import argparse
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();OUT=args.output_dir;OUT.mkdir(parents=True,exist_ok=True)
fig1=read('Figure1 source.csv')
def export(fig,name):
 fig.canvas.draw();r=fig.canvas.get_renderer();bad=[]
 for t in fig.findobj(mpl.text.Text):
  if t.get_visible() and t.get_text():
   b=t.get_window_extent(r)
   if b.x0<-.8 or b.y0<-.8 or b.x1>fig.bbox.x1+.8 or b.y1>fig.bbox.y1+.8:bad.append(t.get_text())
 assert not bad,('clipped text',bad)
 for ext in ['png','pdf','tiff']:
  kw={'dpi':600 if ext=='tiff' else 300}
  if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(OUT/(name+'.'+ext),**kw)
 plt.close(fig)
fig=plt.figure(figsize=(170/25.4,159/25.4))
for fn,box,letter,pos in [(f1a,[.07,.79,.86,.15],'a',[.02,.98]),(f1b,[.10,.43,.25,.24],'b',[.02,.735]),(f1c,[.67,.43,.28,.24],'c',[.47,.735]),(f1d,[.06,.085,.84,.18],'d',[.02,.32])]:
 fn(fig.add_axes(box));fig.text(*pos,letter,fontsize=10,fontweight='bold',va='top')
export(fig,'Figure1')
data=read('S4 diagnostic OOF source.csv')
fig=plt.figure(figsize=(170/25.4,90/25.4))
for j,cohort in enumerate(['GSE236713','GSE63042']):
 ax=fig.add_axes([.19+j*.49,.24,.27,.53]);q=data[data.cohort.eq(cohort)]
 assert len(q)==2
 for y,(_,r),color in zip([1,0],q.iterrows(),[DUSTY_BLUE,DEEP_TEAL]):
  ax.errorbar(r.AUC,y,xerr=[[r.AUC-r.AUC_CI_low],[r.AUC_CI_high-r.AUC]],fmt='o',color=color,mec='white',ms=4,capsize=2,lw=.9)
  ax.text(.98,(y+.35)/2,f'n = {int(r.n)}',transform=ax.transAxes,ha='right',fontsize=6.3)
 ax.axvline(.5,color=AXIS_GRAY,lw=.7,ls=(0,(4,3)));ax.set(xlim=(.25,1),ylim=(-.6,1.6),yticks=[1,0],yticklabels=['Diagnosis','Mortality'],xlabel='Unchanged diagnostic OOF score: AUC',title=cohort)
 ax.xaxis.label.set_fontsize(6.3);ax.tick_params(axis='y',length=0);clean(ax)
 fig.text(.02+j*.49,.96,chr(97+j),fontsize=10,fontweight='bold',va='top')
export(fig,'Figure S4')
