from pathlib import Path
import argparse,json,urllib.request,hashlib
p=argparse.ArgumentParser(description="Acquire public prerequisites only; no analysis or raw-to-paper reconstruction.")
p.add_argument('--destination',required=True);p.add_argument('--plan',action='store_true');p.add_argument('--cts',action='store_true');p.add_argument('--geo-metadata',action='store_true');a=p.parse_args()
root=Path(a.destination);tasks=[]
commit="467fa897ba33b8ec355a8f33a1d3c4f4486846bd"
if a.cts or a.plan:
 for name in ['exp.core.df.txt','coreSamples_CTS_markup.csv','exp.vanish.class.txt']:
  tasks.append((name,'https://raw.githubusercontent.com/bpscicluna/Consensus-Transcriptomic-Subtypes-of-Sepsis/'+commit+'/'+name))
if a.geo_metadata or a.plan:
 for g in ['GSE272769','GSE65682','GSE185263','GSE236713','GSE63042','GSE95233','GSE154918','GSE216009','GSE279452','GSE290679']:
  tasks.append((g+'_family.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/'+g[:-3]+'nnn/'+g+'/soft/'+g+'_family.soft.gz'))
if a.plan or not tasks:
 print(json.dumps(tasks,indent=2));raise SystemExit(0)
root.mkdir(parents=True,exist_ok=True);manifest=[]
for name,url in tasks:
 target=root/name
 if target.exists():raise FileExistsError('Refusing to overwrite '+str(target))
 req=urllib.request.Request(url,headers={'User-Agent':'Score16 prerequisite acquisition'})
 with urllib.request.urlopen(req,timeout=120) as h:data=h.read()
 if not data:raise RuntimeError('Empty response: '+url)
 target.write_bytes(data);manifest.append({'file':name,'source':url,'sha256':hashlib.sha256(data).hexdigest()})
(root/'Public acquisition manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
