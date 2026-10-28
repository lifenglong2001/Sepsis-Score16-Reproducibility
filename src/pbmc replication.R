# Independent donor-level PBMC replication in GSE290679 (T9).
# Donor is the experimental unit; no cell-level pseudo-replication.
# Frozen inputs: donor x lineage pseudobulk counts, donor metadata, and the
# predefined direction registry (GSE279452 original Cliff's deltas).
# Upstream boundary: the author-annotated h5ad -> pseudobulk aggregation is NOT
# part of this stage; it is provided separately as scripts/upstream/.
#
# Scope: five prespecified metabolic/stress Hallmarks x five PBMC lineages
# (25 lineage-pathway comparisons per comparator).
suppressPackageStartupMessages({library(edgeR);library(GSVA);library(BiocParallel)})
stopifnot(as.character(packageVersion('GSVA'))=='2.4.9')

set.seed(20260928)
con<-gzfile(cfg$inputs$gse290679_pseudobulk_counts,'rt')
counts<-read.delim(con,check.names=FALSE)
close(con)
rownames(counts)<-counts[[1]]
counts<-as.matrix(counts[,-1])
storage.mode(counts)<-'double'
meta<-read.csv(cfg$inputs$gse290679_pseudobulk_metadata,stringsAsFactors=FALSE)
stopifnot(identical(colnames(counts),meta$pseudobulk_id))
stopifnot(setequal(unique(meta$group),c('CI-Sep','CI-NS','NHC')))
# 19 sepsis ICU, 19 non-septic ICU and 9 healthy donors, each as donor x lineage pseudobulks.
donors<-unique(meta[,c('donor_id','group')])
stopifnot(identical(as.integer(table(factor(donors$group,levels=c('CI-Sep','CI-NS','NHC')))),c(19L,19L,9L)))

gmt<-strsplit(readLines(cfg$inputs$gmt,warn=FALSE),'\t')
sets<-setNames(lapply(gmt,function(x)unique(x[-c(1,2)])),vapply(gmt,'[',character(1),1))
# The five prespecified metabolic/stress Hallmarks for GSE290679 replication.
stress<-paste0('HALLMARK_',c('HYPOXIA','MTORC1_SIGNALING','GLYCOLYSIS',
  'CHOLESTEROL_HOMEOSTASIS','REACTIVE_OXYGEN_SPECIES_PATHWAY'))
stopifnot(all(stress%in%names(sets)))
sets<-sets[stress]

types<-c('Monocyte','CD4 T','CD8 T','B','NK')
scores<-list()
for(ct in types){
 ii<-meta$cell_type==ct & meta$cell_count>=20
 m<-meta[ii,]
 y<-DGEList(counts[,ii,drop=FALSE])
 keep<-rowSums(cpm(y)>=1)>=min(5,ncol(y))
 y<-y[keep,,keep.lib.sizes=FALSE]
 y<-calcNormFactors(y,method='TMM')
 expr<-cpm(y,log=TRUE,prior.count=2)
 variable<-rownames(expr)[apply(expr,1,var)>0]
 below<-vapply(sets,function(g)length(intersect(g,variable))<10,logical(1))
 stopifnot(!any(below))
 par<-ssgseaParam(exprData=expr,geneSets=sets,minSize=10,maxSize=Inf,alpha=.25,normalize=TRUE)
 s<-gsva(par,verbose=FALSE,BPPARAM=SerialParam())
 stopifnot(nrow(s)==5,all(is.finite(s)))
 z<-cbind(m,as.data.frame(t(s)))
 scores[[ct]]<-z
}
tab<-do.call(rbind,scores)

delta<-function(x,y){d<-outer(x,y,'-');mean(d>0)-mean(d<0)}
rows<-list()
metrics<-names(sets)
for(ct in types)for(control in c('CI-NS','NHC'))for(metric in metrics){
 z<-tab[tab$cell_type==ct,]
 x<-z[z$group=='CI-Sep',metric]
 y<-z[z$group==control,metric]
 enough<-length(x)>=3 & length(y)>=3
 if(enough){
  ci<-quantile(replicate(2000,delta(sample(x,replace=TRUE),sample(y,replace=TRUE))),
               c(.025,.975),type=6)
  tt<-wilcox.test(x,y,exact=FALSE)
  effect<-delta(x,y)
 }else{ci<-c(NA,NA);tt<-list(p.value=NA);effect<-NA}
 rows[[length(rows)+1]]<-data.frame(
  contrast=if(control=='CI-NS')'Sepsis vs ICU non-sepsis' else 'Sepsis vs Healthy',
  cell_type=ct,metric=metric,n_sepsis=length(x),n_control=length(y),
  mean_sepsis=mean(x),mean_control=mean(y),median_difference=median(x)-median(y),
  cliff_delta=effect,ci_lower=ci[1],ci_upper=ci[2],P=tt$p.value,test_performed=enough,
  stringsAsFactors=FALSE)
}
res<-do.call(rbind,rows)
res$FDR<-NA_real_
for(cn in unique(res$contrast)){
 ii<-res$contrast==cn
 res$FDR[ii]<-p.adjust(res$P[ii],'BH')
}
# Read the direction registry and keep only the five stress pathways
target<-read.csv(cfg$inputs$gse290679_direction_registry,stringsAsFactors=FALSE)
target<-target[target$metric%in%stress,c('contrast','cell_type','metric','cliff_delta')]
names(target)[4]<-'original_cliff_delta'
res<-merge(res,target,by=c('contrast','cell_type','metric'),all.x=TRUE,sort=FALSE)
stopifnot(nrow(res)==50L,!anyNA(res$original_cliff_delta))
res$direction_concordant<-sign(res$cliff_delta)==sign(res$original_cliff_delta)
res$same_direction_supported<-res$direction_concordant & res$FDR<.05 & (res$ci_lower>0 | res$ci_upper<0)

summ<-do.call(rbind,lapply(unique(res$contrast),function(cn){
 z<-res[res$contrast==cn,]
 data.frame(contrast=cn,predefined_targets=nrow(z),same_direction=sum(z$direction_concordant),
            same_direction_supported=sum(z$same_direction_supported),
            supported_lineages=length(unique(z$cell_type[z$same_direction_supported])),
            stringsAsFactors=FALSE)
}))

fwrite(as.data.table(res),file.path(OUT,'gse290679_replication_effects.csv'))
fwrite(as.data.table(summ),file.path(OUT,'gse290679_replication_summary.csv'))

# ---- numeric comparison against the frozen T9 tables -------------------------
e<-as.data.table(expected('replication'))
stopifnot(nrow(e)==50L)
stopifnot(nrow(res)==50L)
setkey(e,contrast,cell_type,metric)
r<-as.data.table(res);setkey(r,contrast,cell_type,metric)
stopifnot(identical(e[,.(contrast,cell_type,metric)],r[,.(contrast,cell_type,metric)]))
for(nm in c('n_sepsis','n_control','mean_sepsis','mean_control','median_difference',
            'cliff_delta','ci_lower','ci_upper','P','FDR','original_cliff_delta')){
 for(i in seq_len(nrow(r))){
  check('C_GSE290679_REPLICATION',paste(r$contrast[i],r$cell_type[i],sep='|'),nm,
        r[[nm]][i],e[[nm]][i],
        cfg$tolerance,inputs=c('gse290679_pseudobulk_counts','gse290679_pseudobulk_metadata',
                               'gse290679_direction_registry','gmt'),
        n=r$n_sepsis[i]+r$n_control[i],
        estimand='donor-level Cliff delta; 2000 unstratified group bootstraps; BH over 25 prespecified comparisons')
 }
}
for(nm in c('direction_concordant','same_direction_supported')){
 for(i in seq_len(nrow(r))){
  check('C_GSE290679_REPLICATION',paste(r$contrast[i],r$cell_type[i],sep='|'),nm,
        as.integer(r[[nm]][i]),as.integer(e[[nm]][i]),0,
        inputs='gse290679_direction_registry',n=r$n_sepsis[i]+r$n_control[i],
        estimand='direction vs frozen GSE279452 Cliff delta; FDR<0.05 and CI excludes 0')
 }
}
es<-as.data.table(expected('replication_summary'))
stopifnot(nrow(es)==2L)
for(i in seq_len(nrow(es))){
 q<-as.data.table(summ)[contrast==es$contrast[i]]
 stopifnot(nrow(q)==1L)
 for(nm in c('predefined_targets','same_direction','same_direction_supported','supported_lineages')){
  check('C_GSE290679_SUMMARY',es$contrast[i],nm,q[[nm]],es[[nm]][i],0,
        inputs=c('gse290679_pseudobulk_counts','gse290679_direction_registry'),
        estimand='five frozen metabolic/stress Hallmarks x five PBMC lineages')
 }
}
cat('GSE290679_REPLICATION',nrow(res),'rows; 14/25 healthy and 0/25 ICU expected\n');flush.console()
