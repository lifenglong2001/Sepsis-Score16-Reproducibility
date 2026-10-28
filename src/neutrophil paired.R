# GSE216009 paired immature/progenitor vs mature/activated neutrophil pseudobulk (T7 paired subset).
# Donor is the unit; 17 acute-sepsis donors with >=50 cells in both maturity states.
#
# Scope note (deliberate and documented):
#  * This stage reproduces ONLY the 'neutrophil paired states' /
#    'Aggregated immature/progenitor vs mature/activated' block (51 metrics).
#  * The frozen T7 file also holds 714 unpaired cell-subtype rows and 48 paired
#    subtype-vs-mature rows whose cell-state pseudobulk inputs are NOT part of
#    this release; those remain frozen targets only.
#  * The frozen BH families (fdr_global, fdr_within_metric) span all 381 tested
#    rows of the 813-row table, so they are NOT re-derived here. Only effect
#    sizes, medians and p values are recomputed.
#  * The 'Score16' metric here is the 16-member composite (14 positive minus 2
#    negative Hallmarks).
d<-as.data.table(readinput('gse216009_neutrophil_maturity_scores'))
need<-c('donor_id','phase','clinical_group','cell_type','cell_count','cell_annotation_level')
stopifnot(all(need%in%names(d)))
stopifnot(unique(d$cell_annotation_level)=='neutrophil_maturity')
stopifnot(setequal(unique(d$cell_type),c('Immature/progenitor','Mature/activated')))

sep<-d[clinical_group=='Acute Sepsis' & phase=='Acute']
both<-sep[, .N, by=donor_id][N==2L, donor_id]
stopifnot(length(both)==17L)
sep<-sep[donor_id%in%both]
imm<-sep[cell_type=='Immature/progenitor'][order(donor_id)]
mat<-sep[cell_type=='Mature/activated'][order(donor_id)]
stopifnot(nrow(imm)==17L,nrow(mat)==17L,identical(imm$donor_id,mat$donor_id))
stopifnot(all(imm$cell_count>=50),all(mat$cell_count>=50))

metrics<-c(grep('^HALLMARK_',names(d),value=TRUE),'Score16')
stopifnot(length(metrics)==51L)
res<-rbindlist(lapply(metrics,function(m){
 if(m=='Score16'){
  x<-rowMeans(as.matrix(imm[,..pos]))-rowMeans(as.matrix(imm[,..neg]))
  y<-rowMeans(as.matrix(mat[,..pos]))-rowMeans(as.matrix(mat[,..neg]))
 }else{
  x<-imm[[m]];y<-mat[[m]]
 }
 w<-suppressWarnings(wilcox.test(x,y,paired=TRUE,exact=FALSE))
 data.table(frozen_metric=m,
            metric=m,n_sepsis=length(x),n_control=length(y),
            immature_median=median(x),mature_median=median(y),
            median_paired_difference=median(x-y),P=w$p.value)
}))
fwrite(res,file.path(OUT,'neutrophil_paired_immature_vs_mature.csv'))

e<-as.data.table(expected('neutrophil_paired'))
e<-e[annotation_level=='neutrophil paired states' &
     cell_annotation=='Aggregated immature/progenitor vs mature/activated']
stopifnot(nrow(e)==51L,all(e$minimum_rule_met))
setkey(e,metric);setkey(res,frozen_metric)
stopifnot(identical(e$metric,res$frozen_metric))
for(nm in c('n_sepsis','n_control')){
 for(i in seq_len(nrow(res))){
  check('C_NEUTROPHIL_PAIRED',res$frozen_metric[i],nm,res[[nm]][i],e[[nm]][i],0,
        inputs='gse216009_neutrophil_maturity_scores',n=res$n_sepsis[i],
        estimand='17 paired acute-sepsis donors; immature/progenitor vs mature/activated')
 }
}
for(nm in c('immature_median','mature_median','median_paired_difference','P')){
 en<-switch(nm,immature_median='sepsis_median',mature_median='control_median',
            median_paired_difference='median_difference_sepsis_minus_control',P='p_value')
 for(i in seq_len(nrow(res))){
  check('C_NEUTROPHIL_PAIRED',res$frozen_metric[i],nm,res[[nm]][i],e[[en]][i],cfg$tolerance,
        inputs='gse216009_neutrophil_maturity_scores',n=res$n_sepsis[i],
        estimand='paired Wilcoxon signed-rank exact=FALSE; median of within-donor differences; frozen FDR family not re-derived')
 }
}
cat('NEUTROPHIL_PAIRED',nrow(res),'metrics reproduced (51/99 paired rows); frozen FDR family spans 381 rows and is not re-derived\n')
flush.console()
