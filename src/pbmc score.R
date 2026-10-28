# Recompute ssGSEA from provenance-recorded normalized donor pseudobulk.
# This is intermediate-input reproduction, not raw single-cell reconstruction.
suppressPackageStartupMessages({library(GSVA);library(BiocParallel)})
stopifnot(as.character(packageVersion('GSVA'))=='2.4.9')
o<-readRDS(cfg$inputs$pbmc_object)
gs<-strsplit(readLines(cfg$inputs$gmt),'\t');sets<-setNames(lapply(gs,function(z)unique(z[-(1:2)])),vapply(gs,'[',character(1),1))
metrics<-unique(expected('pbmc')$metric);stopifnot(setequal(metrics,c(pos,neg)))
out<-list();saved<-readinput('pbmc_scores')
for(ct in names(o$normalized)){
 nn<-o$normalized[[ct]]
 sc<-gsva(ssgseaParam(nn$logcpm,sets[metrics],minSize=10,maxSize=Inf,alpha=.25,normalize=TRUE),verbose=FALSE,BPPARAM=SerialParam())
 tab<-cbind(nn$metadata,as.data.table(t(sc)));tab[,Score16:=raw_score16(tab)];out[[ct]]<-tab
 e<-saved[cell_type==ct];ii<-match(tab$donor_id,e$donor_id);stopifnot(!anyNA(ii),!anyDuplicated(tab$donor_id))
 for(p in metrics)check('C_PBMC_SSGSEA',ct,p,max(abs(tab[[p]]-e[[p]][ii])),0,cfg$score_tolerance,inputs=c('pbmc_object','gmt','pbmc_scores'),n=nrow(tab),estimand='GSVA 2.4.9 ssGSEA from frozen donor logCPM; not raw single cells')
 cat('PBMC_SSGSEA_REPRODUCED',ct,nrow(tab),'\n');flush.console()
}
fwrite(rbindlist(out),file.path(OUT,'pbmc_scores_from_logcpm.csv'))
