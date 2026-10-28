# Refit the frozen benchmark formulas and the original paired, stratified bootstrap.
# CTS/SRS/Hi-DEF/MARS labels are frozen inputs; no label classifier is trained here.
bench<-list();draws<-list();ixout<-0L
for(u in c('VANISH','GSE65682')){
 d<-readinput(if(u=='VANISH')'vanish_labels' else 'mars_labels');ss<-readinput(u)
 ii<-match(d$sample_id,ss$sample_id);stopifnot(!anyNA(ii),!anyDuplicated(d$sample_id),identical(as.integer(d$event_death),as.integer(ss$event_death[ii])))
 if(u=='VANISH'){
  stopifnot(identical(as.character(d$CTS),as.character(ss$CTS[ii])))
  d[,`:=`(CTS=factor(CTS),SRS=factor(SRS),sex=factor(sex))]
  models<-list(CTS='CTS',SRS='SRS',HiDEF=c('HiDEF_myeloid','HiDEF_lymphoid'))
 } else {d[,`:=`(MARS=factor(MARS),sex=factor(sex))];models<-list(MARS='MARS')}
 z<-zscore(raw_score16(ss))[ii]
 for(bm in names(models))for(clin in c(FALSE,TRUE)){
  ex<-c(models[[bm]],if(clin)c('age','sex',if(u=='VANISH')'APACHE_II'))
  f<-fit_score(d,z,ex,critical=1.96);x<-f$data
  f0<-glm(reformulate(ex,'event_death'),data=x,family=binomial())
  lr<-anova(f0,f$fit,test='LRT');row<-f$row
  seed<-2026101700+match(bm,c('CTS','SRS','HiDEF','MARS'))*10+as.integer(clin)
  B<-cfg$bootstrap_B;stopifnot(B==2000)
  set.seed(seed);strata<-split(seq_len(nrow(x)),x$event_death);draw<-matrix(NA_real_,B,2);fail<-0L
  for(b in seq_len(B)){
   jj<-unlist(lapply(strata,function(kk)sample(kk,length(kk),TRUE)));xx<-x[jj,,drop=FALSE]
   f0b<-tryCatch(suppressWarnings(glm(reformulate(ex,'event_death'),data=xx,family=binomial())),error=function(e)NULL)
   f1b<-tryCatch(suppressWarnings(glm(reformulate(c('score',ex),'event_death'),data=xx,family=binomial())),error=function(e)NULL)
   if(!is.null(f0b)&&!is.null(f1b)&&isTRUE(f0b$converged)&&isTRUE(f1b$converged))draw[b,]<-c(auc(fitted(f0b),xx$event_death),auc(fitted(f1b),xx$event_death)) else fail<-fail+1L
  }
  ci<-function(v)quantile(v,c(.025,.975),type=6,na.rm=TRUE,names=FALSE)
  lab<-paste(bm,if(clin)'clinical' else 'unadjusted',sep='_')
  row[,`:=`(cohort=u,benchmark=bm,analysis=lab,clinical_adjustment=clin,LRT_1df=lr$Deviance[2],LRT_P=lr$`Pr(>Chi)`[2],AUC_base=auc(fitted(f0),x$event_death),AUC_full=auc(fitted(f$fit),x$event_death),AUC_base_low=ci(draw[,1])[1],AUC_base_high=ci(draw[,1])[2],AUC_full_low=ci(draw[,2])[1],AUC_full_high=ci(draw[,2])[2],delta_AUC=auc(fitted(f$fit),x$event_death)-auc(fitted(f0),x$event_death),delta_AUC_low=ci(draw[,2]-draw[,1])[1],delta_AUC_high=ci(draw[,2]-draw[,1])[2],bootstrap_failures=fail,seed=seed,B=B)]
  e<-expected(if(u=='VANISH')'vanish_benchmarks' else 'mars_benchmarks')[analysis==lab]
  for(nm in c('n','deaths','logOR','SE','OR','CI_low','CI_high','P','AUC','model_AUC','LRT_1df','LRT_P','AUC_base','AUC_full','AUC_base_low','AUC_base_high','AUC_full_low','AUC_full_high','delta_AUC','delta_AUC_low','delta_AUC_high','B','seed','bootstrap_failures'))check(paste0('B_',lab),u,nm,row[[nm]],e[[nm]],inputs=c(u,if(u=='VANISH')'vanish_labels' else 'mars_labels'),n=row$n,deaths=row$deaths,estimand='apparent AUC; 2000 paired stratified refits; percentile type 6; LRT separately recorded')
  ixout<-ixout+1L;bench[[ixout]]<-row;draws[[paste(u,lab)]]<-draw
  cat('BENCHMARK_REFIT',u,lab,'n',nrow(x),'delta_AUC',row$delta_AUC,'failures',fail,'\n');flush.console()
 }
}
rr<-rbindlist(bench);rr[,BH_3_SRS_HiDEF_MARS:=NA_real_]
rr[clinical_adjustment==FALSE&benchmark%in%c('SRS','HiDEF','MARS'),BH_3_SRS_HiDEF_MARS:=p.adjust(LRT_P,'BH')]
for(i in which(!is.na(rr$BH_3_SRS_HiDEF_MARS))){e<-expected(if(rr$cohort[i]=='VANISH')'vanish_benchmarks' else 'mars_benchmarks')[analysis==rr$analysis[i]];check('B_BENCHMARK_BH3',rr$cohort[i],rr$analysis[i],rr$BH_3_SRS_HiDEF_MARS[i],e$BH_3_SRS_HiDEF_MARS,inputs=c('vanish_labels','mars_labels'))}
fwrite(rr,file.path(OUT,'benchmarks_refitted.csv'));saveRDS(draws,file.path(OUT,'benchmarks_bootstrap_draws.rds'))
