datasets<-list();effects<-list();all50<-list();collapsed<-list();patient_scores<-list()
for(u in c(discovery,external)){
 d<-readinput(u);h<-grep('^HALLMARK_',names(d),value=TRUE)
 check('A_INPUT',u,'n',nrow(d),cfg$expected_n[[u]],0,n=nrow(d),deaths=sum(d$event_death))
 check('A_INPUT',u,'deaths',sum(d$event_death),cfg$expected_deaths[[u]],0)
 check('A_INPUT',u,'hallmarks',length(h),50,0)
 check('A_INPUT',u,'duplicate_patient',anyDuplicated(d$patient_id),0,0)
 check('A_INPUT',u,'duplicate_sample',anyDuplicated(d$sample_id),0,0)
 stopifnot(all(d$event_death%in%0:1),all(is.finite(as.matrix(d[,..h]))))
 datasets[[u]]<-d;raw<-raw_score16(d);z<-zscore(raw)
 fit<-fit_score(d,z);rr<-cbind(data.table(cohort=u),fit$row);effects[[u]]<-rr
 e<-if(u%in%discovery)expected('discovery_score16')[cohort==u] else expected('external_score16')[cohort==u]
 for(nm in c('logOR','SE','OR','CI_low','CI_high','P','AUC'))check('B_SCORE16',u,nm,rr[[nm]],e[[nm]],if(nm=='AUC')cfg$auc_tolerance else cfg$tolerance,n=rr$n,deaths=rr$deaths,estimand='death=1; OR per within-cohort 1 SD; fixed direction')
 patient_scores[[u]]<-data.table(cohort=u,sample_id=d$sample_id,patient_id=d$patient_id,event_death=d$event_death,Score16_raw=raw,Score16_z=z)
 if(u%in%discovery){
  for(p in sort(h)){
   a<-fit_score(d,zscore(d[[p]]))$row
   q<-expected('hallmark_cohort')[cohort==u&feature==p]
   for(nm in c('logOR','SE','OR','CI_low','CI_high','P')){
    en<-switch(nm,logOR='logOR_per_1SD',SE='SE',OR='OR_per_1SD',CI_low='OR_CI_low',CI_high='OR_CI_high',P='P_value')
    check(p,u,nm,a[[nm]],q[[en]],n=a$n,deaths=a$deaths,estimand='univariable mortality; per 1 SD pathway')
   }
   all50[[paste(u,p)]]<-cbind(data.table(cohort=u,pathway=p),a)
  }
  cyc<-c('HALLMARK_G2M_CHECKPOINT','HALLMARK_MITOTIC_SPINDLE');other<-setdiff(pos,cyc)
  rawc<-(rowSums(as.matrix(d[,..other]))+rowMeans(as.matrix(d[,..cyc])))/(length(other)+1)-rowMeans(as.matrix(d[,..neg]))
  a<-fit_score(d,zscore(rawc),critical=1.96)$row;collapsed[[u]]<-cbind(data.table(cohort=u),a)
  e<-expected('collapsed')[cohort==u]
  for(nm in c('logOR','SE','OR','CI_low','CI_high','P','AUC'))check('B_CELL_CYCLE_COLLAPSED',u,nm,a[[nm]],e[[nm]],n=a$n,deaths=a$deaths)
 }
 cat('REFIT',u,'n',nrow(d),'deaths',sum(d$event_death),'\n');flush.console()
}
a<-rbindlist(effects);t<-rbindlist(all50);ct<-rbindlist(collapsed)
fwrite(a,file.path(OUT,'score16_cohort_effects.csv'));fwrite(t,file.path(OUT,'hallmark50_cohort_effects.csv'));fwrite(ct,file.path(OUT,'cellcycle_collapsed.csv'));fwrite(rbindlist(patient_scores),file.path(OUT,'score16_patient_scores.csv'))
full<-rbindlist(lapply(sort(unique(t$pathway)),function(p){
  dd<-t[pathway==p];f<-rma.uni(yi=dd$logOR,sei=dd$SE,method='REML',test='z');h<-rma.uni(yi=dd$logOR,sei=dd$SE,method='REML',test='knha');ad<-rma.uni(yi=dd$logOR,sei=dd$SE,method='REML',test='adhoc')
  data.table(pathway=p,pooled_logOR=as.numeric(f$b),SE=f$se,P=f$pval,P_HKSJ=h$pval,P_modified_KH=ad$pval,OR=exp(as.numeric(f$b)),CI_low=exp(f$ci.lb),CI_high=exp(f$ci.ub),I2=f$I2,tau2=f$tau2,strict_direction=all(sign(dd$logOR)==sign(as.numeric(f$b))))
 }));full[,`:=`(FDR=p.adjust(P,'BH'),FDR_HKSJ=p.adjust(P_HKSJ,'BH'),FDR_modified_KH=p.adjust(P_modified_KH,'BH'))];full[,selected:=FDR<.05&strict_direction]
ee<-expected('hallmark_meta')
for(p in full$pathway){r<-full[pathway==p];e<-ee[feature==p];for(nm in c('pooled_logOR','SE','P','OR','CI_low','CI_high','I2','tau2','FDR')){
 en<-switch(nm,SE='pooled_SE',P='P_value',OR='pooled_OR',I2='I2_percent',FDR='FDR_BH',nm)
 check('B_HALLMARK_META',p,nm,r[[nm]],e[[en]],inputs=discovery,estimand='REML-z; BH over exactly 50 pathways')
}}
check('B_SELECTION','discovery','z_BH_count',sum(full$FDR<.05),19,0,inputs=discovery)
check('B_SELECTION','discovery','concordant_count',sum(full$selected),16,0,inputs=discovery)
check('B_SELECTION','discovery','member_set_difference',length(setdiff(full[selected==TRUE,pathway],members$pathway))+length(setdiff(members$pathway,full[selected==TRUE,pathway])),0,0,inputs=discovery)
check('B_SELECTION','discovery','direction_mismatch',sum(sign(full[match(members$pathway,pathway),pooled_logOR])!=ifelse(members$direction=='positive',1,-1)),0,0,inputs=discovery)
fwrite(full,file.path(OUT,'hallmark50_primary_meta.csv'))
metas<-list()
for(g in c('discovery','external')){
 dd<-a[cohort%in%if(g=='discovery')discovery else external]
 q<-metafit(dd);metas[[g]]<-cbind(data.table(group=g,analysis='primary',omitted='none'),q)
 if(g=='discovery'){
  e<-expected('adult_meta')[axis=='Score16'&group=='discovery'&analysis=='primary_discovery_four_units']
  for(nm in names(q))check('B_DISCOVERY_META','discovery',nm,q[[nm]],e[[nm]],inputs=discovery,n=q$n,deaths=q$deaths,estimand='four discovery cohorts; REML-HKSJ')
 }
 if(g=='discovery'){
  # Fixed original Score16 scores: omit an effect estimate, never a pathway.
  for(o in discovery){
   q<-metafit(dd[cohort!=o]);e<-expected('discovery_loco')[omitted==o]
   stopifnot(nrow(e)==1)
   metas[[paste(g,o)]]<-cbind(data.table(group=g,analysis='LOCO',omitted=o),q)
   for(nm in names(q))check('B_DISCOVERY_FIXED_SCORE_LOCO',o,nm,q[[nm]],e[[nm]],inputs=setdiff(discovery,o),n=q$n,deaths=q$deaths,estimand='Original fixed Score16; leave-one-discovery-cohort-out REML-HKSJ')
  }
 }
 if(g=='external'){
  e<-expected('external_meta')[analysis=='primary_three_units']
  for(nm in names(q))check('B_EXTERNAL_META','external',nm,q[[nm]],e[[nm]],inputs=external,n=q$n,deaths=q$deaths,estimand='REML-HKSJ; default prediction interval; three external units')
  for(o in external){
   q<-metafit(dd[cohort!=o]);e<-expected('external_meta')[analysis=='leave_one_study_out'&omitted==o];metas[[paste(g,o)]]<-cbind(data.table(group=g,analysis='LOSO',omitted=o),q)
   for(nm in names(q))check('B_EXTERNAL_LOSO',o,nm,q[[nm]],e[[nm]],inputs=setdiff(external,o),n=q$n,deaths=q$deaths)
  }
  for(test in c('z','adhoc')){q<-metafit(dd,test);e<-expected('external_meta')[analysis==paste0('alternate_REML_',test)];metas[[test]]<-cbind(data.table(group=g,analysis=test,omitted='none'),q);for(nm in names(q))check('B_EXTERNAL_ALTERNATE',test,nm,q[[nm]],e[[nm]],inputs=external)}
 }
}
cm<-metafit(ct);ce<-expected('collapsed_meta');for(nm in names(cm))check('B_COLLAPSED_META','discovery',nm,cm[[nm]],ce[[nm]],inputs=discovery)
fwrite(cm,file.path(OUT,'cellcycle_collapsed_meta.csv'));fwrite(rbindlist(metas),file.path(OUT,'score16_meta_LOSO.csv'))
saveRDS(datasets,file.path(OUT,'patient_inputs_checked.rds'))
