# Recompute maturity variance and projected-composition associations from frozen intermediates.
o<-readRDS(cfg$inputs$mixture_object);ms<-as.data.table(o$neutrophil_maturity_scores);ms[,Score16:=raw_score16(ms)]
ee<-expected('maturity');vr<-list()
# The frozen maturity FDR family includes 16 pathways AND the composite (17).
# Preserve that source rule; T7 displays only the 16 pathway rows.
for(p in c(pos,neg,'Score16')){
 a<-ms[reconstruction=='Observed'][order(donor_id)];b<-ms[reconstruction=='Maturity_standardized'][order(donor_id)];stopifnot(identical(a$donor_id,b$donor_id),nrow(a)==17)
 ro<-suppressWarnings(cor.test(a$actual_immature_progenitor_weight,a[[p]],method='spearman',exact=FALSE));rs<-suppressWarnings(cor.test(b$actual_immature_progenitor_weight,b[[p]],method='spearman',exact=FALSE))
 r<-data.table(metric=p,n=nrow(a),observed_variance=var(a[[p]]),standardized_variance=var(b[[p]]),variance_reduction=1-var(b[[p]])/var(a[[p]]),rho_observed=unname(ro$estimate),P_observed=ro$p.value,rho_standardized=unname(rs$estimate),P_standardized=rs$p.value);vr[[p]]<-r
}
v<-rbindlist(vr);v[,FDR_observed:=p.adjust(P_observed,'BH')];v[,FDR_standardized:=p.adjust(P_standardized,'BH')]
for(p in setdiff(v$metric,'Score16')){r<-v[metric==p];e<-ee[metric==p];for(nm in setdiff(names(v),'metric'))check('C_NEUTROPHIL_MATURITY',p,nm,r[[nm]],e[[nm]],inputs='mixture_object',n=r$n,estimand='paired 17-donor maturity standardization; original 17-metric BH family')}
fwrite(v,file.path(OUT,'maturity_variance_recomputed.csv'))
pr<-readinput('deconv_estimates')[cohort=='GSE185263'];stopifnot(!anyDuplicated(pr[,.(sample_id,method,sensitivity)]))
# Two frozen scenarios: the primary reference matrix and the reference built without
# the 16 score genes. Map the second one onto the package naming.
scen<-unique(pr$sensitivity);stopifnot(length(scen)==2L,'primary'%in%scen)
pr[,sensitivity:=ifelse(sensitivity=='primary','primary','exclude_Score16')]
dt<-readinput('GSE185263');corr<-list();models<-list();j<-0L
for(mm in unique(pr$method))for(ss in unique(pr$sensitivity))for(target in c('immature_within_neut','overall_immature')){
 d<-pr[method==mm&sensitivity==ss];if(!nrow(d))next
 ii<-match(d$sample_id,dt$sample_id);stopifnot(!anyNA(ii),identical(as.integer(d$event),as.integer(dt$event_death[ii])))
 raw<-raw_score16(dt)[ii];keep<-is.finite(d[[target]])&is.finite(raw);dd<-d[keep];dd[,event_death:=event];az<-zscore(raw[keep]);cz<-zscore(dd[[target]])
 ct<-suppressWarnings(cor.test(raw[keep],dd[[target]],method='spearman',exact=FALSE));r<-data.table(cohort='GSE185263',method=mm,sensitivity=ss,target=target,n=nrow(dd),rho=unname(ct$estimate),P=ct$p.value,VIF=1/(1-cor(az,cz)^2));corr[[paste(mm,ss,target)]]<-r
 e<-expected('deconv_corr')[method==mm&sensitivity==ss&target==r$target]
 for(nm in c('n','rho','P','VIF'))check('C_DECONV_CORRELATION',paste(mm,ss,target),nm,r[[nm]],e[[nm]],inputs=c('deconv_estimates','GSE185263'),n=nrow(dd),estimand='frozen inferred fractions; no new MuSiC/DWLS projection')
 for(joint in c(FALSE,TRUE)){
  dd[,composition_z:=cz];f<-fit_score(dd,az,if(joint)'composition_z' else character(),critical=1.96);r<-f$row
  lab<-paste(ss,mm,target,if(joint)'joint' else 'Score16_only',sep=';');e<-expected('deconv_models')[analysis==lab]
  for(nm in c('n','deaths','logOR','SE','OR','CI_low','CI_high','P','AUC','model_AUC'))check('C_DECONV_MODEL',lab,nm,r[[nm]],e[[nm]],inputs=c('deconv_estimates','GSE185263'),n=r$n,deaths=r$deaths,estimand='axis association conditioned on frozen inferred composition')
  j<-j+1L;models[[j]]<-cbind(data.table(analysis=lab),r)
 }
}
fwrite(rbindlist(corr),file.path(OUT,'deconv_correlation_recomputed.csv'));fwrite(rbindlist(models),file.path(OUT,'deconv_models_recomputed.csv'))
