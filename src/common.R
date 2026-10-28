options(stringsAsFactors=FALSE,warn=1)
# Best-effort locale: do not emit a warning on hosts where this locale is absent.
invisible(suppressWarnings(Sys.setlocale('LC_ALL','English_United States.utf8')))
suppressPackageStartupMessages({library(data.table);library(metafor);library(jsonlite)})
argv<-commandArgs(TRUE)
if(length(argv)<2)stop('Usage: Rscript run release.R CONFIG.json core|benchmarks|pbmc_score|pbmc|scoring_portable|cellular|cross_task|clinical_task|cts_reconstruction')
cfg<-fromJSON(argv[1],simplifyVector=TRUE)
if(length(cfg$r_libraries)) .libPaths(c(cfg$r_libraries,.libPaths()))
code_root<-normalizePath(Sys.getenv('SEPSIS_CODE_ROOT'),winslash='/',mustWork=TRUE)
OUT<-cfg$output;dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
stage<-argv[2];ledger<-list();li<-0L
expected<-function(k)fread(file.path(code_root,cfg$checks[[k]]))
readinput<-function(k){
 p<-cfg$inputs[[k]]
 if(is.null(p)||!file.exists(p))stop('Missing explicit input ',k,': ',p)
 # data.table::fread() only decompresses .gz when R.utils is installed. Fall back
 # to a base-R gzfile() connection so gzip inputs never become a hidden dependency.
 if(grepl('\\.gz$',p,ignore.case=TRUE)&&!requireNamespace('R.utils',quietly=TRUE)){
  con<-gzfile(p,'rt');on.exit(close(con),add=TRUE)
  return(fread(text=readLines(con,warn=FALSE)))
 }
 fread(p)
}
members<-expected('score16_members');pos<-members[direction=='positive',pathway];neg<-members[direction=='negative',pathway]
stopifnot(length(pos)==14,length(neg)==2,!anyDuplicated(c(pos,neg)),!('HALLMARK_E2F_TARGETS'%in%c(pos,neg)))
discovery<-cfg$discovery;external<-cfg$external
stopifnot(identical(external,c('GSE63042','GSE95233','VANISH')),length(discovery)==4)
check<-function(id,cohort,field,value,frozen,tol=cfg$tolerance,inputs=cohort,n=NA_integer_,deaths=NA_integer_,estimand=''){
 li<<-li+1L
 ok<-length(value)==1&&length(frozen)==1&&is.finite(value)&&is.finite(frozen)&&abs(value-frozen)<=tol
 ledger[[li]]<<-data.table(result_id=id,cohort=cohort,n=n,deaths=deaths,input_keys=paste(inputs,collapse=';'),script=paste0('src/',gsub('_',' ',stage,fixed=TRUE),'.R'),estimand=estimand,field=field,reproduced_value=if(length(value)==1)value else NA_real_,frozen_value=if(length(frozen)==1)frozen else NA_real_,tolerance=tol,status=if(ok)'PASS' else 'FAIL',reason=if(ok)'independent computation from listed input' else 'numeric or structural mismatch')
}
finish<-function(){
 x<-rbindlist(ledger,fill=TRUE);fwrite(x,file.path(OUT,paste0(stage,'_comparison.tsv')),sep='\t',na='NA')
 capture.output(sessionInfo(),file=file.path(OUT,paste0(stage,'_sessionInfo.txt')))
 cat(stage,'CHECKS',nrow(x),'FAIL',sum(x$status!='PASS'),'\n');flush.console()
 if(any(x$status!='PASS'))stop('Reproduction failed; frozen outputs unchanged; see comparison ledger')
}
raw_score16<-function(d)rowMeans(as.matrix(d[,..pos]))-rowMeans(as.matrix(d[,..neg]))
zscore<-function(x){stopifnot(all(is.finite(x)),sd(x)>0);as.numeric(scale(x))}
auc<-function(s,y){stopifnot(all(y%in%0:1));n1<-sum(y==1);n0<-sum(y==0);stopifnot(n1>0,n0>0);(sum(rank(s,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)}
fit_score<-function(d,z,extra=character(),critical=qnorm(.975)){
 x<-as.data.frame(d);x$score<-z;fm<-reformulate(c('score',extra),'event_death');vars<-all.vars(fm);x<-x[complete.cases(x[,vars,drop=FALSE]),,drop=FALSE]
 warns<-character();f<-withCallingHandlers(glm(fm,data=x,family=binomial()),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 sm<-coef(summary(f))['score',];b<-unname(sm[1]);se<-unname(sm[2]);stopifnot(f$converged,f$rank==length(coef(f)),is.finite(b),is.finite(se))
 list(fit=f,data=x,row=data.table(n=nrow(x),deaths=sum(x$event_death),logOR=b,SE=se,OR=exp(b),CI_low=exp(b-critical*se),CI_high=exp(b+critical*se),P=unname(sm[4]),AUC=auc(x$score,x$event_death),model_AUC=auc(fitted(f),x$event_death),converged=f$converged,warnings=paste(unique(warns),collapse=';')))
}
metafit<-function(d,test='knha'){
 f<-rma.uni(yi=d$logOR,sei=d$SE,method='REML',test=test);pr<-predict(f,transf=exp)
 data.table(k=nrow(d),n=sum(d$n),deaths=sum(d$deaths),pooled_logOR=as.numeric(f$b),pooled_SE=f$se,OR=exp(as.numeric(f$b)),CI_low=exp(f$ci.lb),CI_high=exp(f$ci.ub),P=f$pval,I2=f$I2,tau2=f$tau2,Q=f$QE,PI_low=pr$pi.lb,PI_high=pr$pi.ub)
}
