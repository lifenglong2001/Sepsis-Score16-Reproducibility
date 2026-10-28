suppressPackageStartupMessages({library(GSVA);library(BiocParallel)})
stopifnot(as.character(packageVersion('GSVA'))=='2.4.9')
gmt<-strsplit(readLines(cfg$inputs$gmt),'\t');sets<-setNames(lapply(gmt,function(x)unique(x[-c(1,2)])),vapply(gmt,'[',character(1),1));stopifnot(length(sets)==50)
ranksall<-list();ssall<-list();coverage<-list()
for(u in c(discovery,external)){
 tab<-readinput(paste0('expression_',u));gn<-as.character(tab[[1]]);m<-as.matrix(tab[,-1]);storage.mode(m)<-'double';rownames(m)<-gn
 stopifnot(!anyDuplicated(gn),!anyDuplicated(colnames(m)),!anyNA(m),all(is.finite(m)))
 d<-readinput(u);key<-d$sample_id
 if(!all(key%in%colnames(m))){if('sample_label'%in%names(d)&&all(d$sample_label%in%colnames(m)))key<-d$sample_label else stop('Unresolved exact expression sample mapping for ',u)}
 ii<-match(key,colnames(m));stopifnot(!anyNA(ii))
 # The original discovery ssGSEA script used binary RDS matrices. The portable
 # ranking script used text exports. Their tied ranks are not interchangeable.
 # Keep BOTH source choices exactly as documented by those source scripts.
 score_m<-if(u%in%discovery)as.matrix(readRDS(cfg$inputs[[paste0('ssgsea_expression_',u)]])) else m
 score_key<-if(all(d$sample_id%in%colnames(score_m)))d$sample_id else key
 score_ii<-match(score_key,colnames(score_m));stopifnot(!anyNA(score_ii))
 sc<-gsva(ssgseaParam(score_m,sets,minSize=10,maxSize=Inf,alpha=.25,normalize=TRUE),verbose=FALSE,BPPARAM=SerialParam())
 hs<-names(sets);sc<-sc[hs,score_ii,drop=FALSE]
 err<-max(abs(t(sc)-as.matrix(d[,..hs])))
 score_input<-if(u%in%discovery)paste0('ssgsea_expression_',u) else paste0('expression_',u)
 check('C_SSGSEA50',u,'max_score_difference',err,0,cfg$score_tolerance,inputs=c(score_input,u,'gmt'),n=nrow(d),deaths=sum(d$event_death),estimand='GSVA 2.4.9 ssGSEA from the original processed expression input')
 ss<-as.data.table(t(sc));zz<-zscore(raw_score16(ss));check('C_SSGSEA_SCORE16',u,'max_z_difference',max(abs(zz-zscore(raw_score16(d)))),0,cfg$score_tolerance,inputs=c(score_input,u,'gmt'))
 ssall[[u]]<-cbind(data.table(cohort=u,sample_id=d$sample_id),ss)
 N<-nrow(m);ix<-lapply(sets[c(pos,neg)],function(g)which(gn%in%g));stopifnot(all(lengths(ix)>=10),all(lengths(ix)<N))
 rankmat<-apply(m,2,rank,ties.method='average')
 sc_rank<-t(vapply(ix,function(j){k<-length(j);(colMeans(rankmat[j,,drop=FALSE])-(k+1)/2)/(N-k)-.5},numeric(ncol(m))))
 rownames(sc_rank)<-c(pos,neg);colnames(sc_rank)<-colnames(m)
 rs<-as.data.table(t(sc_rank[,ii,drop=FALSE]));rs[,sample_id:=d$sample_id];ranksall[[u]]<-rs
 single<-rank(m[,1],ties.method='average');one<-vapply(ix,function(j)(mean(single[j])-(length(j)+1)/2)/(N-length(j))-.5,numeric(1))
 check('C_PORTABLE_INVARIANCE',u,'max_difference',max(abs(one-sc_rank[,1])),0,1e-12,inputs=c(paste0('expression_',u),'gmt'))
 coverage[[u]]<-data.table(cohort=u,pathway=names(ix),universe=N,matched=lengths(ix),original=lengths(sets[names(ix)]))
 cat('EXPRESSION_SCORED',u,nrow(m),ncol(m),'max_original_difference',err,'\n');flush.console();rm(tab,m,score_m,rankmat,sc,sc_rank);gc()
}
refs<-rbindlist(lapply(discovery,function(u){v<-raw_score16(ranksall[[u]]);data.table(cohort=u,n=length(v),mean=mean(v),variance=var(v))}))
mu<-mean(refs$mean);sdev<-sqrt(mean(refs$variance+refs$mean^2)-mu^2)
fwrite(refs,file.path(OUT,'portable_reference_cohorts.csv'));fwrite(data.table(mean=mu,SD=sdev,external_used=FALSE),file.path(OUT,'portable_reference_lock.csv'))
saveRDS(list(ranks=ranksall,mean=mu,SD=sdev),file.path(OUT,'portable_recomputed.rds'));fwrite(rbindlist(ssall),file.path(OUT,'ssgsea50_recomputed.csv'));fwrite(rbindlist(coverage),file.path(OUT,'scoring_coverage.csv'))
bootauc<-function(y,p0,p1,seed){set.seed(seed);ss<-split(seq_along(y),y);b<-t(replicate(1000,{jj<-unlist(lapply(ss,function(i)sample(i,length(i),TRUE)));c(auc(p0[jj],y[jj]),auc(p1[jj],y[jj]))}));ci<-function(x)quantile(x,c(.025,.975),type=6,names=FALSE);c(AUC_CI_low=ci(b[,2])[1],AUC_CI_high=ci(b[,2])[2],delta_AUC=auc(p1,y)-auc(p0,y),delta_AUC_low=ci(b[,2]-b[,1])[1],delta_AUC_high=ci(b[,2]-b[,1])[2])}
por<-list()
for(u in external){
 d<-readinput(u);r<-ranksall[[u]];ix<-match(d$sample_id,r$sample_id);score<-(raw_score16(r)[ix]-mu)/sdev;z<-zscore(raw_score16(d));q<-fit_score(d,score,critical=1.96)$row
 seed<-2026101750+match(u,external);bb<-bootauc(d$event_death,z,score,seed);q<-cbind(data.table(cohort=u),q,as.data.table(as.list(bb)));por[[u]]<-q;e<-expected('portable')[cohort==u]
 for(nm in c('n','deaths','logOR','SE','OR','CI_low','CI_high','P','AUC','AUC_CI_low','AUC_CI_high','delta_AUC','delta_AUC_low','delta_AUC_high'))check('C_PORTABLE_MODEL',u,nm,q[[nm]],e[[nm]],inputs=c(paste0('expression_',c(discovery,u)),u,'gmt'),n=q$n,deaths=q$deaths,estimand='rank score per discovery-only reference SD; direction never reselected')
}
fwrite(rbindlist(por),file.path(OUT,'portable_refitted.csv'))
