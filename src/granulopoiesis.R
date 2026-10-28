# Five-gene granulopoiesis module (AZU1 / CTSG / ELANE / MPO / PRTN3) and its
# relation to the fixed Score16 (Figure 4c-e).
#
# The within-sample rank score, the bootstrap counts (2000 for the cohort
# correlations, 1000 for the subtype pairwise Cliff deltas) and the percentile
# interval rule follow the archived screen.
#
# Scope note: this stage covers the gene-module analyses only. It does not fit a
# mortality model; the MARS mortality adjustments belong to `benchmarks`.
genes<-c('AZU1','CTSG','ELANE','MPO','PRTN3')
STAGE_SEED<-2026101961
set.seed(STAGE_SEED)
ci<-function(x)as.numeric(quantile(x,c(.025,.975),type=6,na.rm=TRUE,names=FALSE))
cliff<-function(a,b){n<-length(a);m<-length(b);2*(sum(rank(c(a,b))[seq_len(n)])-n*(n+1)/2)/(n*m)-1}
score<-function(m,g){
 stopifnot(!anyDuplicated(rownames(m)),all(g%in%rownames(m)),nrow(m)>=5000,all(is.finite(m)))
 ix<-match(g,rownames(m))
 vapply(seq_len(ncol(m)),function(j)(sum(rank(m[,j],ties.method='average')[ix])-length(g)*.5)/(length(g)*nrow(m)),numeric(1))
}

allpatients<-list();res<-list();cover<-list();index<-0L
for(u in c(discovery,external)){
 tab<-readinput(paste0('expression_',u))
 gn<-as.character(tab[[1]]);mat<-as.matrix(tab[,-1]);storage.mode(mat)<-'double';rownames(mat)<-gn
 dd<-copy(readinput(u))
 key<-if(all(dd$sample_id%in%colnames(mat)))dd$sample_id else
      if('sample_label'%in%names(dd)&&all(dd$sample_label%in%colnames(mat)))dd$sample_label else
      stop('ExactIDmatchfail ',u)
 ii<-match(key,colnames(mat));mm<-mat[,ii,drop=FALSE]
 stopifnot(!anyNA(ii),!anyDuplicated(key),ncol(mm)==nrow(dd))
 present<-genes%in%gn
 cover[[u]]<-data.table(cohort=u,gene=genes,present=present,universe=nrow(mm),n=nrow(dd),
                        events=sum(dd$event_death))
 stopifnot(all(present))   # frozen marker set must be complete in every cohort
 primary<-score(mm,genes)
 raw<-rowMeans(as.matrix(dd[,..pos]))-rowMeans(as.matrix(dd[,..neg]));az<-as.numeric(scale(raw))
 dd[,`:=`(cohort=u,Score16_raw=raw,Score16_z=az,granulopoiesis5=primary)]
 x<-dd[['granulopoiesis5']]
 ct<-suppressWarnings(cor.test(x,az,method='spearman',exact=FALSE))
 bt<-replicate(2000,{i<-sample.int(length(x),length(x),TRUE);cor(x[i],az[i],method='spearman')})
 index<-index+1L
 res[[index]]<-data.table(cohort=u,definition='granulopoiesis5',n=length(x),events=sum(dd$event_death),
                          rho=unname(ct$estimate),CI_low=ci(bt)[1],CI_high=ci(bt)[2],
                          P=ct$p.value,universe=nrow(mm),marker_n=5)
 allpatients[[u]]<-dd
 cat('MODULE_SCORED',u,'n',nrow(dd),'rho',res[[index]]$rho,'\n');flush.console()
 rm(tab,mat,mm);gc()
}
rr<-rbindlist(res)
rr[,FDR_BH7:=p.adjust(P,'BH'),by=definition]

# ---- MARS (GSE65682) subtype comparison -------------------------------------
subsummary<-subpairs<-subglobal<-list();si<-pi<-gi<-0L
for(u in intersect('GSE65682',names(allpatients))){
 x<-copy(allpatients[[u]])
 labels<-readinput('mars_joined')
 ii<-match(x$sample_id,labels$sample_id)
 stopifnot(!anyNA(ii),identical(as.integer(x$event_death),as.integer(labels$event_death[ii])))
 x[,subtype:=as.character(labels[['MARS']][ii])]
 for(g in sort(unique(x$subtype))){
  a<-x[subtype==g];si<-si+1L
  subsummary[[si]]<-data.table(cohort=u,label='MARS',subtype=g,definition='granulopoiesis5',n=nrow(a),
                               events=sum(a$event_death),median=median(a[['granulopoiesis5']]),
                               q25=quantile(a[['granulopoiesis5']],.25),q75=quantile(a[['granulopoiesis5']],.75))
 }
 kw<-kruskal.test(x[['granulopoiesis5']]~x$subtype);gi<-gi+1L
 subglobal[[gi]]<-data.table(cohort=u,label='MARS',definition='granulopoiesis5',
                             H=unname(kw$statistic),P=kw$p.value)
 for(g in combn(sort(unique(x$subtype)),2,simplify=FALSE)){
  a<-x[subtype==g[1]][['granulopoiesis5']];b<-x[subtype==g[2]][['granulopoiesis5']]
  bt<-replicate(1000,cliff(sample(a,length(a),TRUE),sample(b,length(b),TRUE)))
  pi<-pi+1L
  subpairs[[pi]]<-data.table(cohort=u,definition='granulopoiesis5',g1=g[1],g2=g[2],
                             n1=length(a),n2=length(b),
                             Cliff_delta=cliff(a,b),CI_low=ci(bt)[1],CI_high=ci(bt)[2],
                             P=wilcox.test(a,b,exact=FALSE)$p.value)
 }
}
sp<-rbindlist(subpairs);sp[,Holm_within_cohort_definition:=p.adjust(P,'holm'),by=.(cohort,definition)]
sg<-rbindlist(subglobal)

# ---- GSE216009 immature vs mature module localization ----------------------
o<-readRDS(cfg$inputs$gse216009_pseudobulk_expression)
ob<-o$neutrophil_maturity;m<-ob$logcpm;meta<-as.data.table(ob$metadata)
stopifnot(identical(as.character(meta$pseudobulk_id),colnames(m)),all(genes%in%rownames(m)))
meta[,granulopoiesis5:=score(m,genes)]
stopifnot('cell_count'%in%names(meta))
eligible<-meta[phase=='Acute'&clinical_group=='Acute Sepsis'&cell_count>=50]
w<-dcast(eligible,donor_id~cell_type,value.var='granulopoiesis5');w<-w[complete.cases(w)]
stopifnot(all(c('Immature/progenitor','Mature/activated')%in%names(w)))
delta<-w[['Immature/progenitor']]-w[['Mature/activated']]
test<-wilcox.test(w[['Immature/progenitor']],w[['Mature/activated']],paired=TRUE,exact=FALSE)
bt<-replicate(2000,median(sample(delta,length(delta),TRUE)))
ll<-data.table(cohort='GSE216009',definition='granulopoiesis5',n_pairs=nrow(w),
               median_difference=median(delta),CI_low=ci(bt)[1],CI_high=ci(bt)[2],P=test$p.value)

fwrite(rr,file.path(OUT,'granulopoiesis_cohort_associations.csv'))
fwrite(rbindlist(cover),file.path(OUT,'granulopoiesis_marker_coverage.csv'))
fwrite(rbindlist(subsummary),file.path(OUT,'granulopoiesis_subtype_summary.csv'))
fwrite(sg,file.path(OUT,'granulopoiesis_subtype_overall.csv'))
fwrite(sp,file.path(OUT,'granulopoiesis_subtype_pairwise.csv'))
fwrite(ll,file.path(OUT,'granulopoiesis_localization.csv'))

# ---- numeric comparison against the frozen tables --------------------------
# The frozen tables carry the manuscript-reported targets. The point estimates and
# non-interval quantities enumerated below are re-derived here and asserted against them.
#
# Bootstrap percentile intervals are recomputed and written to the stage outputs
# (granulopoiesis_cohort_associations.csv, granulopoiesis_subtype_pairwise.csv,
# granulopoiesis_localization.csv), but their endpoints are NOT asserted: they have not been
# shown to reproduce the manuscript values exactly in this environment. The interval endpoints
# in checks/ remain the manuscript-reported targets and are deliberately excluded from the
# assertion loops below.
e<-as.data.table(expected('module_assoc'))
stopifnot(nrow(e)==7L,nrow(rr)==7L)
setkey(e,cohort,definition);r<-copy(rr);setkey(r,cohort,definition)
stopifnot(identical(e[,.(cohort,definition)],r[,.(cohort,definition)]))
for(nm in c('n','events','rho','P','universe','marker_n','FDR_BH7')){
 for(i in seq_len(nrow(r))){
  check('C_GRANULOPOIESIS_ASSOC',r$cohort[i],nm,
        r[[nm]][i],e[[nm]][i],cfg$tolerance,
        inputs=c(paste0('expression_',r$cohort[i]),r$cohort[i]),
        n=r$n[i],deaths=r$events[i],
        estimand='within-sample rank score vs z(Score16); Spearman; 2000 patient bootstraps')
 }
}

e2<-as.data.table(expected('module_localization'))
stopifnot(nrow(e2)==1L,nrow(ll)==1L)
for(nm in c('n_pairs','median_difference','P')){
 check('C_MODULE_LOCALIZATION',ll$definition,nm,ll[[nm]][1],e2[[nm]][1],cfg$tolerance,
       inputs='gse216009_pseudobulk_expression',n=ll$n_pairs,
       estimand='paired acute-sepsis donors; median of within-donor differences; 2000 bootstraps')
}

e3<-as.data.table(expected('subtype_summary'))
ss<-rbindlist(subsummary)
stopifnot(nrow(e3)==nrow(ss))
setkey(e3,cohort,label,subtype,definition);s3<-copy(ss);setkey(s3,cohort,label,subtype,definition)
stopifnot(identical(e3[,.(cohort,label,subtype,definition)],s3[,.(cohort,label,subtype,definition)]))
for(nm in c('n','events','median','q25','q75')){
 for(i in seq_len(nrow(s3))){
  check('C_SUBTYPE_MODULE_SUMMARY',paste(s3$subtype[i],s3$definition[i],sep='|'),nm,
        s3[[nm]][i],e3[[nm]][i],cfg$tolerance,
        inputs=c(s3$cohort[i],'mars_joined'),
        n=s3$n[i],estimand='module score by published subtype group')
 }
}

e4<-as.data.table(expected('subtype_overall'))
setkey(e4,cohort,definition);s4<-copy(sg);setkey(s4,cohort,definition)
stopifnot(nrow(e4)==nrow(s4),identical(e4[,.(cohort,definition)],s4[,.(cohort,definition)]))
for(nm in c('H','P')){
 for(i in seq_len(nrow(s4))){
  check('C_SUBTYPE_MODULE_OVERALL',paste(s4$cohort[i],s4$definition[i],sep='|'),nm,
        s4[[nm]][i],e4[[nm]][i],cfg$tolerance,
        inputs=c(s4$cohort[i],'mars_joined'),
        estimand='Kruskal-Wallis across published subtype groups')
 }
}

e5<-as.data.table(expected('subtype_pairwise'))
setkey(e5,cohort,definition,g1,g2);s5<-copy(sp);setkey(s5,cohort,definition,g1,g2)
stopifnot(nrow(e5)==nrow(s5),identical(e5[,.(cohort,definition,g1,g2)],s5[,.(cohort,definition,g1,g2)]))
for(nm in c('n1','n2','Cliff_delta','P','Holm_within_cohort_definition')){
 for(i in seq_len(nrow(s5))){
  check('C_SUBTYPE_MODULE_PAIRWISE',paste(s5$g1[i],s5$g2[i],sep='|'),nm,
        s5[[nm]][i],e5[[nm]][i],cfg$tolerance,
        inputs=c(s5$cohort[i],'mars_joined'),
        n=s5$n1[i]+s5$n2[i],
        estimand='pairwise Cliff delta; 1000 group bootstraps; Holm within cohort x definition')
 }
}

cat('GRANULOPOIESIS',nrow(rr),'cohort associations,',nrow(ss),'subtype summaries,',
    nrow(s5),'pairwise,',nrow(ll),'localization\n');flush.console()
