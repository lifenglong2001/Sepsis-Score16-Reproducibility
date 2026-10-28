# Donor is the experimental unit. No cell-level pseudo-replication.
st<-readinput('pbmc_scores');stopifnot(!anyDuplicated(st[,.(donor_id,cell_type)]))
check('A_PBMC_DONORS','GSE279452','unique_donors',uniqueN(st$donor_id),108,0,inputs='pbmc_scores')
ct_order<-c('Monocyte','CD4 T','CD8 T','B','NK');stopifnot(setequal(unique(st$cell_type),ct_order))
eall<-expected('pbmc');metrics<-c(unique(eall$metric),'Score16');stopifnot(setequal(setdiff(metrics,'Score16'),c(pos,neg)))
st[,Score16:=raw_score16(st)]
cliff<-function(x,y)mean(sign(outer(x,y,'-')))
seed<-2026101701;B<-cfg$bootstrap_B;set.seed(seed);results<-list();i<-0L
for(ct in ct_order)for(control in c('Healthy control','ICU non-sepsis control'))for(m in metrics){
 d<-st[cell_type==ct];x<-d[disease_class=='Sepsis'][[m]];y<-d[disease_class==control][[m]]
 ok<-length(x)>=3&&length(y)>=3
 if(ok){wt<-suppressWarnings(wilcox.test(x,y,exact=FALSE));dd<-cliff(x,y);bt<-replicate(B,cliff(sample(x,length(x),TRUE),sample(y,length(y),TRUE)));ci<-quantile(bt,c(.025,.975),type=6,names=FALSE)}else{wt<-NULL;dd<-NA_real_;ci<-c(NA_real_,NA_real_)}
 i<-i+1L;results[[i]]<-data.table(cell_type=ct,contrast=if(control=='Healthy control')'Sepsis vs Healthy' else 'Sepsis vs ICU non-sepsis',metric=m,n_sepsis=length(x),n_control=length(y),cliff_delta=dd,CI_low=ci[1],CI_high=ci[2],P=if(ok)wt$p.value else NA_real_,median_difference=median(x)-median(y),B=B,seed=seed)
}
r<-rbindlist(results);r[,FDR_16:=NA_real_];r[metric!='Score16',FDR_16:=p.adjust(P,'BH'),by=contrast];r[,FDR_with_Score16:=p.adjust(P,'BH'),by=contrast];r[,clear_signal:=is.finite(FDR_16)&FDR_16<.05&abs(cliff_delta)>=.33]
fwrite(r,file.path(OUT,'pbmc_primary_refitted.csv'))
for(j in which(r$metric!='Score16')){
 q<-r[j];e<-eall[cell_type==q$cell_type&contrast==q$contrast&metric==q$metric]
 for(nm in c('n_sepsis','n_control','cliff_delta','CI_low','CI_high','P','median_difference','FDR_16','FDR_with_Score16'))check(paste(q$cell_type,q$metric,sep=':'),q$contrast,nm,q[[nm]],e[[nm]],inputs='pbmc_scores',n=q$n_sepsis+q$n_control,estimand='donor-level Wilcoxon; Cliff delta; 80-pathway FDR per contrast')
}
for(cc in unique(r$contrast)){
 q<-r[contrast==cc&metric!='Score16'];check('B_PBMC_FDR_FAMILY',cc,'tests',nrow(q),80,0,inputs='pbmc_scores')
 check('B_PBMC_FDR_COUNT',cc,'FDR_lt_005',sum(q$FDR_16<.05),if(cc=='Sepsis vs Healthy')48 else 0,0,inputs='pbmc_scores')
}
