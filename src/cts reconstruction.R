# Exact accepted author workflow; fixed reference-only RF fit, no target-outcome inputs.
# Core author commit 467fa897ba33b8ec355a8f33a1d3c4f4486846bd; original typo
# exp.c.core corrected to preceding exp.core.m, as in author's own printed call.
suppressPackageStartupMessages(library(randomForest))
exp.core.df<-read.csv(cfg$inputs$cts_core_expression,row.names=1,sep='',stringsAsFactors=TRUE)
exp.core.m<-data.matrix(exp.core.df)
core.cts<-read.csv(cfg$inputs$cts_core_labels,stringsAsFactors=TRUE);cts<-as.factor(core.cts$CTS)
exp.vanish.class<-read.csv(cfg$inputs$cts_target_expression,row.names=1,sep='',stringsAsFactors=TRUE)
stopifnot(nrow(exp.core.df)==18,ncol(exp.core.df)==460,nrow(core.cts)==460,identical(colnames(exp.core.df),make.names(as.character(core.cts$chipID))),setequal(levels(cts),c('1','2','3')),nrow(exp.vanish.class)==18,ncol(exp.vanish.class)==176,identical(rownames(exp.vanish.class),rownames(exp.core.df)),!anyNA(exp.core.m),!anyNA(exp.vanish.class))
set.seed(20260927)
rf<-randomForest::randomForest(cts~.,data=t(exp.core.m),proximity=TRUE)
pred<-predict(rf,t(exp.vanish.class),type='response')
stopifnot(length(pred)==176,identical(names(pred),paste0('VN',seq_len(176))),rf$ntree==500,rf$mtry==4)
fwrite(data.table(sample_id=names(pred),CTS_author_workflow=as.character(pred)),file.path(OUT,'cts_reconstructed176.csv'))
# Validation labels are loaded only after reference fitting and predictions are written.
# This minimized source has sample IDs and published CTS labels; no mortality/clinical fields.
truth<-fread(cfg$inputs$cts_public117_labels,select=c('sample_id','CTS_reference'))
stopifnot(nrow(truth)==117,!anyDuplicated(truth$sample_id),all(truth$sample_id%in%names(pred)))
matched<-as.character(pred[truth$sample_id]);correct<-sum(matched==as.character(truth$CTS_reference))
check('CTS_RECONSTRUCTION','published117','n',nrow(truth),117,0,inputs=c('cts_core_expression','cts_core_labels','cts_target_expression','cts_public117_labels'),estimand='Published CTS validation after reference-only randomForest fit; no target outcomes')
check('CTS_RECONSTRUCTION','published117','correct',correct,116,0,inputs='cts_public117_labels')
cm<-table(factor(truth$CTS_reference,levels=1:3),factor(matched,levels=1:3));e<-expected('cts_confusion')
for(i in 1:3)for(j in 1:3)check('CTS_CONFUSION',paste(i,j),paste0('pred_',j),as.integer(cm[i,j]),e[[paste0('pred_',j)]][i],0,inputs='cts_public117_labels')
cts_counts<-table(factor(pred,levels=1:3));e<-expected('cts_counts')
for(i in 1:3)check('CTS_COUNTS',paste0('CTS',i),'n',as.integer(cts_counts[i]),e$N[i],0,inputs='cts_target_expression')
fwrite(data.table(CTS=1:3,n=as.integer(cts_counts)),file.path(OUT,'cts_counts176.csv'))
fwrite(data.table(reference_CTS=1:3,pred1=cm[,1],pred2=cm[,2],pred3=cm[,3]),file.path(OUT,'cts_confusion117.csv'))
fwrite(data.table(seed=20260927,ntree=rf$ntree,mtry=rf$mtry,n_core=460,n_target=176,n_published=117,correct=correct,accuracy=correct/117,target_outcomes_used=FALSE),file.path(OUT,'cts_model_qc.csv'))
