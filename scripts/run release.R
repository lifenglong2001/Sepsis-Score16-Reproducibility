# Pure analysis entry: explicit inputs, new output directory, no network or installs.
root<-Sys.getenv('SEPSIS_CODE_ROOT')
if(!nzchar(root))stop('SEPSIS_CODE_ROOT must name this code release directory')
source(file.path(root,'src/common.R'),encoding='UTF-8')
if(!stage%in%c('core','benchmarks','pbmc','pbmc_score','scoring_portable','cellular','cross_task','clinical_task','cts_reconstruction','pbmc_replication','neutrophil_paired','granulopoiesis'))stop('Unknown stage: ',stage)
source(file.path(root,'src',paste0(gsub('_',' ',stage,fixed=TRUE),'.R')),encoding='UTF-8')
finish()
