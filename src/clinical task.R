# Restored original standardized pathway models; no pathway selection or figures.
d185<-as.data.frame(readinput('clinical_GSE185263'));d154<-as.data.frame(readinput('clinical_GSE154918'))
hallmark_185<-grep('^HALLMARK_',names(d185),value=TRUE);hallmark_154<-grep('^HALLMARK_',names(d154),value=TRUE)
stopifnot(length(hallmark_185)==50,setequal(hallmark_185,hallmark_154),nrow(d185)==345,nrow(d154)==91,!anyDuplicated(d185$patient_id),!anyDuplicated(d154$patient_id))
d185$age<-as.numeric(d185$age);d185$sofa_24h<-as.numeric(d185$sofa_24h);d185$sex<-factor(d185$sex);d185$center<-factor(d185$collection_location)
stopifnot(!anyNA(d185[,c('age','sex','center','sofa_24h',hallmark_185)]))
extract_pathway_lm <- function(df, pathway, model_id) {
  dat <- data.frame(
    sofa_24h = df$sofa_24h,
    pathway_z = as.numeric(scale(df[[pathway]])),
    age = df$age,
    sex = df$sex,
    center = df$center
  )
  if (model_id == "unadjusted") {
    formula_obj <- sofa_24h ~ pathway_z
    formula_text <- "SOFA_24h ~ pathway_z"
  } else if (model_id == "age_sex_adjusted") {
    formula_obj <- sofa_24h ~ pathway_z + age + sex
    formula_text <- "SOFA_24h ~ pathway_z + age + sex"
  } else if (model_id == "full_adjusted") {
    formula_obj <- sofa_24h ~ pathway_z + age + sex + center
    formula_text <- "SOFA_24h ~ pathway_z + age + sex + center"
  } else stop("Unknown model_id")
  fit <- lm(formula_obj, data = dat, na.action = na.omit)
  co <- summary(fit)$coefficients["pathway_z", ]
  ci <- confint(fit, "pathway_z", level = 0.95)
  data.frame(
    pathway = pathway,
    model = model_id,
    model_formula = formula_text,
    effect_scale = "beta in SOFA points per 1 SD higher pathway score",
    beta_per_1SD = unname(co["Estimate"]),
    SE = unname(co["Std. Error"]),
    CI95_low = unname(ci[1]),
    CI95_high = unname(ci[2]),
    P_value = unname(co["Pr(>|t|)"]),
    N_complete = nobs(fit),
    residual_df = df.residual(fit),
    adjusted_R2 = summary(fit)$adj.r.squared,
    effect_direction = ifelse(co["Estimate"] > 0, "higher pathway score associated with higher SOFA", "higher pathway score associated with lower SOFA"),
    model_status = ifelse(is.finite(co["Estimate"]), "converged", "failed"),
    stringsAsFactors = FALSE
  )
}


sofa<-rbindlist(lapply(hallmark_185,function(p)as.data.table(extract_pathway_lm(d185,p,'full_adjusted'))))
sofa[,P_BH_FDR_Hallmark50:=p.adjust(P_value,'BH')];fwrite(sofa,file.path(OUT,'clinical_SOFA50.csv'))
eall<-expected('clinical_sofa')
for(p in hallmark_185){q<-sofa[pathway==p];e<-eall[pathway==p];stopifnot(nrow(e)==1)
 for(nm in c('beta_per_1SD','SE','CI95_low','CI95_high','P_value','N_complete','P_BH_FDR_Hallmark50'))check('CLINICAL_SOFA',p,nm,q[[nm]],e[[nm]],inputs='clinical_GSE185263',n=345,estimand='SOFA24h ~ standardized pathway + age + sex + collection_location factor; BH50')
}
severity_map<-c('Healthy'=0,'Uncomplicated bacterial infection'=1,'Sepsis'=2,'Septic shock'=3)
d154$severity_rank<-unname(severity_map[as.character(d154$clinical_group)]);stopifnot(!anyNA(d154$severity_rank))
stopifnot(identical(as.integer(table(factor(d154$severity_rank,levels=0:3))),c(40L,12L,20L,19L)))
ordered<-rbindlist(lapply(hallmark_154,function(p){test<-suppressWarnings(cor.test(d154[[p]],d154$severity_rank,method='spearman',exact=FALSE));data.table(pathway=p,rho=unname(test$estimate),P_value=test$p.value,N=sum(complete.cases(d154[[p]],d154$severity_rank)))}))
ordered[,P_BH_FDR_Hallmark50:=p.adjust(P_value,'BH')];fwrite(ordered,file.path(OUT,'clinical_ordered50.csv'))
eall<-expected('clinical_ordered')
for(p in hallmark_154){q<-ordered[pathway==p];e<-eall[pathway==p];stopifnot(nrow(e)==1);for(nm in c('rho','P_value','N','P_BH_FDR_Hallmark50'))check('CLINICAL_ORDERED',p,nm,q[[nm]],e[[nm]],inputs='clinical_GSE154918',n=91,estimand='Spearman; fixed cross-sectional ordered groups; BH50')}
