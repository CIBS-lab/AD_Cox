#load library
library(survminer)
library(survival)
library(MatchIt)
library(tibble)
library(sjPlot)
library(naniar) #setting NA values
library(finalfit)
library(broom) #for tidy()
library(coxme) #frailty model
library(stargazer)
library(survcomp) #for log partial likelihood of an extendedcox model
library(arsenal) #for table by
library(dplyr) # for filter()
library(ggplot2)
library(scater) #for multiplot()
library(ggstatsplot) #for bar plot with statistical test


#dataset pre matching (with comorbidities info)
dataset <- load("yourData.RData")



#-----------------matching

#match 3:1 with ALL the comorbidities confounders = age at recruitment, education, comorbidieties, CENTER
m.out1 <- matchit(ad ~ age + education + 
                    MI + CHF + PVD + Stroke + Pulmonary +
                    Rheumatic + PUD + LiverMild + Paralysis +
                    Renal + LiverSevere + centre, data = dataset, method = "nearest", distance = "glm", ratio = 3)
dataset_matched <- match.data(m.out1)
save(dataset_matched, file = 'dataset_matched.RData')


#statistics pre and post matching 
variables <- c('age_disc','sex', 'hyperlipidaemia', 'apoe_discrete', 'ad', 'ad_age', 'hypertention','diabetes', 'obesity', 'depression', 'flag_healthy', 'time')
pathway <- 'your_pathway/statistics_prematching.csv'
table_one <- tableby(ad~ ., data=dataset[ , (names(dataset) %in% variables)]) 
summary(table_one, text=TRUE,pfootnote=TRUE)
write.csv(as.data.frame(summary(table_one)),pathway)

pathway <- 'your_pathway/statistics_postmatching.csv'
table_one <- tableby(ad~ ., data=dataset_matched[ , (names(dataset_matched) %in% variables)]) 
summary(table_one, text=TRUE,pfootnote=TRUE)
write.csv(as.data.frame(summary(table_one)),pathway)



#-------------------------------------------------------------------------------all risk factors


#create time-depending variables
df <- subset(dataset_matched, select = c('eid', 'sex', 'apoe_discrete', 'ad',  'ad_age','hyperlipidaemia_age', 'status', 'hyperlipidaemia', 'diabetes', 'diabetes_age',
                                       'hypertention', 'hypertention_age','obesity', 'obesity_age','depression', 'depression_age', 'sex_apoe', 'time', 'centre', 'subclass', 'flag_healthy'))
newcgd <- tmerge(df, df, id=eid, tstop=time)
df_final <- tmerge(newcgd, df, id=eid, exposure_hlp = event(hyperlipidaemia_age) )
df_final <- tmerge(df_final, df, id=eid, exposure_t2d = event(diabetes_age))
df_final <- tmerge(df_final, df, id=eid, exposure_htn = event(hypertention_age))
df_final <- tmerge(df_final, df, id=eid, exposure_obesity = event(obesity_age))
df_final <- tmerge(df_final, df, id=eid, exposure_depr= event(depression_age))
df_final <-tmerge(df_final, df_final, eid, enum=cumtdc(tstart))

levels(df_final$sex) <- c('Female', 'Male')
df_final$exposure_hlp <- as.factor(df_final$exposure_hlp)
df_final$exposure_htn <- as.factor(df_final$exposure_htn)
df_final$exposure_t2d <- as.factor(df_final$exposure_t2d)
df_final$exposure_obesity<- as.factor(df_final$exposure_obesity)
df_final$exposure_depr <- as.factor(df_final$exposure_depr)
save(df_final, file = '/df_cphm_.RData')

#CPHM stratified with counting process
model_multi<-  coxph(Surv(tstart, tstop, status==1) ~ exposure_htn+exposure_t2d+exposure_hlp + exposure_obesity+ exposure_depr +sex + apoe_discrete + strata(subclass)  , data = df_final)
summary(model_multi)
exp(confint(model_multi))
save(model_multi, file = 'multi_model.RData')
tab_model(model_multi,file = "/result_multi_model.doc")
cph.hr <- broom::tidy(model_multi, conf.int = TRUE, exponetiante = TRUE)
write.csv(cph.hr, "/result_multi_model.csv")

#plotting HAZARD RATIO COEFFICIENTS
plot_models(model_multi,
            show.values = TRUE, show.p = TRUE, p.shape = TRUE, vline.color = "red",
axis.labels = c('Sex[Male]','APOE [e4]','Depression', 'Obesity', 'Hyperlipidemia', 'Diabetes', 'Hypertension'),
m.labels = c( '3:1 matched'),
colors = c( 'darkblue')
)+ ylim(.5,5) +
  labs(title = "Estimated Hazard Ratio ") +
  set_theme(
    geom.outline.size = 5, 
    geom.label.size = 6,
    geom.label.color = "black",
    title.color = "black", 
    title.size = 2, 
    axis.textcolor = "black", 
    base = theme_bw()
  )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file = 'HRs_multi_model.svg',
       width =8, height = 7, dpi = 300)



#PUT all the survival curvesl in a graph
t2d <- survfit( Surv(tstart, tstop, status==1) ~ exposure_t2d, data = filter(df_final,exposure_t2d ==1))
hlp <- survfit( Surv(tstart, tstop, status==1) ~ exposure_hlp, data = filter(df_final,exposure_hlp == 1))
htn <- survfit( Surv(tstart, tstop, status==1) ~ exposure_htn , data = filter(df_final,exposure_htn == 1))
obesity <- survfit( Surv(tstart, tstop, status==1) ~ exposure_obesity , data = filter(df_final,exposure_obesity == 1))
depr <- survfit( Surv(tstart, tstop, status==1) ~ exposure_depr , data = filter(df_final,exposure_depr == 1))
healthy <- survfit( Surv(tstart, tstop,status==1)~ flag_healthy, data = filter(df_final,flag_healthy == 1))
# Combine on the same plot
fit <- list(T2D_t = t2d, HLP_t = hlp, HTN_t = htn, Obesity_t = obesity, DEPR_t= depr, control = healthy)
ggsurvplot_combine(fit, df_final,
                   xlim = c(45,88),
                   break.time.by = 2,
                   palette=c('brown','gold', 'coral', 'darkviolet','red', 'black'),
                   legend.labs=c( "T2D", "HLP", "HTN", 'Obesity', 'Depression', 'Control'),
                   risk.table=FALSE
)

ggsave(file = 'surv_plot.svg',
       width =6.5, height = 5, dpi = 300)




#BUT betas(t)?
cox.zph(model_multi)
zph<- cox.zph(model_multi)$table
write.csv(zph, "/zph.csv")
par(mfrow=c(4,2))
plot_residuals1 <- ggcoxzph(cox.zph(model_multi)[1],ylim = c(-1,2), res = FALSE)
plot_residuals2 <- ggcoxzph(cox.zph(model_multi)[2],ylim = c(-1,2), res = FALSE)
plot_residuals3 <- ggcoxzph(cox.zph(model_multi)[3],ylim = c(-1,2), res = FALSE)
plot_residuals4 <- ggcoxzph(cox.zph(model_multi)[4],ylim = c(-1,2), res = FALSE)
plot_residuals5 <- ggcoxzph(cox.zph(model_multi)[5],ylim = c(-1,2), res = FALSE)
plot_residuals6 <- ggcoxzph(cox.zph(model_multi)[6],ylim = c(-1,2), res = FALSE)
plot_residuals7 <- ggcoxzph(cox.zph(model_multi)[7],ylim = c(-1,2), res = FALSE)

ggarrange(plot_residuals1$`1`,plot_residuals2$`1`,plot_residuals3$`1`,plot_residuals4$`1`,
          plot_residuals5$`1`, plot_residuals6$`1`, plot_residuals7$`1`,
          labels = c("A", "B", "C", "D", "E", "F", 'G'),
          ncol = 2, nrow = 4)
ggsave(file = '/residuals.pdf',
       width =8, height = 12, dpi = 300)


#if betas(t)--> age discretization
#BUt first convert all the exposures from factors to numerical
df_final$status <- as.factor(df_final$status)
df_final$exposure_htn <- as.numeric(df_final$exposure_htn) - 1
df_final$exposure_t2d <- as.numeric(df_final$exposure_t2d) - 1
df_final$exposure_obesity <- as.numeric(df_final$exposure_obesity) - 1
df_final$apoe_discrete <- as.numeric(df_final$apoe_discrete) - 1
df_final$exposure_hlp <- as.numeric(df_final$exposure_hlp) - 1
df_final$exposure_depr <- as.numeric(df_final$exposure_depr) - 1


#then split the dataset
data_splitted <- survSplit(Surv(tstart, tstop, status) ~ exposure_htn+exposure_t2d+exposure_hlp + exposure_obesity + exposure_depr
                           + sex + apoe_discrete , data= df_final, cut=c(62,72),episode= "tgroup", id="eid")
model_splitted <- coxph(Surv(tstart, tstop, status==1) ~ exposure_hlp:strata(tgroup)  + sex  + apoe_discrete:strata(tgroup) +    exposure_htn:strata(tgroup)+exposure_t2d:strata(tgroup)+ exposure_obesity:strata(tgroup) + exposure_depr:strata(tgroup), data=data_splitted)

summary(model_splitted)
cox.zph(model_splitted)
ggcoxzph(cox.zph(model_splitted) , ylim = c(-1,2))
tab_model(model_splitted,file = "/result_multi_model_stratified.doc")
cph.hr <- broom::tidy(model_splitted, conf.int = TRUE, exponetiante = TRUE)
write.csv(cph.hr, "/result_multi_model_stratified.csv")
save(model_splitted_m, file = 'multi_model_splitted.RData')

#check if betas(t) after age stratification?
cox.zph(model_splitted)
zph<- cox.zph(model_splitted)$table
write.csv(zph, "/zph_agestratified.csv")
plot_residuals1 <- ggcoxzph(cox.zph(model_splitted)[1],ylim = c(-1,2), res = FALSE)
plot_residuals2 <- ggcoxzph(cox.zph(model_splitted)[2],ylim = c(-1,2), res = FALSE)
plot_residuals3 <- ggcoxzph(cox.zph(model_splitted)[3],ylim = c(-1,2), res = FALSE)
plot_residuals4 <- ggcoxzph(cox.zph(model_splitted)[4],ylim = c(-1,2), res = FALSE)
plot_residuals5 <- ggcoxzph(cox.zph(model_splitted)[5],ylim = c(-1,2), res = FALSE)
plot_residuals6 <- ggcoxzph(cox.zph(model_splitted)[6],ylim = c(-1,2), res = FALSE)
plot_residuals7 <- ggcoxzph(cox.zph(model_splitted)[7],ylim = c(-1,2), res = FALSE)

ggarrange(plot_residuals1$`1`,plot_residuals2$`1`,plot_residuals3$`1`,plot_residuals4$`1`,
          plot_residuals5$`1`, plot_residuals6$`1`, plot_residuals7$`1`,
          labels = c("A", "B", "C", "D", "E", "F", 'G'),
          ncol = 2, nrow = 4)
ggsave(file = '/residuals_agestratified.pdf',
       width =8, height = 12, dpi = 300)


#comparing all the age discretized-model
plot_models(model_splitted_m1,
            show.values = TRUE, show.p = TRUE, p.shape = TRUE, vline.color = "red" ,
            m.labels = c( '3:1 matched'),
            colors = c('darkgreen')) + ylim(0.05,4) +
  labs(title = "Age -discretized Hazard Ratio ") +
  set_theme(
    geom.outline.size = 5, 
    geom.label.size = 6,
    geom.label.color = "black",
    title.color = "black", 
    title.size = 2.5, 
    axis.textcolor = "black", 
    base = theme_bw()
  )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file = 'HR_stratification_sex_apoe.svg',
       width = 10, height = 3, dpi = 300)


#plot survivals grouped by sex and APOE

t2d <- survfit( Surv(tstart, tstop, status==1) ~ exposure_t2d, data = filter(df_final4,exposure_t2d ==1 & sex_apoe == "female apoe4"))
hlp <- survfit( Surv(tstart, tstop, status==1) ~ exposure_hlp, data = filter(df_final4,exposure_hlp == 1& sex_apoe == "female apoe4"))
htn <- survfit( Surv(tstart, tstop, status==1) ~ exposure_htn , data = filter(df_final4,exposure_htn == 1& sex_apoe == "female apoe4"))
obesity <- survfit( Surv(tstart, tstop, status==1) ~ exposure_obesity , data = filter(df_final4,exposure_obesity == 1& sex_apoe == "female apoe4"))
depr <- survfit( Surv(tstart, tstop, status==1) ~ exposure_depr , data = filter(df_final4,exposure_depr == 1& sex_apoe == "female apoe4"))
healthy <- survfit( Surv(tstart, tstop,status==1)~ flag_healthy, data = filter(df_final4,flag_healthy == 1& sex_apoe == "female apoe4"))
# Combine on the same plot
fit <- list(T2D_t = t2d, HLP_t = hlp, HTN_t = htn, Obesity_t = obesity, Depr_t =depr, control = healthy)
plot_female_e4 <-ggsurvplot_combine(fit, df_final4,#fun = 'event',
                                    xlim = c(55,85),
                                    break.time.by = 5,
                                    palette=c('#1a9e76', '#f2cb27', '#dd9fc1', "#e7288a", '#a6761c',"#666666"),
                                    legend.labs=c( "T2D", "HLP", "HTN", 'BMI','DEP', 'Control Female e4'),
                                    surv.scale = 'percent',
                                    risk.table=FALSE,
                                    censor = FALSE,
                                    linewidth = 10
)

t2d <- survfit( Surv(tstart, tstop, status==1) ~ exposure_t2d, data = filter(df_final4,exposure_t2d ==1 & sex_apoe == "male apoe4"))
hlp <- survfit( Surv(tstart, tstop, status==1) ~ exposure_hlp, data = filter(df_final4,exposure_hlp == 1& sex_apoe == "male apoe4"))
htn <- survfit( Surv(tstart, tstop, status==1) ~ exposure_htn , data = filter(df_final4,exposure_htn == 1& sex_apoe == "male apoe4"))
obesity <- survfit( Surv(tstart, tstop, status==1) ~ exposure_obesity , data = filter(df_final4,exposure_obesity == 1& sex_apoe == "male apoe4"))
depr <- survfit( Surv(tstart, tstop, status==1) ~ exposure_depr , data = filter(df_final4,exposure_depr == 1& sex_apoe == "male apoe4"))
healthy <- survfit( Surv(tstart, tstop,status==1)~ flag_healthy, data = filter(df_final4,flag_healthy == 1& sex_apoe == "male apoe4"))
fit <- list(T2D_t = t2d, HLP_t = hlp, HTN_t = htn, Obesity_t = obesity, Depr_t =depr, control = healthy)

plot_male_e4 <- ggsurvplot_combine(fit1, df_final4,#fun = 'event',
                                   xlim = c(55,85),
                                   break.time.by = 5,
                                   palette=c('#1a9e76', '#f2cb27', '#dd9fc1', "#e7288a", '#a6761c',"#666666"),
                                   legend.labs=c( "T2D", "HLP", "HTN", 'Obesity','Depression', 'Control Male e4'),
                                   surv.scale = 'percent',
                                   risk.table=FALSE,
                                   censor = FALSE
)

t2d <- survfit( Surv(tstart, tstop, status==1) ~ exposure_t2d, data = filter(df_final4,exposure_t2d ==1 & sex_apoe == "female no apoe4"))
hlp <- survfit( Surv(tstart, tstop, status==1) ~ exposure_hlp, data = filter(df_final4,exposure_hlp == 1& sex_apoe == "female no apoe4"))
htn <- survfit( Surv(tstart, tstop, status==1) ~ exposure_htn , data = filter(df_final4,exposure_htn == 1& sex_apoe == "female no apoe4"))
obesity <- survfit( Surv(tstart, tstop, status==1) ~ exposure_obesity , data = filter(df_final4,exposure_obesity == 1& sex_apoe == "female no apoe4"))
depr <- survfit( Surv(tstart, tstop, status==1) ~ exposure_depr , data = filter(df_final4,exposure_depr == 1& sex_apoe == "female no apoe4"))
healthy <- survfit( Surv(tstart, tstop,status==1)~ flag_healthy, data = filter(df_final4,flag_healthy == 1& sex_apoe == "female no apoe4"))
fit <- list(T2D_t = t2d, HLP_t = hlp, HTN_t = htn, Obesity_t = obesity, Depr_t =depr, control = healthy)

plot_female_noe4 <-ggsurvplot_combine(fit2, df_final4,#fun = 'event',
                                      xlim = c(55,85),
                                      break.time.by = 5,
                                      palette=c('#1a9e76', '#f2cb27', '#dd9fc1', "#e7288a", '#a6761c',"#666666"),
                                      legend.labs=c( "T2D", "HLP", "HTN", 'Obesity','Depression', 'Control female noe4'),
                                      surv.scale = 'percent',
                                      risk.table=FALSE,
                                      censor = FALSE
)
t2d <- survfit( Surv(tstart, tstop, status==1) ~ exposure_t2d, data = filter(df_final4,exposure_t2d ==1 & sex_apoe == "male no apoe4"))
hlp <- survfit( Surv(tstart, tstop, status==1) ~ exposure_hlp, data = filter(df_final4,exposure_hlp == 1& sex_apoe == "male no apoe4"))
htn <- survfit( Surv(tstart, tstop, status==1) ~ exposure_htn , data = filter(df_final4,exposure_htn == 1& sex_apoe == "male no apoe4"))
obesity <- survfit( Surv(tstart, tstop, status==1) ~ exposure_obesity , data = filter(df_final4,exposure_obesity == 1& sex_apoe == "male no apoe4"))
depr <- survfit( Surv(tstart, tstop, status==1) ~ exposure_depr , data = filter(df_final4,exposure_depr == 1& sex_apoe == "male no apoe4"))
healthy <- survfit( Surv(tstart, tstop,status==1)~ flag_healthy, data = filter(df_final4,flag_healthy == 1& sex_apoe == "male no apoe4"))
fit <- list(T2D_t = t2d, HLP_t = hlp, HTN_t = htn, Obesity_t = obesity, Depr_t =depr, control = healthy)

plot_male_noe4 <-ggsurvplot_combine(fit3, df_final4,#fun = 'event',
                                    xlim = c(55,85),
                                    break.time.by = 5,
                                    palette=c('#1a9e76', '#f2cb27', '#dd9fc1', "#e7288a", '#a6761c',"#666666"),
                                    legend.labs=c( "T2D", "HLP", "HTN", 'Obesity', 'Depression', 'Control Male noe4'),
                                    surv.scale = 'percent',
                                    risk.table=FALSE,
                                    censor = FALSE
)
km.plot.list <- list(plot_female_e4,  plot_male_e4,plot_female_noe4, plot_male_noe4)

km.plot.list[[1]]$plot <- km.plot.list[[1]]$plot + labs(tag='A')
km.plot.list[[2]]$plot <- km.plot.list[[2]]$plot + labs(tag='C')
km.plot.list[[3]]$plot <- km.plot.list[[3]]$plot + labs(tag='B')
km.plot.list[[4]]$plot <- km.plot.list[[4]]$plot + labs(tag='D')
plot <- arrange_ggsurvplots(km.plot.list, print=TRUE, ncol=2, nrow=2)
ggsave("/survival_curve.svg", plot, dpi = 300, height = 8, width = 8.5)


#-------------------------------------------- survival by sex and APOE for each risk factor

#diabetes
graph <- survfit( Surv(tstart, tstop, status==1) ~ exposure_t2d + sex + apoe_discrete, data = filter(df_final4,exposure_t2d ==1))
p <-ggsurvplot(graph,
               linetype = c('sex'),
               palette = c(rep(c('#0066CC','#FF9933'),2),rep(c('blue4','blueviolet'),2)),
               #fun = 'event',
               xlim = c(50,85),
               break.time.by = 5,
               surv.scale = 'percent',
               censor =FALSE)
p_t2d <- p$plot + scale_linetype_discrete(name = "Sex", labels = c("Female", "Male")) 
ggsave(file = '/diabetes_curve.svg',
       width = 7, height =5, dpi = 300)

#hyertension
graph <- survfit( Surv(tstart, tstop, status==1) ~ exposure_htn + sex + apoe_discrete, data = filter(df_final4,exposure_htn ==1))
p <-ggsurvplot(graph,
               linetype = c('sex'),
               palette = c(rep(c('#0066CC','#FF9933'),2),rep(c('blue4','blueviolet'),2)),
               #fun = 'event',
               xlim = c(50,85),
               break.time.by = 5,
               surv.scale = 'percent',
               censor =FALSE)
p_htn <- p$plot + scale_linetype_discrete(name = "Sex", labels = c("Female", "Male")) 
ggsave(file = '/hypertension_curve.svg',
       width = 7, height =5, dpi = 300)

#hyperlipidemia
graph <- survfit( Surv(tstart, tstop, status==1) ~ exposure_hlp + sex + apoe_discrete, data = filter(df_final4,exposure_hlp ==1))
p <-ggsurvplot(graph,
               linetype = c('sex'),
               palette = c(rep(c('#0066CC','#FF9933'),2),rep(c('blue4','blueviolet'),2)),
               #fun = 'event',
               xlim = c(50,85),
               break.time.by = 5,
               surv.scale = 'percent',
               censor =FALSE)
p_hlp <- p$plot + scale_linetype_discrete(name = "Sex", labels = c("Female", "Male")) 
ggsave(file = '/hyperlipidemiq_curve.svg',
       width = 7, height =5, dpi = 300)

#depression
graph <- survfit( Surv(tstart, tstop, status==1) ~ exposure_depr + sex + apoe_discrete, data = filter(df_final4,exposure_depr ==1))
p <-ggsurvplot(graph,
               linetype = c('sex'),
               palette = c(rep(c('#0066CC','#FF9933'),2),rep(c('blue4','blueviolet'),2)),
               #fun = 'event',
               xlim = c(50,85),
               break.time.by = 5,
               surv.scale = 'percent',
               censor =FALSE)
p_depr <- p$plot + scale_linetype_discrete(name = "Sex", labels = c("Female", "Male")) 
ggsave(file = '/depression_curve.svg',
       width = 7, height =5, dpi = 300)

#obesity
graph <- survfit( Surv(tstart, tstop, status==1) ~ exposure_obesity + sex + apoe_discrete, data = filter(df_final4,exposure_obesity ==1))
p <-ggsurvplot(graph,
               linetype = c('sex'),
               palette = c(rep(c('#0066CC','#FF9933'),2),rep(c('blue4','blueviolet'),2)),
               #fun = 'event',
               xlim = c(50,85),
               break.time.by = 5,
               surv.scale = 'percent',
               censor =FALSE)
p_obesity <- p$plot + scale_linetype_discrete(name = "Sex", labels = c("Female", "Male")) 
ggsave(file = '/obesity_curve.svg',
       width = 7, height =5, dpi = 300)

plots <- multiplot(p_t2d, p_htn, p_hlp, p_depr, p_obesity, cols=2)
ggsave("/survivals_byDisease_sex_apoe.svg", plots, dpi = 300, height = 8, width = 7)



