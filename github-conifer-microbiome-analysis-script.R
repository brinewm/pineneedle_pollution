#packages
library(tidyverse)
library(tidytacos)
library(plotly)
library(RColorBrewer)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(Hmisc) #neutral model
library(stats4) #neutral model
library(minpack.lm) #neutral model
library(reshape2)
library(vegan)
library(svglite)
library(car)
outdir="./figures/"
today=str_remove_all(today(),"-")


# load dataset using tidytacos
taco <- read_tidytacos('~/Wenke-VOC/')

#add in metadata
info <- read.csv('sampledata.csv')

taco <- add_metadata(taco, info, table_type = 'sample') %>% 
  filter_taxa(kingdom == "Bacteria") %>% 
  select_taxa(-taxon,-species) %>%
  add_total_count() %>% 
  add_prevalence() %>%
  add_mean_rel_abundance() %>%
  add_taxon_name()

rm(info)

#filter by blanks so we can visualize them
taco_blank <- filter_samples(taco, pot == 'blank')
tacoplot_stack(taco_blank) +
  theme(text = element_text(size = 20)) + 
  ylab("Relative Abundance") + xlab('Sample ID')


#looks good!

#next we want to filter out all of the blanks & add in some useful information

taco<-taco%>%
  filter_samples(pot != "blank")

ta <- filter_samples(taco, pot != "blank") %>% 
  add_total_count() %>% 
  add_prevalence() %>% 
  add_mean_rel_abundance() %>% 
  mutate_taxa(log_ra = log10(mean_rel_abundance)) %>% 
  add_alphas() %>% 
  add_taxon_name()

tacosum(ta) #this gives an overview of the number of samples, taxa & reads we have

#add in some labels to make the figure prettier 
label_treatment <- c("Ambient Air", "Diesel", "Ozone", "Diesel + Ozone")
names(label_treatment) <- c('A', 'D', 'O', 'OD')

stack <- tacoplot_stack(ta) + 
xlab('Samples') + ylab("Relative Abundance") + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.text = element_text(face='bold'))+
  facet_grid(~treatment_order, scales = "free_x", labeller = labeller(treatment_order= label_treatment)) +
theme(text = element_text(size = 26)) 

#fitting neutral model to our occupancy vs. abundance data ---------
#this is code modified from Shade & Stopnisek, 2019 
#code found at https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/Core_prioritizing_script.R
#requires tidyverse, reshape2, vegan
N <- 1500
nReads=1500
spp <- ta %>% 
  filter_samples(total_count > N) %>% 
  tidytacos::rarefy(N, replace = TRUE)

otu <- spp %>% 
  counts_matrix(sample_name = sample_id,
                taxon_name = taxon_id,) %>% 
  t() #my counts matrix is the wrong way for the next steps so I transpose 

spp.bi <- 1*(otu >0) 
otu_occ <- rowSums(spp.bi)/ncol(spp.bi) #occupancy calculation
otu_rel <- apply(decostand(otu, method='total', MARGIN=2),1,mean) #relative abundance & standardization
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%  #combining occupancy and abundance data frame
  rownames_to_column('taxon_id')

PresenceSum <- data.frame(taxon_id = as.factor(row.names(otu)), otu) %>% 
  gather(sample_id, abun, -taxon_id) %>%
  left_join(ta$samples, by = 'sample_id') %>%
  group_by(taxon_id, ring) %>% #ring chosen instead of site 
  summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
            coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
            detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
  group_by(taxon_id) %>%
  summarise(sumF=sum(plot_freq),
            sumG=sum(coreSite),
            nS=length(ring)*2, #ring chosen instead site in original code
            Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='taxon_id') %>%
  transmute(otu=taxon_id,
            rank=Index) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start=otu_ranked$otu[1]
start_matrix <- as.matrix(otu[otu_start,])
start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)

for(i in 2:150){
  otu_add=otu_ranked$otu[i]
  add_matrix <- as.matrix(otu[otu_add,])
  add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))  
}

x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
elbow <- which.max(BC_ranked$fo_diffs)

  lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])))
s <- ggplot(BC_ranked[1:200,], aes(x=factor(BC_ranked$rank[1:200], levels=BC_ranked$rank[1:200]))) +
  geom_point(aes(y=proportionBC), pch=21, fill='#B8B8B8', size=3, colour='#666666', stroke = 0.2) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col= '#DB3069', cex=0.7) +
  geom_vline(xintercept=lastCall, lty=3, col='#765CEB', cex=0.7) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="#DB3069")+    
  annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="#765CEB") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
s

occ_abun$core <- 'no'
occ_abun$core[occ_abun$taxon_id %in% otu_ranked$otu[1:lastCall]] <- 'core'
occ_abun$core[occ_abun$taxon_id %in% otu_ranked$otu[1:elbow]] <- 'elbow'

spp=t(otu)
taxon=as.vector(rownames(otu))

#from burns
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

occ_abun_2 <- left_join(occ_abun, ta$taxa, by = 'taxon_id')

p <- ggplot(occ_abun_2, aes(x=otu_rel, y=otu_occ)) +
  geom_point(aes(colour=core, shape=core), size = 4, alpha=0.75, stroke = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c("#765CEB", '#DB3069', '#8DCEC3')) +
  geom_label_repel(aes(label=taxon_name), data=subset(occ_abun_2, otu_occ > 0.8), size=5, colour='azure4', alpha=0.5, force=10,box.padding = 0.6)+
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=obs.np$p), alpha=.5) +
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=obs.np$p), alpha=.3)+
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=obs.np$p), alpha=.3)+
  theme_classic() + 
  theme(text=element_text(size =23)) +
  xlab('Mean Relative Abundance (log10)') + ylab('Prevalence') +
  theme(axis.text = element_text(face='bold'), legend.position = 'none')

p

#betadisp to test for significant dispersion --------

abundances_matrix <- ta %>% 
  purrr::modify_at("samples", drop_na, treatment) %>% 
  rel_abundance_matrix()
metadata <- tibble(sample_id = rownames(abundances_matrix)) %>% 
  left_join(ta$samples, by = "sample_id")
abundances_matrix[1:5,1:5]
dist<-vegan::vegdist(abundances_matrix)

mod<-vegan::betadisper(dist,metadata$treatment)
anova(mod) #betadisp not significant

perform_adonis(ta, c("treatment", "ring", "total_count", "needle_weight"), by = 'terms') 
#this is significant! 

#total_count is significant. what happens if we remove samples with low read counts (<1000)
ta_filter <- ta %>% 
  filter_samples(total_count >=1500)

perform_adonis(ta_filter, c("treatment", "ring", "total_count", "needle_weight"), by = 'terms') 
#total_count is no longer significant 

#let's test the variance of the alpha diversity 
leveneTest(invsimpson ~ treatment, data = ta$samples)
# levene test is not significant so the data variance is okay

#let's check the normality using Shapiro
shapiro.test(ta$samples$invsimpson) #this is not significant, therefore normaly distributed

#we can do an ANOVA to see if the inverse simpson values are significantly different
res <- aov(data=ta$samples, invsimpson ~ treatment)
anova(res) #not significant

#let's plot our inverse simpson values ----

alpha_plot <- ggplot(ta$samples, aes(x=treatment_order, y=invsimpson, colour = treatment_order)) +
  geom_boxplot(size = 1.5, outliers=FALSE) + 
  geom_jitter(size=6) +
  theme_classic() +
  scale_colour_manual(values=c("#d95f02", "#7570b3", "#e7298a", "#1b9e77"))+
  theme(text = element_text(size = 23)) +
  theme(axis.text = element_text(face='bold'))+
  theme(legend.position = 'none') +
  xlab('Treatment') + ylab('Inverse Simpson Index') + 
  scale_x_discrete(labels=c("Ambient Air", "Diesel", "Ozone", "Ozone+Diesel")) 
alpha_plot

#species richness 
leveneTest(total_count ~ treatment, data = ta$samples) #data has normal variance
shapiro.test(ta$samples$obs) #data is not normally distributed 

#lets compare the means using kruskal-wallis
compare_means(obs ~ treatment, data = ta$samples, method = "kruskal.test") #not significant

#species richness plot
ggplot(ta$samples, (aes(x=treatment_order, y=obs, colour=treatment_order)))+
  geom_boxplot(size = 1.5, outliers=FALSE) +
  geom_jitter(size = 6) +
  theme_classic() +
  theme(text = element_text(size=23)) + 
  scale_colour_manual(values=c("#d95f02", "#7570b3", "#e7298a", "#1b9e77"))+
  ylab("Species Richness") + xlab("Treatment") +
  scale_x_discrete(labels=c("Ambient Air", "Diesel", "Ozone", "Ozone+Diesel")) +
 scale_y_continuous(limits = c(0,75), n.breaks = 6) +
  theme(axis.text = element_text(face = 'bold')) +
  theme(legend.position = 'none')

##betadisp testing ----

abundances_matrix <- ta %>%
  purrr::modify_at("samples", drop_na, treatment) %>%
  rel_abundance_matrix()
metadata <- tibble(sample_id = rownames(abundances_matrix)) %>%
  left_join(ta$samples, by = "sample_id")
abundances_matrix[1:5,1:5]
dist<-vegan::vegdist(abundances_matrix)

mod<-vegan::betadisper(dist,metadata$treatment)
anova(mod) #betadisp is not significant

# PERMANOVAs including pairwise comparisons ----

perform_adonis(ta, c("treatment", "ring", "total_count", "needle_weight")) #significant
t1<-ta %>% 
  filter_samples(total_count >= 1500) # consider removing low-read samples because total_count doesn't affect alpha diversity
perform_adonis(t1, c("treatment", "ring", "total_count", "needle_weight"), by='terms')#,"Emitter_VOCs"))
#perform_adonis(t1,c("Emitter_VOCs","treatment", "ring"))
t2<-filter_samples(t1,treatment%in%c("O","C")) 
perform_adonis(t2,c("treatment", "ring")) # ozone different from control

t2<-filter_samples(t1,treatment%in%c("D","C"))
perform_adonis(t2,c("treatment", "ring")) # diesel not different from control

t2<-filter_samples(t1,treatment%in%c("B","C")) 
perform_adonis(t2,c("treatment", "ring")) # ozone+diesel are not different from control

t2<-filter_samples(t1,treatment%in%c("B","O")) 
perform_adonis(t2,c("treatment", "ring")) # ozone+diesel vs. ozone alone = not different

t2<-filter_samples(t1,treatment%in%c("B","D","C"))%>%
  mutate_samples(treatment=ifelse(treatment=="B","B","Other")) 
perform_adonis(t2,c("treatment", "ring")) # Diesel & Ozone very different from diesel and control!

rm(t1,t2)

### Which bacteria? Differential abundance analysis with radEmu -----

sort(ta$samples$total_count)
y1<-ta%>%
  aggregate_taxa(rank = "genus") %>% 
  add_taxon_name() %>%
  add_mean_rel_abundance()%>%
  mutate_samples(both=ifelse(treatment=="B","1","0")) %>% 
  mutate_samples(diesel=ifelse(treatment=="D","1","0")) %>% 
  mutate_samples(ozone=ifelse(treatment=="O","1","0")) %>% 
  add_prevalence()

##diesel 
samp<-dplyr::select(y1$samples,sample_id,treatment_order,ring,diesel,ozone,both)
rownames(samp)<-samp$sample_id
#y<-tidy_count_to_matrix(y2$counts)[samp$sample_id,]
y <- counts_matrix(y1)[samp$sample_id,]
radEmu::make_design_matrix(formula=~diesel, data=samp) #treatment D, k=2

fit2 <- radEmu::emuFit(formula = ~ diesel, data=samp, Y = y, 
                       cluster = samp$ring, 
                       test_kj = data.frame(k=2, j=c(1:96)))

#ozone
radEmu::make_design_matrix(formula=~ozone, data=samp) #treatment D, k=2
fit3 <- radEmu::emuFit(formula = ~ ozone, data=samp, Y = y, 
                       cluster = samp$ring, 
                       test_kj = data.frame(k=2, j=c(1:96)))

#both
radEmu::make_design_matrix(formula=~both, data=samp) #treatment D, k=2
fit4 <- radEmu::emuFit(formula = ~ both, data=samp, Y = y, 
                       cluster = samp$ring, 
                       test_kj = data.frame(k=2, j=c(1:96)))


diesel <- fit2$coef %>%
  filter(pval<0.1)%>%
  mutate(taxon_id=category)%>%
  left_join(y1$taxa)%>%
  mutate(name=taxon_name)

diesel$signif <- ifelse(diesel$pval < 0.05, "*", "") #add in some asterisks to indicate significance

plot_diesel <- diesel %>% 
  ggplot(aes(x=estimate, y=name))+ 
  geom_errorbar(aes(xmin = estimate-se, xmax = estimate+se), linewidth=1.2, width = 0.1, colour = "#ADADAD") +
  geom_point(colour = "#7570B3", aes(size=30))+
  theme_bw()+ 
  theme(text = element_text(size = 20)) +
  geom_text(
    data = subset(diesel),
    aes(label = signif),
    size = 7, 
    hjust=-3,vjust=0.7)+
  xlab('Estimate Statistic') + ylab('Taxon') +
  theme(legend.position = 'none') +
  geom_vline(xintercept = 0, linetype='solid', colour = "#E0E0E0", size=1) 

plot_diesel


ozone <- fit3$coef %>%
  filter(pval<0.1)%>%
  mutate(taxon_id=category)%>%
  left_join(y1$taxa)%>%
  mutate(name=taxon_name)

ozone$signif <- ifelse(ozone$pval < 0.05, "*", "")

plot_ozone <- ozone %>% 
  ggplot(aes(x=estimate, y=name))+ 
  geom_errorbar(aes(xmin = estimate-se, xmax = estimate+se), linewidth = 1.2, width = 0.1, colour = "#ADADAD") +
  geom_point(colour = "#E7298A", aes(size=30))+
  theme_bw()+ 
  theme(text = element_text(size = 10)) +
  geom_text(
    data = subset(ozone),
    aes(label = signif),
    size = 7, 
    hjust=-1,vjust=0.7)+
  xlab('Estimate Statistic') + ylab('Taxon') +
  theme(legend.position = 'none') +
  geom_vline(xintercept = 0, linetype='solid', colour = "#E0E0E0", size=1) 

plot_ozone


both <- fit4$coef %>%
  filter(pval<0.1)%>%
  mutate(taxon_id=category)%>%
  left_join(y1$taxa) %>%
  mutate(name=taxon_name)

both$signif <- ifelse(both$pval < 0.05, "*", "")

plot_both <- both %>% 
  ggplot(aes(x=estimate, y=name))+ 
  geom_errorbar(aes(xmin = estimate-se, xmax = estimate+se), size = 1.2, width = 0.1, colour = "#ADADAD") +
  geom_point(colour = "#1B9E77", aes(size=30))+
  theme_bw()+ 
  theme(text = element_text(size = 20)) +
  geom_text(
    data = subset(both),
    aes(label = signif),
    size = 7, 
    hjust=-3,vjust=0.7)+
  xlab('Estimate Statistic') + ylab('Taxon') +
  theme(legend.position = 'none') +
  geom_vline(xintercept = 0, linetype='solid', colour = "#E0E0E0", size=1) 

plot_both

#### it would be better if we had all of these side-by-side especially for a figure
# let's combine into a single object


all <- full_join(ozone, diesel)
all <- full_join(all, both) %>% 
  mutate(name=taxon_name)

#if we plot everything it's a bit overwhelming, let's try and filter the data a bit to just show a subset

sigtax <- all %>% 
  filter(pval<0.05)
gvl <- all %>%
  filter(category %in% sigtax$category) %>%
  mutate(taxon_id=category)%>%
  left_join(all)%>%
  filter(name!="other")


label_emu <- c("Ozone", "Diesel", "Ozone+Diesel")
names(label_emu) <- c('treatmentO', 'treatmentD', 'treatmentC')

plot <- gvl %>% 
  ggplot(aes(x=estimate, y=name, colour=covariate))+ 
  geom_errorbar(aes(xmin = estimate-se, xmax = estimate+se), size = 1.5, width = 0.2, colour = "#8F8F8F") +
  geom_point(aes(size=14))+
  theme_bw()+ 
  scale_colour_manual(values=c("#1B9E77","#7570B3","#E7298A")) +
  theme(text = element_text(size = 20),
        axis.text.y=element_text(face='italic')) +
  geom_text(
    data = subset(gvl),
    aes(label = signif, colour=covariate),
    size = 10, 
    hjust=-3,vjust=0.7, fontface="bold")+
  xlab('Estimate Statistic') + ylab('Taxon') +
 theme(legend.position = 'none') +
  geom_vline(xintercept = 0, linetype='solid', colour = "#E0E0E0", size=1) +
  facet_wrap(~covariate,labeller = labeller(covariate= label_emu)) 

plot 

#in illustrator I changed the order of the panels and removed the minor gridlines :) 

###### DAtest 
#inputs: count table (samples as columns and features as rows) - assumed to be raw & compositional. 
#Advised to trim abundant features from dataset. Second input is the predictor of interest (e.g. control/treatment). 
#Vector in same order as the columns in the count table 

devtools::install_github("Russel88/DAtest")
library(DAtest)
library(DESeq2)
library(limma)
library(edgeR)
library(metagenomeSeq)
library(baySeq)
library(ALDEx2)
library(impute)
library(ANCOMBC)
library(samr)
library(pscl)
library(statmod)
library(mvabund)
library(eulerr)
library(emmeans)
BiocManager::install(c("DESeq2","limma","edgeR","metagenomeSeq","baySeq","ALDEx2","impute","ANCOMBC"))

#prepare our predictor data for input into the formula
conds <- y1$samples$treatment_order
diesel <- y1$samples$diesel
ring <- y1$samples$ring
ozone <- y1$samples$ozone
both <- y1$samples$both

## note using the same input as from radEmu above 
# this is our tidytacos object with the taxa aggregating at the genus level 
counts <- counts_matrix(y1, taxon_name = taxon_name) %>% 
  t() %>% 
  as.data.frame(rownames=sample_id) %>% 
  relocate(s9, .before = s10) #relocate because otherwise s10 is the first sample

test_new <- testDA(counts, predictor = ozone, R = 75, k=c(1,2,3)) 
#swap out the predictor to include ozone, diesel & both
summary(test_new)

summary_plot <- plot(test_new, sort = 'Score') #we can also visualize the scores for the tests


testa_new <- allDA(counts, predictor = ozone, 
                   tests = c("qpo", "per", "ere")) #let's run the 3 tests with scores > 0

vennDA(testa_new, tests = c("qpo", "per", "ere")) #we can now visualize which tests have shared significant features
#only ere has significant features so let's just run this test and visualize it

ere <- DA.ere(counts, predictor=ozone)

ere$signif <- ifelse(ere$pval.adj < 0.05, "*", "")

plot <- ggplot(ere, aes(x=logFC, y=Feature)) +
  geom_point(aes(colour=signif)) +
  geom_text(
    data = subset(ere),
    aes(label = signif),
    size = 8, 
    hjust=-1,vjust=0.7, 
    colour="#8075FF") +
  theme_bw() +
  scale_colour_manual(values=c("#BAC9BA","#695CFF"))+
  theme(legend.position = "none")

plot

