################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(viridis)
library(cowplot)
library(wesanderson)
library(GGally)
library(MASS)
library(extrafont)
library(ggtree)
library(tidytree)
library(dplyr)
library(tidyr)
################################################# Read in the data and mash it together #################################################
setwd("/Users/Cliff/Documents/Kosuri Lab/R and python/ancestral tree/data_for_R")

numpars <- read.csv("200929_num_paralogs_till_spec.csv")
anctodiv <- read.csv("200722_anc_to_diverge_ints.csv")
allpars <- read.csv("200722_allpar_score_and_dist.csv")
homos <- read.csv('200930_homodimer_counting.csv')
disttointprof <- read.csv('200930_interactionprofiles_cors.csv')

# g.paralogsbeforespec <- ggplot(numpars, aes(X.numpars..)) + geom_histogram(bins = 4, fill = "#0091F4") + labs(x = "Number of paralogs before specificity is gained", y = "Count", title = "200722 Specificity is not gained immediately") +
#   theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18))
# g.paralogsbeforespec
# 
# ggsave(g.paralogsbeforespec, file = "200722_num_paralogs_till_specificity.pdf", scale = 0.6)
# 
# 

g.specgains <- ggplot(numpars, aes(X.distance.)) + geom_histogram(bins = 8, fill = "#0091F4") + labs(x = "Branch Length", y = "Count", title = "201001 EA66k distance needed to gain specificty") +
  scale_x_continuous(limits = c(0, 1.4))
g.specgains

ggsave(g.specgains, file = '../figures/201001 EA66k distance needed to gain specificty.pdf', scale = 0.6)


colnames(allpars) <- c("xpep", "ypep","intscore","distance")
colnames(homos) <- c("xpep", "ypep","intscore","extant")


#allpars$distance <- round(allpars$distance, digits = 1)
g.allparsinthist <- ggplot(allpars, aes(as.factor(intscore), fill = as.factor(intscore))) + geom_bar(position = "stack", show.legend = F) + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131")) +
  scale_x_discrete(labels = c("Not measured", "None", "Weak", "Strong")) + labs(x = "Interaction type", y = "Count", title = "200930 EA66k All paralogs by interaction type counts")
g.allparsinthist
ggsave(g.allparsinthist, file = "../figures/200930 Paralogs by interaction type.pdf", scale = 0.6)


g.homosinthist <- ggplot(homos, aes(as.factor(intscore), fill = as.factor(intscore))) + geom_bar(position = "stack", show.legend = F) + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#FE3131")) +
  scale_x_discrete(labels = c("Not measured", "None", "Strong")) + labs(x = "Interaction type", y = "Count", title = "200930 EA66k All Homodimers by interaction type counts")
g.homosinthist
ggsave(g.homosinthist, file = "../figures/200930 homodimers by interaction type.pdf", scale = 0.6)

allpars <- data.frame(allpars, "pairtype" = "para")
allpars <- select(allpars, -distance)
homos <- data.frame(homos, "pairtype" = "homo")
homos <- select(homos, -extant)
parshomo <- bind_rows(allpars, homos)

g.hpinthist <- ggplot(parshomo, aes(pairtype, fill = as.factor(intscore))) + geom_bar(position = "stack") + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("Not measured", "None", "Weak", "Strong")) +
  scale_x_discrete(labels = c("Homodimers","Heterodimers")) + labs(x = "Interaction type", y = "Count", title = "200930 EA66k Homodimers and paralogs by interaction type counts")
g.hpinthist
ggsave(g.hpinthist, file = "../figures/200930 homodimers and paralogs by interaction type.pdf", scale = 0.6)

g.hpinthistper <- ggplot(parshomo, aes(pairtype, fill = as.factor(intscore))) + geom_bar(position = "fill") + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("Not measured", "None", "Weak", "Strong")) +
  scale_x_discrete(labels = c("Homodimers","Heterodimers")) + labs(x = "Interaction type", y = "Percentage", title = "200930 EA66k Homodimers and paralogs by interaction type percentages")
g.hpinthistper

ggsave(g.hpinthistper, file = "../figures/200930 homodimers by interaction type percentage.pdf", scale = 0.6)

colnames(anctodiv) <- c('xpep', 'ypep','distance','anc_to_div')
anctodiv <- separate(anctodiv, anc_to_div, into = c("xpep_anc_int","ypep_anc_int"), sep = "\\+")
anctodiv <- data.frame(anctodiv, "ancints" = as.character(c("")))
anctodiv$ancints <- as.character(anctodiv$ancints)
anctodiv$xpep_anc_int <- as.numeric(anctodiv$xpep_anc_int)
anctodiv$ypep_anc_int <- as.numeric(anctodiv$ypep_anc_int)
anctodiv$ancints[anctodiv$xpep_anc_int == 1.0 & anctodiv$ypep_anc_int == 1.0] <- "Both"
anctodiv$ancints[anctodiv$xpep_anc_int == 1.0 & anctodiv$ypep_anc_int == 0.5] <- "Both"
anctodiv$ancints[anctodiv$xpep_anc_int == 1.0 & anctodiv$ypep_anc_int == 0.0] <- "one"
anctodiv$ancints[anctodiv$xpep_anc_int == 0.0 & anctodiv$ypep_anc_int == 1.0] <- "one"

g.changesforspec <- ggplot(anctodiv, aes(ancints, fill = ancints)) + geom_bar(position = "stack", show.legend = F) + scale_fill_manual(values = c("#FE3131", "#0090F3")) +
  scale_x_discrete(labels = c("Both", "Single")) + labs(x = "How many branches specificity changes occured on", y = "Count", title = "200930 EA66k Changes on which branches led to specificity")
g.changesforspec
ggsave(g.changesforspec, file = "../figures/200930 EA66k branch changes to spec.pdf", scale = 0.6)

g.changesspecdist <- ggplot(anctodiv, aes(ancints, distance, fill = ancints)) + geom_dotplot(binaxis = "y", show.legend = F)+ scale_fill_manual(values = c("#FE3131", "#0090F3")) +
  scale_x_discrete(labels = c("Both", "Single")) + labs(x = "How many branches specificity changes occured on", y = "Distance to gain specificity", title = "200930 EA66k Changes on which branches led to specificity")
g.changesspecdist
ggsave(g.changesspecdist, file = "../figures/200930 EA66k branch changes to spec by dist.pdf", scale = 0.6)


colnames(disttointprof) <- c('xpep', 'ypep', 'medRD', 'distance', 'corr','logtype','pearsonlog','spearman','intscore','parentscore','paraparents','ancestor','isparentchild')
disttointprof <- filter(disttointprof, !is.na(corr))


distcor <- data.frame(x = c(0.2,3), y = c(-.2, 0.6), logtype = c("ortholog","paralogs"), label = c(paste("italic(r)==",round(cor(filter(disttointprof, logtype == 'ortholog')$distance, filter(disttointprof, logtype == 'ortholog')$pearsonlog), digits = 2)), 
                                                               paste("italic(r)==",round(cor(filter(disttointprof, logtype == 'paralogs')$distance, filter(disttointprof, logtype == 'paralogs')$pearsonlog), digits = 2))))
g.disttointpro <- ggplot(disttointprof,  aes(distance, pearsonlog, color = logtype)) + geom_point(alpha = 0.1) + geom_smooth(method = 'lm', formula = y~x) + scale_color_manual(values = c("#FE3131", "#0090F3"), name = "Homolog type", labels = c("Ortholog","Paralog")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Correlation of interaction profiles vs evo distance by homolog type') + geom_text(data = distcor, aes(x=x, y=y, label=label), parse = T, color = 'Black')
g.disttointpro
ggsave(g.disttointpro, file = "../figures/200930 EA66k dist vs correlation of interact profile by evo relation.pdf", scale = 0.6)


g.disttointpro <- ggplot(disttointprof,  aes(distance, pearsonlog, color = logtype)) + geom_point(alpha = 0.1) + geom_smooth(method = 'lm', formula = y~x) + scale_color_manual(values = c("#FE3131", "#0090F3"), name = "Homolog type", labels = c("Ortholog","Paralog")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Correlation of interaction profiles vs evo distance by homolog type') + geom_text(data = distcor, aes(x=x, y=y, label=label), parse = T, color = 'Black') +
  facet_wrap(~logtype)
g.disttointpro 
ggsave(g.disttointpro, file = "../figures/200930 EA66k dist vs correlation of interact profile by evo relation wrap.pdf", scale = 0.6)


g.disttointprobystr <- ggplot(disttointprof,  aes(distance, pearsonlog, color = as.factor(intscore))) + geom_point(alpha = 0.1) + geom_smooth(method = 'lm', formula = y~x, se = F) + scale_color_manual(values = c("#A7A7A7","#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("Not measured", "None", "Weak","Strong")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Correlation of interaction profiles vs evo distance by interaction strength')
g.disttointprobystr
ggsave(g.disttointprobystr, file = "../figures/200930 EA66k dist vs correlation of interact profile by evo relation int strength.pdf", scale = 0.6)


g.disttointprobystrpara <- ggplot(filter(disttointprof, logtype == "paralogs"),  aes(distance, pearsonlog, color = as.factor(intscore))) + geom_point(alpha = 0.3) + geom_smooth(method = 'lm', formula = y~x, se = F) + scale_color_manual(values = c("#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("None", "Weak","Strong")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Paralogs correlation of interaction profiles vs evo distance by interaction strength')
g.disttointprobystrpara
ggsave(g.disttointprobystrpara, file = "../figures/200930 EA66k paras dist vs correlation of interact profile by evo relation int strength.pdf", scale = 0.6)


####spearmans
distcors <- data.frame(x = c(0.2,3), y = c(-.2, 0.6), logtype = c("ortholog","paralogs"), label = c(paste("italic(r)==",round(cor(filter(disttointprof, logtype == 'ortholog')$distance, filter(disttointprof, logtype == 'ortholog')$spearman), digits = 2)), 
                                                                                                   paste("italic(r)==",round(cor(filter(disttointprof, logtype == 'paralogs')$distance, filter(disttointprof, logtype == 'paralogs')$spearman), digits = 2))))
g.disttointpros <- ggplot(disttointprof,  aes(distance, spearman, color = logtype)) + geom_point(alpha = 0.1) + geom_smooth(method = 'lm', formula = y~x) + scale_color_manual(values = c("#FE3131", "#0090F3"), name = "Homolog type", labels = c("Ortholog","Paralog")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Spearmans of interaction profiles vs evo distance by homolog type') + geom_text(data = distcors, aes(x=x, y=y, label=label), parse = T, color = 'Black')
g.disttointpros
ggsave(g.disttointpros, file = "../figures/200930 EA66k dist vs spearman of interact profile by evo relation.pdf", scale = 0.6)


g.disttointpros <- ggplot(disttointprof,  aes(distance, spearman, color = logtype)) + geom_point(alpha = 0.1) + geom_smooth(method = 'lm', formula = y~x) + scale_color_manual(values = c("#FE3131", "#0090F3"), name = "Homolog type", labels = c("Ortholog","Paralog")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Spearmans of interaction profiles vs evo distance by homolog type') + geom_text(data = distcors, aes(x=x, y=y, label=label), parse = T, color = 'Black') +
  facet_wrap(~logtype)
g.disttointpros
ggsave(g.disttointpros, file = "../figures/200930 EA66k dist vs spearman of interact profile by evo relation wrap.pdf", scale = 0.6)


g.disttointprobystrs <- ggplot(disttointprof,  aes(distance, spearman, color = as.factor(intscore))) + geom_point(alpha = 0.1) + geom_smooth(method = 'lm', formula = y~x, se = F) + scale_color_manual(values = c("#A7A7A7","#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("Not measured", "None", "Weak","Strong")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Spearman of interaction profiles vs evo distance by interaction strength')
g.disttointprobystrs
ggsave(g.disttointprobystrs, file = "../figures/200930 EA66k dist vs spearman of interact profile by evo relation int strength.pdf", scale = 0.6)


g.disttointprobystrparas <- ggplot(filter(disttointprof, logtype == "paralogs"),  aes(distance, spearman, color = as.factor(intscore))) + geom_point(alpha = 0.3) + geom_smooth(method = 'lm', formula = y~x, se = F) + scale_color_manual(values = c("#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("None", "Weak","Strong")) +
  labs(x = 'Distance', y = 'Correlation coeffecient between interaction profiles', title = '200930 EA66k Paralogs spearman of interaction profiles vs evo distance by interaction strength')
g.disttointprobystrparas
ggsave(g.disttointprobystrparas, file = "../figures/200930 EA66k paras dist vs correlation of interact profile by evo relation int strength.pdf", scale = 0.6)

colnames(homos) <- c('xpep','ypep','homoint','homoextant')
homos$homoint[homos$homoint == -1] <- -0.25
intsbyhomos <- left_join(disttointprof, homos, by = 'xpep')
intsbyhomos <- left_join(intsbyhomos, homos, by = c('ypep.x' = 'ypep'))
intsbyhomos <- mutate(intsbyhomos, 'combohomo' = homoint.x + homoint.y)

g.intsbyhomos <- ggplot(intsbyhomos, aes(as.factor(combohomo), fill= as.factor(intscore))) + geom_bar(position = "stack") + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131"), labels = c('Not measured', 'None', 'Weak','Strong'), name = 'Interaction\nstrength') +
  labs(x = "Number of homodimers interacting", y = "Count", title = "201002 EA66k Interaction strength by number of interacting homodimers") + scale_x_discrete(labels = c('Neither\nmeasured','One unmeasured\none nonhomodimer','None','One unmeasured\none homodimer','One','Both'))
g.intsbyhomos
ggsave(g.intsbyhomos, file ='../figures/201002 Interaction strength by homodimers.pdf')

paraintsbyhomos <- filter(intsbyhomos, logtype == "paralogs") #, homoextant.x == "Yes")
paraintsbyhomos <- separate(paraintsbyhomos, paraparents, into = c('parentx','parenty'), sep = '\\+')
paraintsbyhomos <- left_join(paraintsbyhomos, homos, by = c('parentx' = 'xpep'))
paraintsbyhomos <- left_join(paraintsbyhomos, homos, by = c('parenty' = 'ypep'))
paraintsbyhomos <- select(paraintsbyhomos, -homoextant.x, -ypep.y, -xpep.y, -homoextant.y, -homoextant.x.x, -homoextant.y.y, -ypep, -xpep)
colnames(paraintsbyhomos) <- c("xpep",'ypep','medRD','distance','corr','logtype','pearsonlog','spearman','intscore','parentscore','parentx','parenty','ancestor','homointx','homointy','combohomo','homointparentx','homointparenty')
paraintsbyhomos <- mutate(paraintsbyhomos, 'parentcombohomo' = homointparentx +  homointparenty)
#paraintsbyhomos <- mutate(paraintsbyhomos,)

sumparahomos <- group_by(paraintsbyhomos, as.factor(combohomo)) %>%
  summarise('strongints' = sum(intscore == 1.0), 'weakints' = sum(intscore == 0.5), 'nonints' = sum(intscore == 0.0), 'notmeasured' = sum(is.na(intscore) | intscore == -1.0))

sumparahomos <- data.frame(sumparahomos, "homodimers" = c('NA and none', 'Both none', 'NA and strong','Strong and None', 'Both strong'))
sumparahomos <- select(sumparahomos, homodimers, strongints, weakints, nonints)

thinsumparahomos <- gather(sumparahomos, key = 'inttype', value = 'count', -homodimers)
thinsumparahomos$homodimers <- factor(thinsumparahomos$homodimers, levels = c("Both strong", "Strong and None", "NA and strong", "Both none", "NA and none"))
g.sumparahomos <- ggplot(thinsumparahomos, aes(homodimers, count, fill = inttype)) + geom_bar(position = "dodge", stat = 'Identity') + scale_fill_manual(values = c("#0090F3","#FE3131", "#C300EC"), labels = c("None", "Strong", "Weak"), name = 'Interaction\nstrength') +
  labs(x = "Homodimerization group", y = "Count", title = '201006 EA66k Paralog interactions by homodimerization' )
g.sumparahomos
ggsave(g.sumparahomos, file = '../figures/201006 Paralog interactions by homodimertype.pdf')



sumparahomosparents <- group_by(paraintsbyhomos, as.factor(combohomo), as.factor(parentcombohomo)) %>%
  summarise('strongints' = sum(intscore == 1.0), 'weakints' = sum(intscore == 0.5), 'nonints' = sum(intscore == 0.0), 'notmeasured' = sum(is.na(intscore) | intscore == -1.0))
colnames(sumparahomosparents) <- c('childhomo', 'parenthomo', 'strongints','weakints','nonints','notmeasured')
sumparahomosparents <- ungroup(sumparahomosparents)
sumparahomosparents <- gather(sumparahomosparents, inttype, counts, -childhomo, -parenthomo)
sumparahomosparents$inttype <- factor(sumparahomosparents$inttype, levels = c('strongints','weakints','nonints','notmeasured'))

intlabs <- c('nonints' = 'None','notmeasured' = 'NA', 'strongints' = 'Strong' , 'weakints' = 'Weak')
g.homointtrans <- ggplot(sumparahomosparents, aes(as.factor(childhomo), as.factor(parenthomo))) + geom_tile(aes(fill = counts)) + facet_wrap(~inttype, labeller = labeller(inttype = intlabs)) +
  labs(x = 'Child paralog homodimerization', y = 'Parent paralog homodimerization', title = '201007 Paralog interaction scores by parent/child homodimerization') + scale_fill_viridis(name = "Count") +
  scale_x_discrete(labels = c("NA and\nnone", "Both\nnone", "NA and\nstrong", "Strong\nand none", "Both\nstrong")) + scale_y_discrete(labels = c("NA and\nnone", "Both\nnone", "NA and\nstrong", "Strong\nand none", "Both\nstrong")) +
  theme(strip.background = element_rect(fill = 'White'))
g.homointtrans

ggsave(g.homointtrans, file = '../figures/201007 Paralog interactions by parent and child homodimerization.pdf')


distvsspec <- data.frame(paraintsbyhomos, "isspecgain" = 0)
distvsspec$isspecgain[((distvsspec$xpep %in% anctodiv$xpep & distvsspec$ypep %in% anctodiv$ypep) | (distvsspec$ypep %in% anctodiv$xpep & distvsspec$xpep %in% anctodiv$ypep)) & distvsspec$intscore == 0] <- 1

distcomp <- t.test(filter(distvsspec, isspecgain == 1)$distance, filter(distvsspec, isspecgain == 0)$distance)
g.distvsspec <- ggplot(distvsspec, aes(as.factor(isspecgain), distance)) + geom_boxplot(size = 1.3) + geom_jitter(width = 0.2, alpha = 0.7, size = 2.5, color = "#0090F3") + labs(x = 'Is point at which specificity is gained', y = "Branch length", title = '201007 EA66k paralogs branch length for specgain') +
  scale_x_discrete(labels = c('No', 'Yes')) + geom_text(aes(x = 2, y = 5, label = signif(distcomp$p.value, digits = 3)))
g.distvsspec

ggsave(g.distvsspec, file = '../figures/201007 paralogs by branch length vs specgain.pdf')

g.distvscorrspecgain <- ggplot(filter(distvsspec, isspecgain == 1), aes(distance, pearsonlog)) + geom_point(color = "#0090F3", size = 4) + #scale_color_manual(values = c("#0090F3","#C300EC", "#FE3131"), labels = c("None and missing", "One","Both"), name = 'Homodimerization') +
  labs(x = "Branch length", y = "Correlation Coefficient", title = '201007 EA66k Spec gains branch length vs correlation') + geom_text(aes(label = paste(xpep, ypep)))
g.distvscorrspecgain
ggsave(g.distvscorrspecgain, file = '../figures/201007 specgains branchlength vs correlation.pdf')

specpoints <- left_join(disttointprof, distvsspec) %>%
  filter(isspecgain == 1) %>%
  select(colnames(disttointprof))
ancpoints <- data.frame()
for(row in 1:nrow(specpoints)){
  anc1 <- filter(disttointprof, xpep == specpoints[row,]$xpep, ypep == specpoints[row,]$ancestor)
  anc2 <- filter(disttointprof, ypep == specpoints[row,]$xpep, xpep == specpoints[row,]$ancestor)
  anc3 <- filter(disttointprof, xpep == specpoints[row,]$ancestor, ypep == specpoints[row,]$ypep)
  anc4 <- filter(disttointprof, ypep == specpoints[row,]$ancestor, xpep == specpoints[row,]$ypep)
  ancpoints <- bind_rows(ancpoints, anc1, anc2, anc3, anc4)
}

specpoints <- data.frame(specpoints, 'rel' = 'paralogs')
ancpoints <- data.frame(ancpoints, 'rel' = 'ancdesc')
ancsspecs <- bind_rows(specpoints, ancpoints)

g.distvscorrspecgainancs <- ggplot(ancsspecs, aes(distance, pearsonlog, color = rel)) + geom_point(size = 4) + scale_color_manual(values = c("#0090F3", "#FE3131"), name = 'Relationship') +
  labs(x = "Branch length", y = "Correlation Coefficient", title = '201007 EA66k Spec gains branch length vs correlation') + geom_text(aes(label = paste(xpep, ypep)))
g.distvscorrspecgainancs

g.distvscorrspecgainancs <- ggplot(ancsspecs, aes(rel, pearsonlog, fill = rel)) + geom_violin() + geom_boxplot(fill = 'Black', width = 0.03) + stat_summary(fun.y = median, geom = "point", color = "White", show.legend = FALSE)  + 
  scale_fill_manual(values = c("#0090F3", "#FE3131"),labels = c('Anc/Desc','Descendents'), name = 'Relationship') +
  labs(x = "Relationship", y = "Correlation Coefficient", title = '201007 EA66k correlation at spec gain of intprof') + scale_y_continuous(limits = c(-0.5,1.0)) + scale_x_discrete(labels = c('Ancestor/Descendant','Descendant')) +
  theme(text = element_text(size = 14,family = 'Myriad Web Pro'), axis.title = element_text(size = 18))
g.distvscorrspecgainancs

ggsave(g.distvscorrspecgainancs, file = '../figures/201029 corr at specgain ancdesc vs desc.pdf', scale = 0.6)

sumparahomosparents <- group_by(paraintsbyhomos, as.factor(combohomo), as.factor(parentcombohomo)) %>%
  summarise('strcstrp' = sum(intscore == 1.0 & parentscore == 1.0), 'strcweakp' = sum(intscore == 1.0 & parentscore == 0.5), 'strcnonp' = sum(intscore == 1.0 & parentscore == 0.0), 'weakcstrp' = sum(intscore == 0.5 & parentscore == 1.0),
            'weakcweakp' = sum(intscore == 0.5 & parentscore == 0.5), 'weakcnonp' = sum(intscore == 0.5 & parentscore == 0.0), 'noncstrp' = sum(intscore == 0.0 & parentscore == 1.0), 'noncweakp' = sum(intscore == 0.0 & parentscore == 0.5), 'noncnonp' = sum(intscore == 0.0 & parentscore == 0.0))
colnames(sumparahomosparents)[1:2] <- c('childhomo', 'parenthomo')
sumparahomosparents <- ungroup(sumparahomosparents)
sumparahomosparents <- gather(sumparahomosparents, inttype, counts, -childhomo, -parenthomo)
#sumparahomosparents <- gather(sumparahomosparents, partype, parcounts, -childhomo, -parenthomo, -inttype, -counts)
sumparahomosparents <- mutate(sumparahomosparents, 'childints' = paste(as.character(childhomo),inttype, sep = ''), "parentints" = paste(as.character(parenthomo), inttype, sep = ''))

sumparahomosparents$inttype <- factor(sumparahomosparents$inttype, levels = c('strongints','weakints','nonints','notmeasured'))

g.homointtrans <- ggplot(sumparahomosparents, aes(as.factor(childints), as.factor(parentints))) + geom_tile(aes(fill = counts)) +
  labs(x = 'Child paralog homodimerization', y = 'Parent paralog homodimerization', title = '201007 Paralog interaction scores by parent/child homodimerization') + scale_fill_viridis(name = "Count") +
  theme(axis.text.x = element_text(angle = 90)) #+
  scale_x_discrete(labels = c("NA and\nnone", "Both\nnone", "NA and\nstrong", "Strong\nand none", "Both\nstrong")) + scale_y_discrete(labels = c("NA and\nnone", "Both\nnone", "NA and\nstrong", "Strong\nand none", "Both\nstrong")) +
  theme(strip.background = element_rect(fill = 'White'))
g.homointtrans

sumparahomosparents$inttype <- factor(sumparahomosparents$inttype, levels = c('strcstrp', 'strcweakp', 'strcnonp', 'weakcstrp', 'weakcweakp', 'weakcnonp', 'noncstrp', 'noncweakp','noncnonp'))

intlabs <- c('strcstrp' = 'Strong/Strong', 'strcweakp' = 'Strong/Weak', 'strcnonp' = 'Strong/None', 'weakcstrp'='Weak/Strong', 'weakcweakp' = 'Weak/Weak', 'weakcnonp' = 'Weak/None', 'noncstrp' = 'None/Strong', 'noncweakp' = 'None/Weak','noncnonp' = 'None/None')
g.maybe  <- ggplot(sumparahomosparents, aes(as.factor(childhomo), as.factor(parenthomo))) + geom_tile(aes(fill = counts)) + facet_wrap(~inttype, labeller = labeller(inttype =intlabs)) +
  scale_fill_viridis(name = 'Count') + scale_x_discrete(labels = c("NA and\nnone", "Both\nnone", "NA and\nstrong", "Strong\nand none", "Both\nstrong")) + scale_y_discrete(labels = c("NA and\nnone", "Both\nnone", "NA and\nstrong", "Strong\nand none", "Both\nstrong")) +
  labs(x = 'Child paralog homodimerization', y = 'Parent paralog homodimerization', title = '201007 Paralog interaction scores by child/parent homodimerization')
g.maybe
ggsave(g.maybe, file = '../figures/201009 paralog interactions by parent child homodimer interactions.pdf')



p174 <- filter(disttointprof, xpep == 174 | ypep == 174)
p255 <- filter(disttointprof, xpep == 255 | ypep == 255)
p185 <- filter(disttointprof, xpep == 185 | ypep == 185)

p174 <- data.frame(p174, 'age' = 'ancestor')
p255 <- data.frame(p255, 'age' = 'leg1')
p185 <- data.frame(p185, 'age' = 'leg2')

firstspecgain <- bind_rows(p174, p255, p185) 
g.firstspectgain <- ggplot(firstspecgain, aes(age, fill = as.factor(intscore))) + geom_bar(position = 'fill') +scale_fill_manual(values = c("#0090F3","#C300EC","#FE3131"), name  = "Interaction\nstrength", labels = c("None","Weak","Strong")) +
  labs(x = "Ancestor/Descendant", y = "Percentage", title = "201016 EA66k Interaction types by position in specificity gain") + scale_x_discrete(labels = c("Ancestor","Paralog 1","Paralog 2"))
g.firstspectgain

ancs <- filter(disttointprof, xpep %in% filter(distvsspec, isspecgain == 1)$ancestor | ypep %in% filter(distvsspec, isspecgain == 1)$ancestor)
desc1 <- filter(disttointprof, xpep %in% filter(distvsspec, isspecgain == 1)$xpep | ypep %in% filter(distvsspec, isspecgain == 1)$xpep)
desc2 <- filter(disttointprof, xpep %in% filter(distvsspec, isspecgain == 1)$ypep | ypep %in% filter(distvsspec, isspecgain == 1)$ypep)

ancs <- data.frame(ancs, 'age' = 'ancestor')
desc1 <- data.frame(desc1, 'age' = 'leg1')
desc2 <- data.frame(desc2, 'age' = 'leg1')

allspecsgain <- bind_rows(ancs, desc1, desc2)
allspecsgain <- data.frame(allspecsgain, 'whichpar' = NA)
allspecsgain$whichpar[allspecsgain$xpep %in% c(173, 224, 278) | allspecsgain$ypep %in% c(173, 224, 278)] <- 173
allspecsgain$whichpar[allspecsgain$xpep %in% c(174,255, 185) | allspecsgain$ypep %in% c(174,255, 185)] <- 174
allspecsgain$whichpar[allspecsgain$xpep %in% c(175, 'hlfcrassostrea_gigas1', 'crassostrea_gigas1')  | allspecsgain$ypep %in% c(175, 'hlfcrassostrea_gigas1', 'crassostrea_gigas1') ] <- 175
allspecsgain$whichpar[allspecsgain$xpep %in% c(200, 218, 210) | allspecsgain$ypep %in% c(200, 218, 210)] <- 200
allspecsgain$whichpar[allspecsgain$xpep %in% c(239, 242, 246) | allspecsgain$ypep %in% c(239, 242, 246)] <- 239
allspecsgain$whichpar[allspecsgain$xpep %in% c(310, 331, 312) | allspecsgain$ypep %in% c(310, 331, 312)] <- 310
allspecsgain$whichpar[allspecsgain$xpep %in% c(317, 322, 327) | allspecsgain$ypep %in% c(317, 322, 327)] <- 317
allspecsgain$whichpar[allspecsgain$xpep %in% c(323, 333, 331) | allspecsgain$ypep %in% c(323, 333, 331)] <- 323


g.allspecsints <- ggplot(allspecsgain, aes(age, fill = as.factor(intscore))) + geom_bar(position = 'stack') +scale_fill_manual(values = c("#0090F3","#C300EC","#FE3131"), name  = "Interaction\nstrength", labels = c("None","Weak","Strong")) +
  labs(x = "Ancestor/Descendant", y = "Percentage", title = "201016 EA66k Interaction types by position in specificity gain") + scale_x_discrete(labels = c("Ancestor","Paralog 1","Paralog 2")) + facet_wrap(~whichpar)
g.allspecsints

sumallspecs <- group_by(allspecsgain, age, whichpar) %>%
  summarise('nonint' = sum(intscore == 0)/n(), 'weakint' = sum(intscore == 0.5)/n(), 'strongint' = sum(intscore == 1.0)/n())

sumallspecs <- gather(sumallspecs, 'inttype', 'counts', -age, -whichpar)

sumsumallspecs <- group_by(sumallspecs, age, inttype) %>%
  summarise(intavg = mean(counts), intsd = sd(counts))
sumsumallspecs$inttype <- factor(sumsumallspecs$inttype, levels=  c('nonint','weakint','strongint'))


intavgsigs <- data.frame(x=c(0.2, 0.4, 0.6), y = c(0.9, 0.8,0.7), label = c(round(t.test(filter(sumallspecs, age == 'ancestor', inttype == "nonint")$counts, filter(sumallspecs, age == 'leg1', inttype == "nonint")$counts)$p.value, digits = 3),
                                                          round(t.test(filter(sumallspecs, age == 'ancestor', inttype == "weakint")$counts, filter(sumallspecs, age == 'leg1', inttype == "weakint")$counts)$p.value, digits = 3),
                                                          round(t.test(filter(sumallspecs, age == 'ancestor', inttype == "strongint")$counts, filter(sumallspecs, age == 'leg1', inttype == "strongint")$counts)$p.value, digits = 3)))
intavgerr <- aes(ymin = intavg - intsd, ymax = intavg + intsd)
g.sumsumallspecs <- ggplot(sumsumallspecs, aes(age, intavg, fill = inttype)) + geom_bar(stat = 'identity', position = 'dodge') +scale_fill_manual(values = c("#0090F3","#C300EC","#FE3131"), name  = "Interaction\nstrength", labels = c("None","Weak","Strong")) +
  geom_errorbar(intavgerr, position = position_dodge(width = 0.9), width = 0.3) + labs(x = '', y = 'Percentage of interactions', title = '201016 EA66k percentage of anc/desc all interactions for specgain by intype') + scale_x_discrete(labels = c("Ancestor", "Descendant")) #+
  geom_text(data = intavgsigs, aes(x=x, y=y, label=label, color = 'Black'))
g.sumsumallspecs
ggsave(g.sumsumallspecs, file= '../figures/201028 Interactions percentages by ancdes.pdf')


ancs <- filter(disttointprof, xpep %in% filter(distvsspec, isspecgain == 1)$ancestor | ypep %in% filter(distvsspec, isspecgain == 1)$ancestor)
desc1 <- filter(disttointprof, xpep %in% filter(distvsspec, isspecgain == 1)$xpep | ypep %in% filter(distvsspec, isspecgain == 1)$xpep)
desc2 <- filter(disttointprof, xpep %in% filter(distvsspec, isspecgain == 1)$ypep | ypep %in% filter(distvsspec, isspecgain == 1)$ypep)

ancs <- data.frame(ancs, 'age' = 'ancestor')
desc1 <- data.frame(desc1, 'age' = 'leg1')
desc2 <- data.frame(desc2, 'age' = 'leg1')



tree <- read.tree(file = "/Users/Cliff/Documents/Kosuri Lab/Design docs/181114 Ancestral tree/190918_EA_full_tree1.txt")
tree <- as_tibble(tree$TREE1_with_branch_lengths)
tree$label[is.na(tree$label)] <- as.character(tree$node[is.na(tree$label)])
tree$label <- tolower(tree$label)

p198 <- filter(disttointprof, xpep == '291')
p198y <- filter(disttointprof, ypep == '291') %>%
  select(ypep, xpep, medRD, distance, corr, logtype, pearsonlog, spearman, intscore, parentscore, paraparents, ancestor, isparentchild)
colnames(p198y) <- c('ypep', 'xpep', 'medRD', 'distance', 'corr', 'logtype', 'pearsonlog', 'spearman', 'intscore', 'parentscore', 'paraparents', 'ancestor', 'isparentchild')
p198 <- bind_rows(p198, p198y)
p198 <- distinct(p198)

tree198 <- left_join(tree, p198, by = c("label" = "ypep"))
tree198 <- as.treedata(tree198)

g.p198 <- ggtree(tree198) + labs(title = "191104 protein 198 Interaction Tree") + geom_point(aes(color = as.factor(intscore)),size = 5) + geom_nodelab(vjust = -1, size = 3)+ geom_tiplab(size = 3) +
  scale_color_manual(values = c("#0090F3","#C300EC","#FE3131"), na.value = '#A7A7A7', name  = "Interaction\nstrength", labels = c("None","Weak","Strong",'Not measured')) + theme(legend.position = "right")
g.p198


p176 <- filter(disttointprof, xpep %in% c('194','192','176'))
tree176 <- left_join(tree, p176, by = c("label" = "ypep"))
tree176 <- as.treedata(tree176)

g.p176 <- ggtree(tree176) + labs(title = "191104 protein 176 Interaction Tree") + geom_point(aes(color = as.factor(intscore)),size = 5) + geom_nodelab(vjust = -1, size = 3)+ geom_tiplab(size = 3) +
  scale_color_manual(values = c("#0090F3","#C300EC","#FE3131"), na.value = '#A7A7A7', name  = "Interaction\nstrength", labels = c("None","Weak","Strong",'Not measured')) + theme(legend.position = "right")
g.p176

parentchilds <- filter(disttointprof, isparentchild == 'Yes', pearsonlog != 1)
g.parenchildcounts <- ggplot(filter(parentchilds, pearsonlog < .9999), aes(isparentchild, fill = as.factor(intscore))) + geom_bar(position = 'fill') +
  scale_fill_manual(values = c("#0090F3","#C300EC","#FE3131"), labels = c('None', 'Weak','Strong'), name = 'Interaction\nstrength') +
  labs(x='Parent/Child Interactions', y = 'Percentage', title = "201029 EA66k Parent/child interactions") + scale_x_discrete(label = '') + 
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 16))
g.parenchildcounts
ggsave(g.parenchildcounts, file = '../figures/201029 Parent and child interactions percents.pdf', scale = 0.6)

g.parentchilds <- ggplot(filter(parentchilds, pearsonlog < .9999), aes(pearsonlog, fill = as.factor(intscore))) + geom_histogram() +
  scale_fill_manual(values = c("#0090F3","#C300EC","#FE3131"), labels = c('None', 'Weak','Strong'), name = 'Interaction\nstrength') +
  labs(x = 'Correlation between parent and child', y = 'Count', title = '201031 EA66k Parent/Child Interactions by correlation') +
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 16))
g.parentchilds

ggsave(g.parentchilds, file = '../figures/201031 Parentchild correlation histogram.pdf', scale = 0.6)

g.parentchildcorrbl <- ggplot(filter(parentchilds, pearsonlog < .9999), aes(distance, pearsonlog, color = as.factor(intscore))) + geom_point(alpha = 0.7) + 
  scale_color_manual(values = c("#0090F3","#C300EC","#FE3131"), labels = c('None', 'Weak','Strong'), name = 'Interaction\nstrength') +
  labs(x = 'Branch Length', y = 'Correlation between parent and child', title = '201031 EA66k Parent/Child interactions by correlation, by distance') + 
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 16))
g.parentchildcorrbl

ggsave(g.parentchildcorrbl, file = '../figures/201031 Parentchild correlation by brlen.pdf', scale = 0.6)

parentchildswmatched <- filter(disttointprof, isparentchild == 'Yes')

g.distvsintscoreparent <- ggplot(parentchildswmatched, aes(distance, fill = as.factor(intscore))) + geom_histogram() +
  scale_fill_manual(values = c('#A7A7A7', "#0090F3","#C300EC","#FE3131"), labels = c('Not measured','None', 'Weak','Strong'), name = 'Interaction\nstrength')
g.distvsintscoreparent

