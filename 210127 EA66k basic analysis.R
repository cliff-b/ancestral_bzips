################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(viridis)
library(wesanderson)
library(GGally)
library(MASS)
library(extrafont)
library(ggtree)
library(tidytree)
library(dplyr)
library(tidyr)
library(cowplot)
################################################# Read in the data and mash it together #################################################
setwd("/Users/Cliff/Documents/Kosuri Lab/R and python/ancestral tree/data_for_R")

numpars <- read.csv("200929_num_paralogs_till_spec.csv")
numparsweak <- read.csv("201208_num_paralogs_till_spec_with_weak.csv")
anctodiv <- read.csv("200722_anc_to_diverge_ints.csv")
anctodivweak <- read.csv("201208_anc_to_diverge_with_weak.csv")
allpars <- read.csv("200722_allpar_score_and_dist.csv")
homos <- read.csv('200930_homodimer_counting.csv')
disttointprof <- read.csv('200930_interactionprofiles_cors.csv')

numpars <- data.frame(numpars, 'strength' = 'strong')
numparsweak <- data.frame(numparsweak, 'strength' = 'weak')
numpars <- bind_rows(numpars, numparsweak)
numpars$strength <- factor(numpars$strength, levels = c('weak', 'strong'))
numpars <- mutate(numpars, 'morethanone' = (X.numpars.. > 1))

g.paralogsbeforespec <- ggplot(numpars, aes(morethanone, fill = strength)) + geom_bar(color = "#A7A7A7") + labs(x = "Number of paralogs before specificity is gained", y = "Count", title = "200722 Path to specificity starts quickly but is not finished fast") +
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18))+ scale_x_discrete(labels=c('1','>1'))  +scale_y_continuous(breaks = c(0, 2, 4,6, 8)) + scale_fill_manual(values = c("#FE3131","#0091F4"), name = 'Type of\nloss', labels = c('Weak', 'None'))
g.paralogsbeforespec

#ggsave(g.paralogsbeforespec, file = '../figures/201215_num_paralogs_till_specificity for weak and strong.pdf', scale = 0.6)
# ggsave(g.paralogsbeforespec, file = "200722_num_paralogs_till_specificity.pdf", scale = 0.6)
# 
# 

g.specgains <- ggplot(numpars, aes(X.distance., fill = strength)) + geom_histogram(bins = 8, color = '#A7A7A7') + labs(x = "Branch Length", y = "Count", title = "201001 EA66k distance needed to gain specificty") +
  scale_x_continuous(limits = c(0, 1.4)) + scale_fill_manual(values = c("#FE3131","#0091F4"), name = 'Type of\ninteraction', labels = c('Weak', 'Strong'))
g.specgains

#ggsave(g.specgains, file = '../figures/201001 EA66k distance needed to gain specificty.pdf', scale = 0.6)
#ggsave(g.specgains, file = '../figures/201214 EA66k distance needed to gain spec for weak and strong.pdf', scale = 0.6)

colnames(homos) <- c("xpep", "ypep","intscore","extant")
numpars <- left_join(numpars, homos, by = 'xpep')
numpars <- left_join(numpars, homos, by = c('ypep.x' = 'ypep'))
numpars$intscore.y[numpars$intscore.y == -1] <- -0.5
numpars <- mutate(numpars, 'combohomo' = intscore.x + intscore.y)

g.spechomos <- ggplot(numpars, aes(as.factor(combohomo), fill = strength)) + geom_bar(color = "#A7A7A7") + labs(x = "Homodimers at point specificity gained", y = "Count", title = "200722 Often one protein becomes a weaker interactor") +
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18))+ scale_x_discrete(labels = c('No and NA', 'Neither', 'One', 'Both'))  +scale_y_continuous(breaks = c(0, 2, 4,6, 8)) + scale_fill_manual(values = c("#FE3131","#0091F4"), name = 'Type of\nloss', labels = c('Weak', 'None'))
g.spechomos
#ggsave(g.spechomos, file = '../figures/201215 EA66k homodimerzation at specgain.pdf', scale = 0.6)

numpars <- select(numpars, 'xpep'=xpep.x, 'ypep' = ypep.x, 'ancestor' = X..ancestor., 'distance' = X.distance., 'numpars' = X.numpars.., strength, intscore.x, intscore.y, combohomo)

allintswithspec <- data.frame()
for (i in 1:nrow(numpars)) {
  firstor <- filter(disttointprof, xpep == numpars[i,]$xpep)
  secondor <- filter(disttointprof, ypep == numpars[i,]$xpep)
  thirdor <- filter(disttointprof, xpep == numpars[i,]$ypep)
  fourthor <- filter(disttointprof, ypep == numpars[i,]$ypep)
  allors <- bind_rows(firstor, secondor, thirdor, fourthor)
  allors <- data.frame(allors, 'strength' = numpars[i,]$strength, 'combohomo'= numpars[i,]$combohomo)
  allintswithspec <- bind_rows(allintswithspec, allors)
}

sumallintswithspec <- group_by(allintswithspec, strength, combohomo) %>%
  summarise('nonint' = sum(intscore == 0)/n(), 'weakint' = sum(intscore == 0.5)/n(), 'strongint' = sum(intscore == 1.0)/n())

sumallintswithspec <- gather(sumallintswithspec, key='inttype', value='percent', -strength, -combohomo)
sumallintswithspec$inttype <- factor(sumallintswithspec$inttype, levels = c('strongint','weakint','nonint'))

g.sumallintswithspec <- ggplot(sumallintswithspec, aes(as.factor(combohomo), percent, fill = inttype)) + geom_bar(stat='identity', position = 'fill') +
  scale_fill_manual(values = c( "#FE3131","#C300EC","#0090F3"), name = "Interaction\nstrength", labels = c( "Strong", "Weak", "None")) + 
  scale_x_discrete(labels = c(c('No and NA', 'Neither', 'One', 'Both')))  + labs(x= "Homodimers at point specificity gained", y = 'Percent', title = '201216 EA66k all interactions at time specficity is gained') +
  theme(text = element_text(size = 14, family = 'Myriad Web Pro'), axis.title = element_text(size = 18))+ 
  facet_wrap(~strength)
g.sumallintswithspec

#ggsave(g.sumallintswithspec, filename = '../figures/201216 EA66k interactions at specgain by homotype and strength.pdf', scale = 0.6)


colnames(allpars) <- c("xpep", "ypep","intscore","distance")



#allpars$distance <- round(allpars$distance, digits = 1)
g.allparsinthist <- ggplot(allpars, aes(as.factor(intscore), fill = as.factor(intscore))) + geom_bar(position = "stack", show.legend = F) + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131")) +
  scale_x_discrete(labels = c("Not measured", "None", "Weak", "Strong")) + labs(x = "Interaction type", y = "Count", title = "200930 EA66k All paralogs by interaction type counts")
g.allparsinthist
#ggsave(g.allparsinthist, file = "../figures/200930 Paralogs by interaction type.pdf", scale = 0.6)


g.homosinthist <- ggplot(homos, aes(as.factor(intscore), fill = as.factor(intscore))) + geom_bar(position = "stack", show.legend = F) + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#FE3131")) +
  scale_x_discrete(labels = c("Not measured", "None", "Strong")) + labs(x = "Interaction type", y = "Count", title = "200930 EA66k All Homodimers by interaction type counts")
g.homosinthist
#ggsave(g.homosinthist, file = "../figures/200930 homodimers by interaction type.pdf", scale = 0.6)

allpars <- data.frame(allpars, "pairtype" = "para")
allpars <- select(allpars, -distance)
homos <- data.frame(homos, "pairtype" = "homo")
homos <- select(homos, -extant)
parshomo <- bind_rows(allpars, homos)

g.hpinthist <- ggplot(parshomo, aes(pairtype, fill = as.factor(intscore))) + geom_bar(position = "stack") + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("Not measured", "None", "Weak", "Strong")) +
  scale_x_discrete(labels = c("Homodimers","Heterodimers")) + labs(x = "Interaction type", y = "Count", title = "200930 EA66k Homodimers and paralogs by interaction type counts")
g.hpinthist
#ggsave(g.hpinthist, file = "../figures/200930 homodimers and paralogs by interaction type.pdf", scale = 0.6)

g.hpinthistper <- ggplot(parshomo, aes(pairtype, fill = as.factor(intscore))) + geom_bar(position = "fill") + scale_fill_manual(values = c( "#A7A7A7","#0090F3","#C300EC","#FE3131"), name = "Interaction\nstrength", labels = c("Not measured", "None", "Weak", "Strong")) +
  scale_x_discrete(labels = c("Homodimers","Heterodimers")) + labs(x = "Interaction type", y = "Percentage", title = "200930 EA66k Homodimers and paralogs by interaction type percentages")
g.hpinthistper

#ggsave(g.hpinthistper, file = "../figures/200930 homodimers by interaction type percentage.pdf", scale = 0.6)

colnames(anctodiv) <- c('xpep', 'ypep','distance','anc_to_div')
anctodiv <- separate(anctodiv, anc_to_div, into = c("xpep_anc_int","ypep_anc_int"), sep = "\\+")
anctodiv <- data.frame(anctodiv, "ancints" = as.character(c("")))
anctodiv$ancints <- as.character(anctodiv$ancints)
anctodiv$xpep_anc_int <- as.numeric(anctodiv$xpep_anc_int)
anctodiv$ypep_anc_int <- as.numeric(anctodiv$ypep_anc_int)
anctodiv$ancints[anctodiv$xpep_anc_int == 1.0 & anctodiv$ypep_anc_int == 1.0] <- "Both"
anctodiv$ancints[anctodiv$xpep_anc_int == 1.0 & anctodiv$ypep_anc_int == 0.5] <- "Both"
anctodiv$ancints[anctodiv$xpep_anc_int == 0.5 & anctodiv$ypep_anc_int == 1.0] <- "Both"
anctodiv$ancints[anctodiv$xpep_anc_int == 1.0 & anctodiv$ypep_anc_int == 0.0] <- "one"
anctodiv$ancints[anctodiv$xpep_anc_int == 0.0 & anctodiv$ypep_anc_int == 1.0] <- "one"

##########if taking 0.5 to be noninteracting ##############
colnames(anctodivweak) <- c('xpep', 'ypep','distance','anc_to_div')
anctodivweak <- separate(anctodivweak, anc_to_div, into = c("xpep_anc_int","ypep_anc_int"), sep = "\\+")
anctodivweak <- data.frame(anctodivweak, "ancints" = as.character(c("")))
anctodivweak$ancints <- as.character(anctodivweak$ancints)
anctodivweak$xpep_anc_int <- as.numeric(anctodivweak$xpep_anc_int)
anctodivweak$ypep_anc_int <- as.numeric(anctodivweak$ypep_anc_int)
anctodivweak$ancints[anctodivweak$xpep_anc_int == 1.0 & anctodivweak$ypep_anc_int == 1.0] <- "Both"
anctodivweak$ancints[anctodivweak$xpep_anc_int == 1.0 & anctodivweak$ypep_anc_int == 0.5] <- "one"
anctodivweak$ancints[anctodivweak$xpep_anc_int == 1.0 & anctodivweak$ypep_anc_int == 0.0] <- "one"
anctodivweak$ancints[anctodivweak$xpep_anc_int == 0.0 & anctodivweak$ypep_anc_int == 1.0] <- "one"
anctodivweak$ancints[anctodivweak$xpep_anc_int == 0.5 & anctodivweak$ypep_anc_int == 1.0] <- "one"
anctodivweak$ancints[anctodivweak$xpep_anc_int == 0.5 & anctodivweak$ypep_anc_int == 0.5] <- "none"

anctodiv <- data.frame(anctodiv, 'strength' = 'strong')
anctodivweak <- data.frame(anctodivweak, 'strength' = 'weak')
anctodiv <- bind_rows(anctodiv, anctodivweak)


anctodiv2 <- left_join(anctodiv, disttointprof, by = c("xpep", "ypep"))
anctodiv2 <- left_join(anctodiv2, disttointprof, by = c('xpep'='ypep','ypep'='xpep'))
anctodiv2 <- mutate(anctodiv2, "pearsonlog" = coalesce(pearsonlog.x, pearsonlog.y))
anctodiv2 <- mutate(anctodiv2, "ancestor" = coalesce(ancestor.x, ancestor.y))
anctodiv2 <- select(anctodiv2, xpep, ypep, 'distance' = distance.x, anctype, ancscore, ancints, pearsonlog, ancestor, strength)



anctodiv <- gather(anctodiv, 'anctype', 'ancscore', -xpep, -ypep, -distance, -ancints,-strength )
anctodiv$ancscore <- factor(anctodiv$ancscore, levels = c(1, 0.5, 0))
anctodiv$ancints <- factor(anctodiv$ancints, levels = c("Both", "one", "none"))


g.changesforspec <- ggplot(anctodiv, aes(ancints, fill = ancscore)) + geom_bar(position = "stack") + scale_fill_manual(values = c("#FE3131","#C300EC", "#0090F3", "#A7A7A7"), name = 'Interaction\ntype', labels = c('Strong','Weak','None')) +
  scale_x_discrete(labels = c("Both", "One", "None")) + labs(x = "Branches interacting with ancestor", y = "Count", title = "200930 EA66k Changes on which branches led to specificity") +
  facet_wrap(~strength) + theme(strip.background = element_rect(fill = 'White'), text = element_text(family="Myriad Web Pro", size = 14), axis.title = element_text(size = 18))
g.changesforspec
#ggsave(g.changesforspec, file = "../figures/201208 EA66k branch changes to spec for weak and strong loss.pdf", scale = 0.6 )
#ggsave(g.changesforspec, file = "../figures/200930 EA66k branch changes to spec.pdf", scale = 0.6)

g.changesspecdist <- ggplot(anctodiv, aes(ancints, distance, fill = ancints)) + geom_dotplot(binaxis = "y", show.legend = F)+ scale_fill_manual(values = c("#FE3131", "#0090F3")) +
  scale_x_discrete(labels = c("Both", "Single")) + labs(x = "How many branches specificity changes occured on", y = "Distance to gain specificity", title = "200930 EA66k Changes on which branches led to specificity")
g.changesspecdist
ggsave(g.changesspecdist, file = "../figures/200930 EA66k branch changes to spec by dist.pdf", scale = 0.6)


colnames(disttointprof) <- c('xpep', 'ypep', 'medRD', 'distance', 'corr','logtype','pearsonlog','spearman','intscore','parentscore','paraparents','ancestor','isparentchild','nodesancdesc')
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
distvsspec$isspecgain[((distvsspec$xpep %in% anctodiv$xpep & distvsspec$ypep %in% anctodiv$ypep) | (distvsspec$ypep %in% anctodiv$xpep & distvsspec$xpep %in% anctodiv$ypep)) & (distvsspec$intscore == 0 |distvsspec$intscore == 0.5)] <- 1

distcomp <- t.test(filter(distvsspec, isspecgain == 1)$distance, filter(distvsspec, isspecgain == 0)$distance)
g.distvsspec <- ggplot(distvsspec, aes(as.factor(isspecgain), distance)) + geom_boxplot(size = 1.3) + geom_jitter(width = 0.2, alpha = 0.7, size = 2.5, color = "#0090F3") + labs(x = 'Is point at which specificity is gained', y = "Branch length", title = '201007 EA66k paralogs branch length for specgain') +
  scale_x_discrete(labels = c('No', 'Yes')) + geom_text(aes(x = 2, y = 5, label = signif(distcomp$p.value, digits = 3)))
g.distvsspec

ggsave(g.distvsspec, file = '../figures/201007 paralogs by branch length vs specgain.pdf')

g.distvscorrspecgain <- ggplot(filter(distvsspec, isspecgain == 1), aes(distance, pearsonlog)) + geom_point(color = "#0090F3", size = 4) + #scale_color_manual(values = c("#0090F3","#C300EC", "#FE3131"), labels = c("None and missing", "One","Both"), name = 'Homodimerization') +
  labs(x = "Branch length", y = "Correlation Coefficient", title = '201007 EA66k Spec gains branch length vs correlation') + geom_text(aes(label = paste(xpep, ypep)))
g.distvscorrspecgain
ggsave(g.distvscorrspecgain, file = '../figures/201007 specgains branchlength vs correlation.pdf')

anctodivsames <- group_by(anctodiv2, ancestor) %>%
  summarise(pmax = max(pearsonlog), pmin = min(pearsonlog), dmax = max(distance), dmin = min(distance))

g.distvscorrspecgain <- ggplot(anctodiv2, aes(distance, pearsonlog, color = strength)) + geom_segment(data=anctodivsames, aes(x=dmin, y=pmax, xend=dmax, yend=pmin), color='Black')+ geom_point(size = 4) + scale_color_manual(values = c("#0090F3", "#FE3131"), labels = c("Strong", "Weak"), name = 'Type of\ninteraction') +
  labs(x = "Branch length", y = "Correlation Coefficient", title = '201007 EA66k Spec gains branch length vs correlation') + #+ geom_text(aes(label = paste(xpep, ypep)))
  theme(text = element_text(family="Myriad Web Pro", size = 14), axis.title = element_text(size=18))
g.distvscorrspecgain
ggsave(g.distvscorrspecgain, file = '../figures/201214 specgains branchlength vs correlation for strong and weak losses.pdf', scale=0.6)


specpoints <- anctodiv2
  #left_join(disttointprof, distvsspec) %>%
  #filter(isspecgain == 1) %>%
  #select(colnames(disttointprof))
ancpoints <- data.frame()
for(row in 1:nrow(specpoints)){
  anc1 <- filter(disttointprof, xpep == specpoints[row,]$xpep, ypep == specpoints[row,]$ancestor)
  anc2 <- filter(disttointprof, ypep == specpoints[row,]$xpep, xpep == specpoints[row,]$ancestor)
  anc3 <- filter(disttointprof, xpep == specpoints[row,]$ancestor, ypep == specpoints[row,]$ypep)
  anc4 <- filter(disttointprof, ypep == specpoints[row,]$ancestor, xpep == specpoints[row,]$ypep)
  ancgen <- bind_rows(anc1, anc2, anc3, anc4)
  ancgen <- data.frame(ancgen, 'strength' = specpoints[row,]$strength)
  ancpoints <- bind_rows(ancpoints, ancgen)
}

specpoints <- data.frame(specpoints, 'rel' = 'paralogs')
ancpoints <- data.frame(ancpoints, 'rel' = 'ancdesc')
ancsspecs <- bind_rows(specpoints, ancpoints)

g.distvscorrspecgainancs <- ggplot(ancsspecs, aes(distance, pearsonlog, color = rel)) + geom_point(size = 4) + scale_color_manual(values = c("#0090F3", "#FE3131"), name = 'Relationship') +
  labs(x = "Branch length", y = "Correlation Coefficient", title = '201007 EA66k Spec gains branch length vs correlation') + geom_text(aes(label = paste(xpep, ypep)))
g.distvscorrspecgainancs

g.distvscorrspecgainancs <- ggplot(ancsspecs, aes(rel, pearsonlog, fill = rel)) + geom_violin(show.legend=F) + geom_boxplot(fill = 'Black', width = 0.03, outlier.colour = NA) + stat_summary(fun.y = median, geom = "point", color = "White", show.legend = FALSE)  + 
  scale_fill_manual(values = c("#0090F3", "#FE3131"),labels = c('Anc/Desc','Descendents'), name = 'Relationship') +
  labs(x = "Relationship", y = "Correlation Coefficient", title = '201007 EA66k correlation at spec gain of intprof') + scale_y_continuous(limits = c(-0.5,1.0)) + scale_x_discrete(labels = c('Ancestor/Descendant','Descendant/Descendant')) +
  theme(text = element_text(size = 14,family = 'Myriad Web Pro'), axis.title = element_text(size = 18), strip.background = element_rect(fill='White')) + facet_wrap(~strength)
g.distvscorrspecgainancs

ggsave(g.distvscorrspecgainancs, file = '../figures/201029 corr at specgain ancdesc vs desc for weak and strong loss.pdf', scale = 0.6)
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

extantparas <- filter(intsbyhomos, homoextant.x == 'Yes', logtype == 'paralogs')
g.extantparas <- ggplot(extantparas, aes(as.factor(combohomo),fill =as.factor(intscore))) + geom_bar(position = 'fill') + scale_fill_manual(values = c("#0090F3","#C300EC","#FE3131"), labels = c('None', 'Weak','Strong'), name = 'Interaction\nstrength') +
  labs(x = 'Homodimerization of paralogs', y = 'Percent', title = '201106 EA66k Interactions of extant paralogs by homodimerization') + scale_x_discrete(labels = c('NA/None','None/None','NA/Strong','None/Strong','Strong/Strong'))
g.extantparas
ggsave(g.extantparas, file = '../figures/201106 EA66k interactions of extant paralogs by homodimerization percent.pdf', scale = 0.6 )

paras <- filter(disttointprof, logtype == 'paralogs')
paras$ancestor <- as.character(paras$ancestor) 
paras <- group_by(paras, ancestor)

avgs <- data.frame('number' = c(), 'avg' = c())
for (i in 1:10000) {
  samp <- sample_n(paras, 1)
  avgs <- bind_rows(avgs, data.frame('number' = i, 'avg' = mean(samp$distance)))
}

g.parapermutdistall <- ggplot(avgs, aes(avg)) + geom_density(fill = "#0090F3") +
  labs(x = 'Average branch length', y = "Density", title = '201117 EA66k permutation test on paralogs distance for specgain vs ancestor') + 
  theme(text = element_text(family = 'Myriad Web Pro', size = 14), axis.title = element_text(size = 18)) + geom_vline(xintercept = mean(anctodiv$distance)) +
  annotate("text", x = 0.6, y = 8, label = paste("p = ", round(sum(avgs$avg < mean(anctodiv$distance))/10000, digits = 3)))
g.parapermutdistall

ggsave(g.parapermutdistall, file = '../figures/201117_permutation test of spec gain branches to all paralogs.pdf', scale = 0.6)


para2 <- left_join(paras, filter(anctodiv,strength=='strong'), by = c('xpep', 'ypep'))
para2 <- left_join(para2, filter(anctodiv,strength=='strong'), by = c('xpep' = 'ypep', 'ypep' = 'xpep'))

para2 <- filter(para2, !is.na(distance.y) | !is.na(distance))
paras3 <- filter(paras, ancestor %in% para2$ancestor)

avgs <- data.frame('number' = c(), 'avg' = c())
for (i in 1:10000) {
  samp <- sample_n(paras3, 1)
  avgs <- bind_rows(avgs, data.frame('number' = i, 'avg' = mean(samp$distance)))
}

para2 <- left_join(paras, filter(anctodiv,strength=='weak'), by = c('xpep', 'ypep'))
para2 <- left_join(para2, filter(anctodiv,strength=='weak'), by = c('xpep' = 'ypep', 'ypep' = 'xpep'))

para2 <- filter(para2, !is.na(distance.y) | !is.na(distance))
paras4 <- filter(paras, ancestor %in% para2$ancestor)

avgs2 <- data.frame('number' = c(), 'avg' = c())
for (i in 1:10000) {
  samp <- sample_n(paras4, 1)
  avgs2 <- bind_rows(avgs2, data.frame('number' = i, 'avg' = mean(samp$distance)))
}

avgs <- data.frame(avgs, 'strength'='strong')
avgs2 <- data.frame(avgs2, 'strength' = 'weak')
avgs <- bind_rows(avgs, avgs2)

anctodivsum <- group_by(anctodiv, strength) %>%
  summarise(meandist = mean(distance))
g.parapermutdist <- ggplot(avgs, aes(avg, fill=strength)) + geom_histogram(data=filter(avgs, strength=='weak'),bins = 70) + geom_histogram(data=filter(avgs, strength=='strong'),bins = 70) + 
  labs(x = 'Average branch length', y = "Count", title = '201117 EA66k permutation test on specgaining paralogs distance for specgain vs ancestor') + 
  theme(text = element_text(family = 'Myriad Web Pro', size = 14), axis.title = element_text(size = 18)) + geom_vline(data=filter(anctodivsum, strength =='strong'), aes(xintercept = meandist), color ="#0091F4", size=2, linetype='dashed') +
  geom_vline(data=filter(anctodivsum, strength =='weak'), aes(xintercept = meandist), color ="#FE3131", size=2, linetype='dashed') +
  annotate("text", x = 0.5, y = 600, label = paste("p = ", round(sum(avgs$avg < filter(anctodivsum, strength=='weak')$meandist)/10000, digits = 3))) + 
  annotate("text", x = 0.75, y = 300, label = paste("p = ", round(sum(avgs$avg <  filter(anctodivsum, strength=='strong')$meandist)/10000, digits = 3))) +
  scale_fill_manual(values = c("#0091F4","#FE3131"), name = 'Type of\ninteraction', labels = c('Strong','Weak'))
g.parapermutdist

ggsave(g.parapermutdist, file = '../figures/201208_permutation test of spec gain branches to specgain branch paralogs for weak and strong.pdf', scale = 0.6)
ggsave(g.parapermutdist, file = '../figures/201117_permutation test of spec gain branches to specgain branch paralogs.pdf', scale = 0.6)

parentchilds <- read.csv(file = '201118_children_toparent_ints.csv')
colnames(parentchilds) = c('xpep','ypep', 'parent','child1', 'child2')
parentchilds <- mutate(parentchilds, 'samechild' = child1 == child2, "sameparent1" = child1 == parent, 'sameparent2' = child2 == parent) 
parentchilds <- mutate(parentchilds, 'bothparent' = sameparent1 & sameparent2)

g.parentchild <- ggplot(filter(parentchilds, child1 != -1, child2 != -1), aes(as.factor(bothparent), fill = as.factor(samechild))) + geom_bar(position = "stack", stat = "count")
g.parentchild


changesfromparents <- filter(parentchilds, bothparent == FALSE, child1 != -1, child2 != -1, parent != -1)

avgs <- data.frame('number' = c(), 'samesintsavg' = c())
for (i in 1:10000){
  samp1 <- sample_n(disttointprof, 64)
  samp2 <- sample_n(disttointprof, 64)
  row.names(samp1) <- c(1:64)
  row.names(samp2) <- c(1:64)
  samps <- merge(samp1, samp2, by = "row.names")
  avgs <- bind_rows(avgs, data.frame('number' = i, 'samesintsavg' = mean(samps$intscore.x == samps$intscore.y)))
}


g.changesfromparentsdist <- ggplot(avgs, aes(samesintsavg)) + geom_histogram(fill = "#0090F3", bins = 200) +
  labs(x = 'Average branch length', y = "Count", title = '201117 EA66k permutation test on specgaining paralogs distance for specgain vs ancestor') + 
  theme(text = element_text(family = 'Myriad Web Pro', size = 14), axis.title = element_text(size = 18)) + geom_vline(xintercept = mean(changesfromparents$samechild)) +
  annotate("text", x = 0.6, y = 300, label = paste("p = ", round(sum(avgs$samesintsavg < mean(changesfromparents$samechild)/10000, digits = 3))))
g.changesfromparentsdist

stepsquant <- read.csv(file = '201120number_of_steps_to_root.csv')
stepsquant <- read.csv(file = '201208number_of_steps_to_root.csv')
colnames(stepsquant) <- c('xpep', 'ypep', 'steps', 'curparent')
stepsquant <- mutate(stepsquant, 'junk' = "x")
g.stepsquant <- ggplot(stepsquant, aes(junk, fill = as.factor(steps))) + geom_bar(position = 'fill') + labs(x = "Number of steps back to root", y = "Percent", title = "201120 Number of steps back to root") +
  scale_x_discrete(labels = "") + scale_fill_brewer(name = "Steps back\nto root", palette = "Set1") + theme(text = element_text(family= "Myriad Web Pro", size = 14), axis.title = element_text(size = 18))
g.stepsquant

ggsave(g.stepsquant, file = '../figures/201120 steps to root quanted.pdf', scale = 0.6)
ggsave(g.stepsquant, file = '../figures/201208 steps to root quanted with weak as loss.pdf', scale = 0.6)

orthovspara <- read.csv('201120_orthosandparas.csv')
orthovsparaweak <- read.csv('201208_orthosandparasweakintsasints.csv')
orthovspara <- data.frame(orthovspara, 'strength'='strong')
orthovsparaweak <- data.frame(orthovsparaweak, 'strength' = 'weak')
orthovspara <- bind_rows(orthovspara, orthovsparaweak)
#orthovspara <- read.csv('201208_orthosandparasweakintsasloss.csv')
sumorthosvspara <- group_by(orthovspara, bin_avg, logtype, strength) %>%
  summarise(meanscore = mean(avg_score), medscore = median(avg_score), sdscore = sd(avg_score), me = qnorm(0.975) * sd(avg_score)/min(count), counts = mean(count), sdcount =  sd(count))

orthoerr <- aes(ymin = meanscore - sdscore, ymax = meanscore + sdscore)
g.orthovspara <- ggplot(sumorthosvspara, aes(bin_avg, meanscore, color = logtype)) + geom_point(size = 4) + geom_errorbar(orthoerr, width=0.05) +
  labs(x = "Branch length", y = "Fraction interacting", title = "201120 EA66k Interaction scores vs brlen by logtype") + theme_cowplot(font_family = "Myriad Web Pro") +
  scale_color_manual(labels = c("Ortholog","Paralog"), name = "", values = c("#0090F3","#FE3131")) + theme(text = element_text(size = 14), axis.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(-0.3, 1.05), breaks = c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + facet_wrap(~strength, scales = 'free')
g.orthovspara


ggsave(g.orthovspara, file = '../figures/210126 ortho vs par intscores by brlen weakloss and strongloss.pdf', scale = 0.6)
ggsave(g.orthovspara, file = '../figures/201208 ortho vs para intscores by brlen.pdf', scale = 0.6)
ggsave(g.orthovspara, file = "../figures/201121 Inteaction scores vs brlen.pdf", scale = 0.6)



#stepsquant <- select(stepsquant, -junk)
steps <- left_join(disttointprof, stepsquant, by = 'xpep')
steps <- left_join(steps, stepsquant, by = c('ypep.x' = 'ypep'))
steps <- filter(steps, steps.x == steps.y)
steps <- mutate(steps, 'sameparent' = (curparent.x == curparent.y))
steps$sameparent <- factor(steps$sameparent, levels = c("TRUE", "FALSE"))

g.stepsints <- ggplot(steps, aes(as.factor(steps.x), fill = as.factor(intscore))) + geom_bar(position = 'fill') + facet_wrap(~sameparent, labeller = labeller("FALSE" = "Different parents","TRUE"= "Same parents")) +
  scale_fill_manual(labels = c("Not Measured", "None", "Weak", "Strong"), name = "Interaction\nstrength", values = c('#A7A7A7',"#0090F3","#C300EC","#FE3131")) +
  labs(x = "Number of steps", y = "Fraction of Interactions", title = "201122 EA66k percent of interactions by steps to root") + theme(strip.background = element_rect(fill = "White"), text=  element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size=18))
g.stepsints

ggsave(g.stepsints, file = '../figures/201122 Interaction percents by steps to root.pdf', scale = 0.6)


orthonames <- read.csv('201123_filtered_orthosonly.csv')
colnames(orthonames) <- c('xpep','ypep','intscore', 'junk')
orthonames <- select(orthonames, -junk)
ortho2 <- select(orthonames, ypep, xpep, intscore)
colnames(ortho2) <- c('xpep','ypep','intscore')
orthonames <- bind_rows(orthonames, ortho2)
orthonames <- distinct(orthonames)

orthosandparas <- left_join(disttointprof, orthonames, by = c('xpep', 'ypep'))
orthosandparas <- filter(orthosandparas, logtype == 'paralogs' | !is.na(intscore.y))

ortparcor <- data.frame(x = c(.8, .8), y = c(4.5, 4.5), logtype = c('ortholog', 'paralogs'), label = c(paste('italic(r)==', round(cor(filter(orthosandparas, logtype == 'ortholog')$distance, filter(orthosandparas, logtype == 'ortholog')$pearsonlog), digits = 3)), 
                                                                                                       paste('italic(r)==', round(cor(filter(orthosandparas, logtype == 'paralogs')$distance, filter(orthosandparas, logtype == 'paralogs')$pearsonlog), digits = 3))))  
g.orthopars <- ggplot(orthosandparas, aes(distance, pearsonlog, color = logtype)) + geom_point(alpha = 0.4) + geom_smooth(method = 'lm') + facet_wrap(~logtype) +
  scale_color_manual(values = c("#0090F3","#FE3131")) + 
g.orthopars


medEA66k <- read.csv()
