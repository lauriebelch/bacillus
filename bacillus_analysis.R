
## set working directory 

### packages
library(VennDiagram)
library(wesanderson)
library(scales)
library(ggplot2)
library(plyr)
library(Cairo)
library(rstatix)
library(ggpubr)
library(WRS2)
library(FSA)
library(gridExtra)
library(MASS)
library(car)
library(multcomp)
library(jmv)
### functions
se <- function(x) sd(x)/sqrt(length(x))
belch3 <- function(x, y, z) {eval(parse(text=(paste(x, y, z, sep=""))))}
belch5 <- function(x, y, z, a, b) {eval(parse(text=(paste(x, y, z, a, b, sep=""))))}
### colour palettes
pals15 <- wes_palette(5, name="Zissou1", type='continuous') # colour palette
pals15 <- c("black", pals15[5], pals15[5], pals15[1], pals15[1], pals15[1])
pals1.6 <- wes_palette(5, name="Zissou1", type='continuous') # colour palette
pals1.7 <- wes_palette(20, name="Zissou1", type='continuous') # colour palette
pals16 <- c(pals1.6[5], pals1.6[5], "grey55", pals1.6[1], pals1.6[1])
pals16 <- c(pals1.7[1], pals1.7[6], "grey55", pals1.7[17], pals1.7[20])
pals166 <- c(pals1.6[5], pals1.6[5], "black", pals1.6[1], pals1.6[1])
pals166 <- c(pals1.7[1], pals1.7[6], "black", pals1.7[17], pals1.7[20])
pals22 <- wes_palette(5, name="Zissou1", type='continuous') # colour palette
pals222 <- c("grey55", pals22[3], pals22[1])
pals226 <- c("black", pals22[3], pals22[1])
pals227=alpha(pals222, 0.3)
pals222 <- c(pals22[3], pals22[1])
pals226 <- c(pals22[3], pals22[1])
pals227=alpha(pals222, 0.3)

### load data
sel=read.csv("Bacillus_onlyQS.csv",header=T)
sel1=read.csv("Bacillus_all.csv",header=T)

## do statistical analysis
socpam=56
j=48
param=j
if(j==14){labby= "Tajima's D"}
if(j==43){labby="Divergence (Ka)"}
if(j==44){labby="Divergence (Ks)"}
if(j==45){labby="Ka/Ks"}
if(j==19){labby= "-log(Neutrality Index)"}
if(j==21){labby= "Direction of Selection"}
if(j==22){labby= "Mcdonald-Kreitman P-value"}
if(j==48){labby=bquote("Diversity ("~pi~"/ site)")} # label for graphs of that parameter
if(j==49){labby=bquote(pi[S]~"/ site")}
if(j==50){labby=bquote(pi[N]~"/ site")}
if(j==59){labby=bquote(pi[N]~"/"~pi[S]~"/ site")}

YYLOW=range(na.omit(sel[,j]))[1]
YYHI=range(na.omit(sel[,j]))[2]

allgene=na.omit(sel1[, param])
nonsocial=na.omit(sel1[which(sel1[,34]=="Cytoplasmic"), param]) 
nonsocial=na.omit(sel1[which(sel1[,34]=="Cytoplasmic" & sel1[,53]==0), param]) 
social=na.omit(sel[which(sel[,socpam]==1), param]) 
QSn=na.omit(sel[which(sel[,socpam]==0), param]) 

allgene1 <- data.frame("var"=allgene, "social"=rep("all genes", length(allgene)))
nonsocial1 <- data.frame("var"=nonsocial, "social"=rep("non-social", length(nonsocial)))
social1 <- data.frame("var"=social, "social"=rep("Social QS", length(social)))
QSn1 <- data.frame("var"=QSn, "social"=rep("Non-Social QS", length(QSn)))

vardf11=rbind(nonsocial1, social1, QSn1)
vardf11$social <- factor(vardf11$social, levels=c("non-social", "Non-Social QS", "Social QS"))
vardf1=rbind(social1, QSn1)
vardf1$social <- factor(vardf1$social, levels=c("Non-Social QS", "Social QS"))

#################################################
vardfL=vardf11

### welch anova
anovaOneW(formula = var ~ social,
          data = vardfL,
          welchs = TRUE,
          norm = F,
          eqv = F,
          phMethod = 'gamesHowell')

#### non-parametric ##
kruskal.test(var ~ social, data = vardf11)
dunnTest(var ~ social, data = vardf11, method="bh") 

slow=quantile( vardf1$var[which(vardf1$social=="Social QS")] , 1/4)
shi=quantile( vardf1$var[which(vardf1$social=="Social QS")] , 3/4)
nlow=quantile( vardf1$var[which(vardf1$social=="Non-Social QS")] , 1/4)
nhi=quantile( vardf1$var[which(vardf1$social=="Non-Social QS")] , 3/4)

p7 <- ggplot(vardf1, aes(x=social, y=var, fill=social, color=social)) + 
  geom_hline(yintercept=median( nonsocial1$var), linetype="dashed", color = "grey75")+
  geom_jitter(width=0.3, size=2, shape=21, stroke=.1, aes(colour = social, fill=social))+
  scale_color_manual(values=pals222)+
  scale_fill_manual(values=pals222)+
  labs(x='',y=labby) +
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax=median,
               geom="crossbar", width=0.75, color=pals226)+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey75")) +
  ylim(YYLOW, YYHI)+
  theme(legend.position="none")+
  theme(aspect.ratio=1) + 
  scale_x_discrete(labels= c("Private", "Cooperative"))+
  theme(text = element_text(size=24))+#text 24 default
  theme(axis.text.x = element_text(angle =-0, vjust = 0.5, colour=pals226))+
  #ggtitle('B: Synonymous')+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  annotate(geom="text", x=0.51, y=0.9*median( nonsocial1$var), label="Non-Social",color="grey55", angle=0, size=1.5)
p7


