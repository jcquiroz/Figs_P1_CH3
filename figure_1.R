rm(list=ls())
.libPaths() 
library(ggplot2)
library(ggthemes)

# ------------- plot 1a ------------------

gt <- read.csv2('TAC_fig1.csv', sep=',')

fig1a <- ggplot(gt, aes(x=year, y=tac, fill=source)) + geom_bar(stat="identity",position="dodge")#, alpha = 0.80)
fig1a <- fig1a + labs(y = 'Removal (t)', x = 'Year', fill = 'Source:') 
fig1a <- fig1a + theme_few() + scale_fill_wsj("colors6") #scale_fill_tableau("colorblind10")
fig1a <- fig1a + geom_vline(xintercept=2013.5, linetype = "longdash", colour="Darkgrey") 
#fig1a <- fig1a + geom_hline(yintercept=c(1000,2000,300), linetype = "longdash", colour="Darkgrey", size=0.5) 
fig1a <- fig1a + geom_text(aes(2014.8,3100), label = 'New management\n Framework', size = 3, colour = 'Darkgrey', alpha = 0.5)
fig1a <- fig1a + theme(legend.position="top", axis.text.y=element_text(angle = 90, hjust = .5)) + scale_x_continuous(breaks = c(2005,2007,2009,2011,2013,2015))
fig1a

png('fig1a.png', height=2000, width=2400, res = 400, units = 'px')
par(mfrow = c(1,1), mar=c(4,4,3,1))
print(fig1a)
dev.off()

# ------------- plot 1b ------------------
require(reshape2)
require(reshape2)
require(KernSmooth)

source('per_recruit_model.R')
source('read.admb.r')
fit      <- read.admb('2015')    #  stock assessment output 
yr       <- fit$Yr
ssb      <- fit$SSB
f1       <- apply(fit$TotF,1,max)

S1 <- list(SB=ssb[33:69,2], Ftot=f1, Fpor=c(0.0777, 0.0889), Btar=29567, Blim=29567/2)
fig1b <- .kobe(S1, fig=FALSE)  # use TRUE to plot
fig1b

png('fig1b.png', height=2000, width=2400, res = 400, units = 'px')
par(mfrow = c(1,1), mar=c(4,4,3,1))
print(fig1b)
dev.off()


# ------------- plot 1c ------------------

hcr <- read.table('hcrs.mse')
names(hcr) <- c('Type','B','F')

fig1c <- ggplot(hcr, aes(x=B, y=F, colour=Type)) + geom_line(stat="identity",position="dodge", size=1.1, alpha = 0.80)
fig1c <- fig1c + labs(y = 'Target fishing mortality ratio', x = 'Target spawning biomass ratio', colour = 'Type HCR:') 
fig1c <- fig1c + theme_few() + scale_colour_wsj("colors6") #scale_colour_tableau("colorblind10") #
fig1c <- fig1c + geom_vline(xintercept=c(0.4), linetype = "longdash", colour="Darkgrey") 
fig1c <- fig1c + geom_hline(yintercept=c(1), linetype = "longdash", colour="Darkgrey", size=0.5) 
#fig1c <- fig1c + geom_text(aes(2014.8,3100), label = 'New management\n Framework', size = 3, colour = 'Darkgrey', alpha = 0.5)
fig1c <- fig1c + theme(legend.position="top", axis.text.y=element_text(angle = 90, hjust = .5)) 
fig1c <- fig1c + scale_x_continuous(breaks = c(0.2,0.4), 
                                  labels = c(expression(bold(B[lim])),expression(bold(B[t]/B[MSY]))))
fig1c <- fig1c + scale_y_continuous(breaks = c(0.0,1,1.1), labels = c('',expression(bold(F[t]/F[MSY])),''), limits=c(0,1.1))
fig1c  

png('fig1c.png', height=2000, width=2400, res = 400, units = 'px')
par(mfrow = c(1,1), mar=c(4,4,3,1))
print(fig1c)
dev.off()

# ------------- plot 1all ------------------

png('../images/fig1all.png', height=6000, width=2400, res = 400, units = 'px')
par(mfrow = c(1,1), mar=c(4,4,3,1))
multiplot(fig1a, fig1b, fig1c, cols=1)
dev.off()


# ------------- plot 2 ------------------
library(dplyr)

nb   <- 60
ns   <- 30000
br   <- c(0,.55,.75,.9,1)

hr   <- filter(data.frame(h=rbeta(ns, 5.94, 1.97)),  h>=0.45 & h<=1)
hs   <- hist(hr$h,nb)
hs   <- data.frame(hs$density,hs$counts,hs$mids,up=hs$breaks[-1])
hs   <- mutate(hs, hseg = cut(up, breaks = br, labels = c('Base-case','Low Prod.','Medium Prod.','High Prod.')))
head(hs)
pro  <- round(diff(pbeta(br, 5.94, 1.97, ncp = 0, lower.tail = TRUE, log.p = FALSE)),3)

h     = seq(0.45,1,0.0001)
Ph    <- dbeta(h, 5.94, 1.97, ncp = 0) 
hc    <- data.frame(h,Ph)

br[1] = 0.45
fig2 <- ggplot(hs, aes(x=hs.mids, y=hs.density, fill=hseg)) + geom_bar(stat="identity",position="dodge")#, alpha = 0.80)
fig2 <- fig2 + labs(y = 'Density', x = 'Steepness', fill = 'Productivity\nlevels:')  + ylim(0,3.5)
fig2 <- fig2 + theme_few() + scale_fill_tableau("colorblind10") #scale_fill_wsj("colors6") #
fig2 <- fig2 + geom_vline(xintercept=br, linetype = "longdash", colour="Darkgrey") 
fig2 <- fig2 + theme(legend.position="top") + scale_x_continuous(breaks = br) 
fig2 <- fig2 + geom_text(aes(0.5,3.35), label = paste('P=',pro[1]), size = 3.5, colour = 'Darkgrey') 
fig2 <- fig2 + geom_text(aes(0.65,3.35), label = paste('P=',pro[2]), size = 3.5, colour = 'Darkgrey') 
fig2 <- fig2 + geom_text(aes(0.825,3.35), label = paste('P=',pro[3]), size = 3.5, colour = 'Darkgrey') 
fig2 <- fig2 + geom_text(aes(0.95,3.35), label = paste('P=',pro[4]), size = 3.5, colour = 'Darkgrey') 
fig2

png('../images/fig2.png', height=2000, width=2400, res = 400, units = 'px')
par(mfrow = c(1,1), mar=c(4,4,3,1))
print(fig2)
dev.off()



















