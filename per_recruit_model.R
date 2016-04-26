.perR <- function(fe=0, POP)
{
	with(POP, {
		Fe 	= seq(0,F,delta_Fe)
		dim = c(length(age), length(Fe))

		l <- L    <- array(1, dim)
		selex     <- array(Sel, dim)
		masex     <- array(msex, dim)
		wmed      <- array(Wm, dim)
		Ma        <- l * M
		Fee       <- t(array(Fe, c(length(Fe), length(age))))
		Z         <- Ma + Fee * selex

		# 1. Age-specific total mortality and survival
			for(j in 2:A)
			{
				l[j,] <- l[j-1,] * exp(-M)
			}
			l[A,] <- l[A,] + l[A,]*exp(-M)

		# 2. Survivorship under fished conditions 
			for(j in 2:A)
			{
				L[j,] <- L[j-1,] * exp(-Z[j-1,])
			}
			L[A,] <- L[A,] + L[A,]*exp(-Z[j,])


		# 3. Incidence Functions  l:unfished -- L:fished
			fa      <- wmed * masex
			theta_E <- colSums(l*fa)
			theta_e <- colSums(L*fa)
			theta_B <- colSums(l*wmed*selex)
			theta_b <- colSums(L*wmed*selex)
			theta_q <- colSums(((L*wmed*selex)/(Ma + Fee*selex))*(1 - exp(-Ma - Fee*selex)))


		# 4. Reference points
			spr		<- theta_e / theta_E
			Fpor    = approx(spr, Fe, Fspr/100, method="linear")$y
			SPRpor  = approx(Fe, spr, Fpor, method="linear")$y

		# 5. output 

		POP$Fpor	   <- Fpor
		POP$SPRpor 	 <- SPRpor
		POP$spr  	   <- spr
		POP$Fe 		   <- Fe
		POP$theta_E  <- theta_E
		POP$theta_e  <- theta_e

		return(POP)

			})
}

.kobe <- function(POP, fig=TRUE)
  {
  with(POP, {
    # 1. Arrreglo con las biomasas desovantes, mortalidades por pesca, y puntos de referencia
    dat.cl <- as.data.frame(cbind(SB,Ftot))
    Fmsy  <- Fpor[2]
    Bmsy  <- Btar
    Blim  <- Blim 
    SBlast <- tail(SB, n=1)
    Ftlast <- tail(Ftot, n=1) 
    
    # 2. Montecarlo de las estimacione de Ft y SB para el ultimo ahno    
    # los valores 0.12 y 0.35 deben ser extraidos desde los errores estandar del mcmc
    n.mcmc <- 500
    x.mcmc <- (SBlast/Bmsy)*rlnorm(n.mcmc,0,0.18*(SBlast/Bmsy)/(SBlast/Bmsy)) 
    y.mcmc <- (Ftlast/Fmsy)*rlnorm(n.mcmc,0,0.21*(Ftlast/Fmsy)/(Ftlast/Fmsy)) 
    bw = 15
    levs = seq(0,1,0.1)
    bwx = (max(x.mcmc) - min(x.mcmc))/bw
    bwy = (max(y.mcmc) - min(y.mcmc))/bw
    est <- bkde2D(cbind(x.mcmc, y.mcmc), bandwidth = c(bwx, bwy), gridsize = c(51,51))
    est$fhat = est$fhat/max(est$fhat)
    tmp.1 <- expand.grid(x = est$x1, y = est$x2)
    tmp.2 <- melt(est$fhat)
    tmp.3 <- data.frame(x=tmp.1$x, y=tmp.1$y, z=tmp.2$value)
    
    # 3. Limites del grafico
    FmsyU <- 1.001
    FmsyL <- 0.999    
    BmsyU <- ((Bmsy/SB[1])+0.000001)/(Bmsy/SB[1])
    BmsyL <- ((Bmsy/SB[1])-0.000001)/(Bmsy/SB[1])
    rBD   = Blim/Bmsy
    yymax = 1.1*max(Ftot/Fmsy)
    xxmax = 1.1*max(SB/Bmsy)
    
    # 4. Grafica Kobe
    options(warn=-1)
    kobe <- ggplot(dat.cl) +
      theme_bw() + scale_x_continuous(limits=c(0,xxmax), breaks = sort(c(0.10, rBD, 1.00, xxmax)), 
                                      labels = c('0.10', expression(BD[lim]),expression(bold(BD/BD[MSY])), 
                                                 paste(round(xxmax,digits=0)))) + scale_y_continuous(limits=c(0,yymax),
                         breaks = round(sort(c(0.5,1.0,1.5,yymax)),digits = 3), labels = c('0.5',expression(bold(F/F[MSY])),'1.5',paste(round(yymax,digits=2))))
    
    kobe <- kobe +
      geom_rect(xmin = 0.0, xmax = rBD, ymin = 0.0, ymax = yymax, fill = 'Red', alpha = 0.65) +
      geom_rect(xmin = rBD, xmax = 1, ymin = 0.0, ymax = yymax, fill = 'Goldenrod', alpha = 0.65) +
      geom_rect(xmin = 1, xmax = xxmax, ymin = 0, ymax = 1, fill = 'SeaGreen', alpha = 0.55) +
      geom_rect(xmin = 1, xmax = xxmax, ymin = 1, ymax = yymax, fill = 'Goldenrod', alpha = 0.85) +
      geom_contour(data = tmp.3, aes(x = x, y = y, z = z, colour = ..level..), bins = 12, size = 0.4) +
      scale_colour_gradient(low = 'darkgray', high = 'white') +
      geom_path(aes(x = SB/Bmsy, y = Ftot/Fmsy), linetype = 3, size = 0.4) +
      geom_point(aes(x = SB/Bmsy, y = Ftot/Fmsy, size = yr, fill = yr), shape = 21, show_guide=FALSE) +
      scale_size(range=c(2,5)) +
      labs(y = expression(F[t]~F[MSY]~~Ratio), x = BD[t]~BD[MSY]~~Ratio, colour = 'ic:', size = 'Year:', fill = 'Year:') +
      scale_fill_gradient(limits=c(yr[1], tail(yr,n=1)), high="cornsilk3", low='cornsilk3') + # colour = yr
      geom_text(data = head(dat.cl,1), aes(x = SB/Bmsy, y = Ftot/Fmsy), label = yr[1], colour = 'black', vjust = -2, size = 3) +
      geom_text(data = tail(dat.cl,1), aes(x = SB/Bmsy, y = Ftot/Fmsy), label = yr[length(yr)],
                colour = 'black', vjust = 2, size = 3)  + 
      geom_hline(yintercept=1, linetype = "longdash", colour="Darkgrey") + 
      geom_vline(xintercept=1, linetype = "longdash", colour="Darkgrey") +
      geom_text(data = dat.cl, aes(0.25,3), label = 'Over Fished', size = 5, colour = 'white', alpha = 0.2, angle = 90) +
      geom_text(data = dat.cl, aes(3,3), label = 'Over Fishing', size = 5, colour = 'black', alpha = 0.2) +
      geom_text(data = dat.cl, aes(3,0.5), label = 'Safe-Zone', size = 5, colour = 'black', alpha = 0.2) +
      theme(legend.position="top", axis.text.y=element_text(angle = 90, hjust = .5)) + guides(fill=guide_legend(nrow=1,byrow=FALSE), size=FALSE, colour=FALSE)
    
    if(fig==TRUE) {
      png("kobeplot.png", height=2500, width=3500, res = 350, units = 'px')
      print(plot1)
      dev.off()
    } 
    return(kobe)
    options(warn=0)
  })
}

.mySeg <- function(S1=POP)
{
  with(S1, {
    segments(Fpor[1],0,Fpor[1],SPRpor[1], lty=3, col='gray50')
    segments(0,SPRpor[1],Fpor[1],SPRpor[1], lty=3, col='gray50')
    segments(Fpor[2],0,Fpor[2],SPRpor[2], lty=2, col='gray50')
    segments(0,SPRpor[2],Fpor[2],SPRpor[2], lty=2, col='gray50')
  })
}

.myPBR <- function(S1=POP)
{
  cat("------------------------------------\n")
  cat("Valores de PBR para: \n")
  cat(sprintf("SPR %d ", Fspr[1]), "y", sprintf("SPR %d \n", Fspr[2]))
  print(S1$Fpor)
  cat("------------------------------------")
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
