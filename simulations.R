library(gss)
library(lme4)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plyr)
library(mvtnorm)
library(xtable)

############################################################################################################

# Simulate H0- and H1-consistent data and evaluate models on that
# (ideas and functions from Dave Kleinschmidt)

polynomials <- poly(1:47, 3)

# Define a data-generation-function definer
make.data.generator <- function(true.effects=c(rep(0, times=12)),
                                resid.var=1,
                                ranef.covar=diag(c(rep(1, times=12))),
                                n.subj=20
                               )
{   
  # create design matrix for our made up experiment
  data.str <- data.frame(bin   = rep(1:47, each=3),
  						 time1 = rep(polynomials[,1], each=3), 
  				     	 time2 = rep(polynomials[,2], each=3),
  				     	 time3 = rep(polynomials[,3], each=3),
  				    	 image = factor(rep(c('abig', 'bresponse', 'cother'), each=1, times=nrow(polynomials)))
  				    	)
  contrasts(data.str$image) <- contr.sum(3) #no contrast for time1/2/3, since they're not factors
  model.mat <- model.matrix(~ 1 + (time1+time2+time3)*image, data.str)
  
  generate.data <- function() {
    # sample data set under mixed effects model with random slope/intercepts
    simulated.data <- rdply(n.subj, {
      # create "true" parameters for this participant (by adding participant-specific random/error term to actual true paramenters)
      beta <- t(rmvnorm(n=1, sigma=ranef.covar)) + true.effects
      
      # create data based on those params (same value for each combination of experimental var.s at this point!)
      expected.perc <- model.mat %*% beta
            
      # create one "error value" for each trial
      epsilon <- rnorm(n=length(expected.perc), mean=0, sd=sqrt(resid.var))
      
      data.frame(data.str,
      # combine data and error term to yield "observed" datapoints
                 #percentage = expected.perc + epsilon
                 percentage = ((expected.perc + epsilon)/8) + 0.2
                )
    })
    # rdply() gives as first column the iteration number, which conveniently is subj.ID here
    #names(simulated.data)[1] <- 'subject'		
    simulated.data$subject <- factor(simulated.data[,1])
    simulated.data
  }
}


# Define a model fitter for GCA
fit.models <- function(simulated.data) {
  # fit models and extract coefs
  gca.intonly.coefs <- coefficients(summary(lmer(percentage ~ (time1 + time2 + time3) * image +
  										 	     (1 + time1 + time2 + time3 | subject),
  										 	     data = simulated.data
  										 	    )
  							    	   )
  							   )
  gca.slope.coefs <- coefficients(summary(lmer(percentage ~ (time1 + time2 + time3) * image +
  										 	   (1 + time1 + time2 + time3 | subject) +
  										 	   (1 + time1 + time2 + time3 | subject:image),
  										 	   data = simulated.data
  										 	  )
  							    	   )
  							   )
  # format output all pretty
  rbind(data.frame(model='gca.int', predictor=rownames(gca.intonly.coefs), gca.intonly.coefs),
        data.frame(model='gca.slope', predictor=rownames(gca.slope.coefs), gca.slope.coefs))
}

# Define a model fitter for SSANOVA
fit.ssanova <- function(simulated.data){
	# make ssanova model and calculate model fits
	ssa.df <- simulated.data
	ssa.fitted <- predict(object = ssanova(percentage ~ bin * image,
										   random = ~ 1|subject,
										   data = ssa.df#,
										   #seed = 2332 #this is global, I guess? upshot being, each call to the random data generator after this (i.e. each after the first call of this fitting function) returns the same data.
				   						  ),
						  newdata = subset(ssa.df, select=c('bin', 'image', 'subject')),
						  se = TRUE
						 )
	ssa.df$fit <- ssa.fitted$fit
	ssa.df$err <- ssa.fitted$se.fit
	
	# return a df of bin,image2lower,image1upper
	dummyfun <- function(bin.no){
		foo <- subset(ssa.df, bin==bin.no)
		response.upper = with(subset(foo, image=='bresponse'), mean(fit)+1.96*mean(err))
		response.lower = with(subset(foo, image=='bresponse'), mean(fit)-1.96*mean(err))
		other.upper = with(subset(foo, image=='cother'), mean(fit)+1.96*mean(err))
		other.lower = with(subset(foo, image=='cother'), mean(fit)-1.96*mean(err))
		return(c(bin=bin.no, ru=response.upper, rl=response.lower, ou=other.upper, ol=other.lower))
	}
	bounds <- do.call(rbind, lapply(1:47, dummyfun))
	bounds <- data.frame(bounds)
	return(bounds)
}

# # Loop to test calculation of ssanova fit values (result: it does work like that)
# for(i in 14:60){
	# # calculate lower bound of animate "interval" as in plot
	# anim <- subset(dfAa_gapthm, bin==i & image=='animate')
	# anim.mean <- mean(anim$ssa_fit)
	# anim.se <- mean(anim$ssa_err)
	# anim.lower <- anim.mean - 1.96*anim.se
	# # calculate upper bound of inanimate "interval" as in plot
	# inan <- subset(dfAa_gapthm, bin==i & image=='inanimate')
	# inan.mean <- mean(inan$ssa_fit)
	# inan.se <- mean(inan$ssa_err)
	# inan.upper <-inan.mean + 1.96*inan.se
	# # print handsomely
	# indicator <- ifelse(anim.lower > inan.upper, 'GREATER', 'not greater')
	# print(paste(i, anim.lower, indicator, inan.upper))
# }

# Call the definer-definer to get an actual data generation function (called noeffect.generator() here) 
noeffect.generator <- make.data.generator(true.effects=c(rep(0, times=12)), n.subj=20)
yeseffect.generator <- make.data.generator(true.effects=c(rep(0, times=9), 1, 0, 0), n.subj=20)
#image2:time2 is 10th column in model.mat

# Call the model fitters with those actual data generation functions as its data argument 1000 times.
Sys.time(); print('Starting.')
set.seed(2332)
noeffect.gca.simulations <- rdply(.n=1000, .expr=fit.models(noeffect.generator()), .progress='text')
write.csv(noeffect.gca.simulations, file='gaze_gca_no.csv', row.names=FALSE)
Sys.time(); print('Done with no-effect GCA simulation.')
set.seed(2332)
yeseffect.gca.simulations <- rdply(.n=1000, .expr=fit.models(yeseffect.generator()), .progress='text')
write.csv(yeseffect.gca.simulations, file='gaze_gca_yes.csv', row.names=FALSE)
Sys.time(); print('Done with yes-effect GCA simulation.')

set.seed(2332)
noeffect.ssa.simulations <- rdply(.n=1000, .expr=fit.ssanova(noeffect.generator()), .progress='text')
write.csv(noeffect.ssa.simulations, file='gaze_ssa_no.csv', row.names=FALSE)
Sys.time(); print('Done with no-effect SSANOVA simulation.')
set.seed(2332)
yeseffect.ssa.simulations <- rdply(.n=1000, .expr=fit.ssanova(yeseffect.generator()), .progress='text')
write.csv(yeseffect.ssa.simulations, file='gaze_ssa_yes.csv', row.names=FALSE)
Sys.time(); print('Done with yes-effect SSANOVA simulation.')
#note that there were more than 20 hours between 'Starting.' and 'Done'

# Load saved simulation results back
noeffect.gca.simulations <- read.csv('gaze_gca_no.csv')
yeseffect.gca.simulations <- read.csv('gaze_gca_yes.csv')
noeffect.ssa.simulations <- read.csv('gaze_ssa_no.csv')
yeseffect.ssa.simulations <- read.csv('gaze_ssa_yes.csv')

# Analyze simulation results
# Count false positives and negatives in GCA simulations
length(with(subset(noeffect.gca.simulations, model=="gca.slope" & predictor=="time2:image2"), which(-1.96 < t.value & t.value < +1.96)))
length(with(subset(yeseffect.gca.simulations, model=="gca.slope" & predictor=="time2:image2"), which(-1.96 < t.value & t.value < +1.96)))

# Density plot of GCA simulation t values for the effect of interest
noeffect.gca.simulations$effect <- 'no effect'; yeseffect.gca.simulations$effect <- 'yes effect'
gca.simulations <- rbind(noeffect.gca.simulations, yeseffect.gca.simulations)
gca.sim.plot <- ggplot(subset(gca.simulations, predictor=='time2:image2'), aes(x=t.value, fill=model)) +
	geom_vline(xintercept=c(-1.96,+1.96), linetype='dotted', color='#898989') +
	geom_density(alpha=0.3, color=NA) +
	ggtitle('Density of t-values for the manipulated effect') + ylab('') +
	theme_tufte(ticks=TRUE, base_size=12) +
	scale_color_tableau('colorblind10') + scale_fill_tableau('colorblind10') +
	coord_cartesian(xlim=c(-7.5, 7.5)) +
	facet_grid(. ~ effect) + theme(panel.margin=unit(5, 'mm'), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position='top')

# Count positives in SSANOVA simulations
noeffect.diffs <- vector(); for(n in unique(noeffect.ssa.simulations$.n)){
	intervals.here <- subset(noeffect.ssa.simulations, .n==n)
	diff.here <- with(intervals.here, any((rl > ou) | (ru < ol)))
	noeffect.diffs <- c(noeffect.diffs, diff.here)
}
sum(noeffect.diffs)
yeseffect.diffs <- vector(); for(n in unique(yeseffect.ssa.simulations$.n)){
	intervals.here <- subset(yeseffect.ssa.simulations, .n==n)
	diff.here <- with(intervals.here, any((rl > ou) | (ru < ol)))
	yeseffect.diffs <- c(yeseffect.diffs, diff.here)
}
sum(yeseffect.diffs)

# Check "size" of difference in SSANOVA simulations
diffs.df <- data.frame()
for(i in which(noeffect.diffs)){
	intervals.here <- subset(noeffect.ssa.simulations, .n==i)
	# find maximal difference in this model
	intervals.here$resp.above.other <- intervals.here$rl - intervals.here$ou
	intervals.here$other.above.resp <- intervals.here$ol - intervals.here$ru
	maxdiff <- with(intervals.here, max(c(max(resp.above.other), max(other.above.resp))))
	# determine "length" of difference area
	ndiff <- nrow(subset(intervals.here, (resp.above.other > 0) | (other.above.resp > 0)))
	# put the whole thing into a df
	diffs.df <- rbind(diffs.df, c(0, i, maxdiff, ndiff))
}
for(i in which(yeseffect.diffs)){
	intervals.here <- subset(yeseffect.ssa.simulations, .n==i)
	# find maximal difference in this model
	intervals.here$resp.above.other <- intervals.here$rl - intervals.here$ou
	intervals.here$other.above.resp <- intervals.here$ol - intervals.here$ru
	maxdiff <- with(intervals.here, max(c(max(resp.above.other), max(other.above.resp))))
	# determine "length" of difference area
	ndiff <- nrow(subset(intervals.here, (resp.above.other > 0) | (other.above.resp > 0)))
	# put the whole thing into a df
	diffs.df <- rbind(diffs.df, c(1, i, maxdiff, ndiff))
}; colnames(diffs.df) <- c('effect', 'model.no', 'maxdiff', 'ndiff')
summary(diffs.df)
ssa.maxdiffs.plot <- ggplot(diffs.df, aes(x=maxdiff, fill=factor(effect))) +
	geom_density(color=NA, alpha=0.3) +
	ggtitle('Maximum difference between\nconfidence intervals in SSANOVA simulation') + xlab('difference') +
	theme_tufte(ticks=TRUE, base_size=12) +
	scale_fill_tableau('colorblind10')
ssa.ndiff.plot <- ggplot(diffs.df, aes(x=ndiff, fill=factor(effect))) +
	geom_histogram(color=NA, binwidth=5, position='dodge') +
	ggtitle('Length of area of difference\nbetween confidence intervals (SSANOVA sim.)') + xlab('length, in bins\n(47 max.)') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_fill_tableau('colorblind10')

nrow(subset(diffs.df, effect==0 & ndiff==47))/nrow(subset(diffs.df, effect==0))
nrow(subset(diffs.df, effect==1 & ndiff==47))/nrow(subset(diffs.df, effect==1))

# Plot examples of SSANOVA models
ssa.simno.plot <- ggplot(subset(noeffect.ssa.simulations, .n %in% 131:139), aes(x=bin)) + #100:108 # 635:643
	geom_ribbon(aes(ymin=rl, ymax=ru, fill='response'), alpha=0.3) +
	geom_ribbon(aes(ymin=ol, ymax=ou, fill='other'), alpha=0.3) +
	ggtitle('some ssanova() model fits on\nsimulated data without effect') + xlab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_fill_tableau('colorblind10') +
	facet_wrap( ~ .n, nrow=3, ncol=3) + theme(panel.margin=unit(5, 'mm'), axis.text=element_blank(), legend.position='top')

ssa.simyes.plot <- ggplot(subset(yeseffect.ssa.simulations, .n %in% 131:139), aes(x=bin)) +
	geom_ribbon(aes(ymin=rl, ymax=ru, fill='response'), alpha=0.3) +
	geom_ribbon(aes(ymin=ol, ymax=ou, fill='other'), alpha=0.3) +
	ggtitle('some ssanova() model fits on\nsimulated data with effect') + xlab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_fill_tableau('colorblind10') +
	facet_wrap( ~ .n, nrow=3, ncol=3) + theme(panel.margin=unit(5, 'mm'), axis.text=element_blank(), legend.position='top')

# Generate a plot of the underlying functions for the simulated data
real.effects = c(rep(0, times=9), 1, 0, 0); polynomials <- poly(1:47, 3)
real.data <- data.frame(bin   = rep(1:47, each=3),
  						 time1 = rep(polynomials[,1], each=3), 
  				     	 time2 = rep(polynomials[,2], each=3),
  				     	 time3 = rep(polynomials[,3], each=3),
  				    	 image = factor(rep(c('abig', 'bresponse', 'cother'), each=1, times=nrow(polynomials)))
  				    	)
contrasts(real.data$image) <- contr.sum(3)
real.mat <- model.matrix(~ 1 + (time1+time2+time3)*image, real.data)
real.data$y <- ((real.mat %*% real.effects) / 8) + 0.2
levels(real.data$image)[levels(real.data$image)=='abig'] <- 'big'; levels(real.data$image)[levels(real.data$image)=='bresponse'] <- 'response'; levels(real.data$image)[levels(real.data$image)=='cother'] <- 'other' #rename
real.data$image <- factor(real.data$image, levels=c('other' , 'response', 'big')) #reorder
real.plot <- ggplot(real.data, aes(x=bin, y=y, color=image)) +
	geom_line() +
	#ggtitle('Underlying function for simulations with effect') +
	xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') +
	theme(axis.text=element_blank(), legend.position='top')

noeffect.data <- data.frame(bin = rep(1:47, each=3),
  					 time1 = rep(polynomials[,1], each=3), 
  				     	 time2 = rep(polynomials[,2], each=3),
  				     	 time3 = rep(polynomials[,3], each=3),
  				    	 image = factor(rep(c('big', 'response', 'other'), each=1, times=nrow(polynomials)))
  				    	)
noeffect.data$image <- factor(noeffect.data$image, levels=c('other' , 'response', 'big'))
noeffect.data$y <- 0.2
noeffect.plot <- ggplot(noeffect.data, aes(x=bin, y=y, color=image)) +
	geom_line() +
	xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') +
	theme(axis.text=element_blank(), legend.position='top')

# Generate example plots of one set of simulated data without and one with effect
set.seed(8)
example.no.data <- noeffect.generator()
levels(example.no.data$image)[levels(example.no.data$image)=='abig'] <- 'big'; levels(example.no.data$image)[levels(example.no.data$image)=='bresponse'] <- 'response'; levels(example.no.data$image)[levels(example.no.data$image)=='cother'] <- 'other' #rename
example.no.data$image <- factor(example.no.data$image, levels=c('other' , 'response', 'big')) #reorder
example.no.plot <- ggplot(example.no.data, aes(x=bin, y=percentage, color=image)) +
	geom_smooth(fill=NA) +
	#ggtitle('Example of simulated data without effect\n(averaged across simulated participants)') +
	xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') +
	theme(axis.text=element_blank(), legend.position='top')

set.seed(4)
example.yes.data <- yeseffect.generator()
levels(example.yes.data$image)[levels(example.yes.data$image)=='abig'] <- 'big'; levels(example.yes.data$image)[levels(example.yes.data$image)=='bresponse'] <- 'response'; levels(example.yes.data$image)[levels(example.yes.data$image)=='cother'] <- 'other' #rename
example.yes.data$image <- factor(example.yes.data$image, levels=c('other' , 'response', 'big')) #reorder
example.yes.plot <- ggplot(example.yes.data, aes(x=bin, y=percentage, color=image)) +
	geom_smooth(fill=NA) +
	#ggtitle('Example of simulated data with effect\n(averaged across simulated participants)') +
	xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') +
	theme(axis.text=element_blank(), legend.position='top')

# Save all plots and clear workspace

############################################################################################################

# This part simulates, analyzes, and plots GCA vs SSANOVA behavior on a model based on the adult data response-wise GCA coefficients.

polynomials <- poly(1:47, 3)

# Define a data-generation-function definer
make.data.generator <- function(true.effects=rep(0, times=12),
                                resid.var=1,
                                ranef.covar=diag(c(rep(1, times=12))),
                                n.subj=20
                               )
{   
  # create design matrix for our made up experiment
  data.str <- data.frame(bin   = rep(1:47, each=3),
  						 time1 = rep(polynomials[,1], each=3), 
  				     	 time2 = rep(polynomials[,2], each=3),
  				     	 time3 = rep(polynomials[,3], each=3),
  				    	 image = factor(rep(c('abig', 'bother', 'cresponse'), each=1, times=nrow(polynomials)))
  				    	)
  contrasts(data.str$image) <- contr.sum(3) #no contrast for time1/2/3, since they're not factors
  model.mat <- model.matrix(~ 1 + (time1+time2+time3)*image, data.str)
  
  generate.data <- function() {
    # sample data set under mixed effects model with random slope/intercepts
    simulated.data <- rdply(n.subj, {
      # create "true" parameters for this participant (by adding participant-specific random/error term to actual true paramenters)
      beta <- t(rmvnorm(n=1, sigma=ranef.covar)) + true.effects
      
      # create data based on those params (same value for each combination of experimental var.s at this point!)
      expected.perc <- model.mat %*% beta
      # ...and remove any negative values (which don't make sense if we're still saying it's percentages)
      expected.perc[expected.perc < 0] <- 0
            
      # create one "error value" for each trial
      epsilon <- rnorm(n=length(expected.perc), mean=0, sd=sqrt(resid.var))
      
      # combine data and error term to yield "observed" datapoints
      simulated.perc <- ((expected.perc + epsilon)/8)
      simulated.perc[simulated.perc < 0] <- 0
      
      data.frame(data.str,
                 #percentage = expected.perc + epsilon
                 percentage = simulated.perc
                )
    })
    # rdply() gives as first column the iteration number, which conveniently is subj.ID here
    #names(simulated.data)[1] <- 'subject'		
    simulated.data$subject <- factor(simulated.data[,1])
    simulated.data
  }
}

# Define a model fitter for GCA and SSANOVA
fit.models <- function(simulated.data) {
  # fit GCA model and extract coefs
  gca.coefs <- coefficients(summary(lmer(percentage ~ (time1 + time2 + time3) * image +
  										 	   (1 + time1 + time2 + time3 | subject) +
  										 	   (1 + time1 + time2 + time3 | subject:image),
  										 	   simulated.data
  										 	  )
  							    	   )
  							   )
  
  # fit SSANOVA model and calculate fits
  ssa.df <- simulated.data
  ssa.model <- ssanova(percentage ~ bin * image,
							random = ~ 1|subject,
							data = ssa.df
					   )
  ssa.fitted <- predict(ssa.model, newdata = subset(ssa.df, select=c('bin', 'image', 'subject')), se = TRUE)
  ssa.df$fit <- ssa.fitted$fit
  ssa.df$err <- ssa.fitted$se.fit 
  
  # make a df of bin and relevant CI boundaries of the three splines
  dummyfun <- function(bin.no){
		foo <- subset(ssa.df, bin==bin.no)
		response.lower = with(subset(foo, image=='cresponse'), mean(fit) - 1.96*mean(err))
		other.upper    = with(subset(foo, image=='bother'),    mean(fit) + 1.96*mean(err))
		big.upper      = with(subset(foo, image=='abig'),      mean(fit) + 1.96*mean(err))
		return(c(bin=bin.no, rl=response.lower, ou=other.upper, bu=big.upper))
  }
  bounds <- do.call(rbind, lapply(1:47, dummyfun))
  bounds <- data.frame(bounds)
  # calculate the two differences (resp. above other, resp. above big) for each bin from that df of CIs
  bounds$rao <- bounds$rl - bounds$ou
  bounds$rab <- bounds$rl - bounds$bu
  # count the rows where both diffs are positive (i.e. where resp. is significantly above other AND big)
  whitespace.length <- nrow(subset(bounds, rao > 0 & rab > 0))

  # extract interesting cosine diagnostics
  ssa.cos <- summary(ssa.model, diagnostics=TRUE)$cosines
  interesting.cos.names <- c('bin.cos.y', 'image.cos.y', 'interaction.cos.y', 'bin.cos.e', 'image.cos.e', 'interaction.cos.e')
  interesting.cos.values <- c(ssa.cos[1,1], ssa.cos[1,2], ssa.cos[1,3], ssa.cos[2,1], ssa.cos[2,2],ssa.cos[2,3])

  # format output all pretty
  data.frame(model = c(rep('gca', 12), rep('ssa', 7)),
  			 thing = c(rownames(gca.coefs),
  			 		   'resp.above.both.length',
  			 		   interesting.cos.names
  			 		  ),
  			 value = c(unname(gca.coefs[,3]), whitespace.length, interesting.cos.values)
  			)
}

############################################################################################################

noeffect.generator <- make.data.generator(true.effects=c(0.6, rep(0, times=11)), n.subj=20)
yeseffect.generator <- make.data.generator(true.effects=c(0.6, 0, 0, 1.6, 0, 0, 0, 3, 2.8, 0, 0, -2.6),  n.subj=20)

# Generate plot of "underlying" data function without effect
# Yes, this is three perfectly overlapping lines. Firstly, that's the point (to show they turn very much non-overlapping with randomization brought by simulated participants); secondly, use different linetypes so the lines don't overlap everywhere.
no.effects=c(0.6, rep(0, times=11))
noeff.data <- data.frame(bin   = rep(1:47, each=3),
  						 time1 = rep(polynomials[,1], each=3), 
  				     	 time2 = rep(polynomials[,2], each=3),
  				     	 time3 = rep(polynomials[,3], each=3),
  				    	 image = factor(rep(c('abig', 'bother', 'cresponse'), each=1, times=nrow(polynomials)))
  				    	)
contrasts(noeff.data$image) <- contr.sum(3)
noeff.mat <- model.matrix(~ 1 + (time1+time2+time3)*image, noeff.data)
noeff.data$y <- ((noeff.mat %*% no.effects) / 8) + 0.2
noeff.data$y[noeff.data$y < 0] <- 0
levels(noeff.data$image)[levels(noeff.data$image)=='abig'] <- 'big'; levels(noeff.data$image)[levels(noeff.data$image)=='bother'] <- 'other'; levels(noeff.data$image)[levels(noeff.data$image)=='cresponse'] <- 'response' #rename
noeff.data$image <- factor(noeff.data$image, levels=c('other', 'response', 'big')) #reorder
noeff.plot <- ggplot(noeff.data, aes(x=bin, y=y, color=image, linetype=image)) +
	geom_line() +
	ggtitle('Underlying function for simulations without effect') + xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') + scale_linetype_manual(values=c('11', '88', 'FF')) +
	theme(axis.text=element_blank(), legend.position='top')


# Generate plot of "underlying" data simulation function
real.effects = c(0.6, 0, 0, 1.6, 0, 0, 0, 3, 2.8, 0, 0, -2.6); polynomials <- poly(1:47, 3)
real.data <- data.frame(bin   = rep(1:47, each=3),
  						 time1 = rep(polynomials[,1], each=3), 
  				     	 time2 = rep(polynomials[,2], each=3),
  				     	 time3 = rep(polynomials[,3], each=3),
  				    	 image = factor(rep(c('abig', 'bother', 'cresponse'), each=1, times=nrow(polynomials)))
  				    	)
contrasts(real.data$image) <- contr.sum(3)
real.mat <- model.matrix(~ 1 + (time1+time2+time3)*image, real.data)
real.data$y <- ((real.mat %*% real.effects) / 8) + 0.2
real.data$y[real.data$y < 0] <- 0
levels(real.data$image)[levels(real.data$image)=='abig'] <- 'big'; levels(real.data$image)[levels(real.data$image)=='bother'] <- 'other'; levels(real.data$image)[levels(real.data$image)=='cresponse'] <- 'response' #rename
real.data$image <- factor(real.data$image, levels=c('other', 'response', 'big')) #reorder
real.plot <- ggplot(real.data, aes(x=bin, y=y, color=image)) +
	geom_line() +
	ggtitle('Underlying function for simulations with effect') + xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') +
	theme(axis.text=element_blank(), legend.position='top')

# Generate example plot of one set of simulated data with effect
set.seed(47) #3 or 47 are good
example.yes.data <- yeseffect.generator()
levels(example.yes.data$image)[levels(example.yes.data$image)=='abig'] <- 'big'; levels(example.yes.data$image)[levels(example.yes.data$image)=='bother'] <- 'other'; levels(example.yes.data$image)[levels(example.yes.data$image)=='cresponse'] <- 'response' #rename
example.yes.data$image <- factor(example.yes.data$image, levels=c('other', 'response', 'big')) #reorder
example.yes.plot <- ggplot(example.yes.data, aes(x=bin, y=percentage, color=image)) +
	geom_smooth(fill=NA) +
	ggtitle('Example of simulated data with effect\n(averaged across simulated participants)') + xlab('') + ylab('') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_color_tableau('colorblind10') +
	theme(axis.text=element_blank(), legend.position='top')

############################################################################################################

# Run simulations
Sys.time(); print('Starting.')
set.seed(2332)
noeffect.simulations <- rdply(.n=1000, .expr=fit.models(noeffect.generator()), .progress='text')
write.csv(noeffect.simulations, file='gaze_complex_noeffect.csv', row.names=FALSE)
Sys.time(); print('Done with no-effect simulations.')
set.seed(2332)
yeseffect.simulations <- rdply(.n=1000, .expr=fit.models(yeseffect.generator()), .progress='text')
write.csv(yeseffect.simulations, file='gaze_complex_yeseffect.csv', row.names=FALSE)
Sys.time(); print('Done with yes-effect simulations.')
# Note that this took about 18h to run.

# Load saved simulation results back
noeffect.simulations <- read.csv('gaze_complex_noeffect.csv')
yeseffect.simulations <- read.csv('gaze_complex_yeseffect.csv')

############################################################################################################

# Analyze simulation results

# Count false positives and false negatives in GCA
with(subset(noeffect.simulations, model=='gca' & !thing=='(Intercept)' & abs(value) > 1.96), length(unique(.n)))
with(subset(noeffect.simulations, model=='gca' & !thing=='(Intercept)' & abs(value) > 1.96), table(thing))
with(subset(yeseffect.simulations, model=='gca' & thing %in% c('time3', 'time1:image2', 'time2:image1', 'time3:image2') & abs(value) > 1.96), length(unique(.n)))
with(subset(yeseffect.simulations, model=='gca' & !thing %in% c('(Intercept)', 'time3', 'time1:image2', 'time2:image1', 'time3:image2') & abs(value) > 1.96), length(unique(.n)))
with(subset(yeseffect.simulations, model=='gca' & !thing %in% c('(Intercept)', 'time3', 'time1:image2', 'time2:image1', 'time3:image2') & abs(value) > 1.96), table(thing))

# Density plot of GCA simulation t-values for the effects of interest
noeffect.simulations$effect <- 'no effect'; yeseffect.simulations$effect <- 'yes effect'
both.gca.simulations <- rbind(subset(noeffect.simulations, model=='gca'), subset(yeseffect.simulations, model=='gca'))
sim.plot <- ggplot(subset(both.gca.simulations, thing %in% c('time3', 'time1:image2', 'time2:image1', 'time3:image2')), aes(x=value, fill=effect)) +
	geom_vline(xintercept=c(-1.96,+1.96), linetype='dotted', color='#898989') +
	geom_density(alpha=0.3, color=NA) +
	ggtitle('Density of t-values for the manipulated effects\n(simulation no.3)') + ylab('') + xlab('t-value') +
	theme_tufte(ticks=TRUE, base_size=12) +
	scale_fill_tableau('colorblind10') +
	facet_grid(thing ~ effect) + theme(panel.margin.x=unit(1, 'cm'), legend.position='top', axis.text.y=element_blank(), axis.ticks.y=element_blank())

# Count false positives and false negatives in SSANOVA
with(subset(noeffect.simulations, thing=='resp.above.both.length'), hist(value))
with(subset(noeffect.simulations, thing=='resp.above.both.length'), length(which(value > 0)))
with(subset(noeffect.simulations, thing=='resp.above.both.length'), length(which(value == 47)))
with(subset(noeffect.simulations, thing=='resp.above.both.length'), length(which(value > 0 & value < 47)))
with(subset(yeseffect.simulations, thing=='resp.above.both.length'), hist(value))
with(subset(yeseffect.simulations, thing=='resp.above.both.length'), length(which(value == 0)))
with(subset(yeseffect.simulations, thing=='resp.above.both.length'), max(value))

# Histogram of length of difference area
both.ssa.simulations <- rbind(subset(noeffect.simulations, model=='ssa'), subset(yeseffect.simulations, model=='ssa'))
ssa.ndiff.plot <- ggplot(subset(both.ssa.simulations, thing=='resp.above.both.length'), aes(x=value, fill=effect)) +
	geom_histogram(color=NA, binwidth=5, position='dodge') +
	ggtitle('Length of area of difference\nbetween confidence intervals\n(SSANOVA, simulation no. 3)') + xlab('length, in bins\n(47 max.)') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_fill_tableau('colorblind10')
ssa.ndiff.plot.noextremes <- ggplot(subset(both.ssa.simulations, thing=='resp.above.both.length' & !value %in% c(0,47)), aes(x=value, fill=effect)) +
	geom_histogram(color=NA, binwidth=5, position='dodge') +
	ggtitle('Length of area of difference between\nconfidence intervals, excluding 0 and 47\n(SSANOVA, simulation no. 3)') + xlab('length, in bins\n(47 max.)') +
	theme_tufte(ticks=FALSE, base_size=12) +
	scale_fill_tableau('colorblind10')

# Plots of the cosine diagnostics (divided up between cos.e and cos.y because their values are far apart and the automatic scaling in ggplot makes it hard to see anything in this case)
cos.es <- subset(both.ssa.simulations, thing %in% c('bin.cos.e', 'image.cos.e', 'interaction.cos.e'))
cos.e.plot <- ggplot(cos.es, aes(x=value, fill=effect)) +
	geom_density(alpha=0.3, color=NA) +
	ggtitle('Densities of cosines between\neffect and error terms\n(SSANOVA, simulation no.3)') + ylab('') + xlab('cos') +
	theme_tufte(ticks=TRUE, base_size=12) +
	scale_fill_tableau('colorblind10') +
	facet_grid(thing ~ .) + theme(legend.position='top', axis.text.y=element_blank(), axis.ticks.y=element_blank())

cos.ys <- subset(both.ssa.simulations, thing %in% c('bin.cos.y', 'image.cos.y', 'interaction.cos.y'))
cos.y.plot <- ggplot(cos.ys, aes(x=value, fill=effect)) +
	geom_density(alpha=0.3, color=NA) +
	ggtitle('Densities of cosines between\neffect and response terms\n(SSANOVA, simulation no.3)') + ylab('') + xlab('cos') +
	theme_tufte(ticks=TRUE, base_size=12) +
	scale_fill_tableau('colorblind10') +
	facet_grid(thing ~ .) + theme(legend.position='top', axis.text.y=element_blank(), axis.ticks.y=element_blank())


# Check whether the models that had a significant t-value (in GCA) also had a difference (in SSANOVA)
noeffect.gca.diffs <- unique(noeffect.simulations$.n[with(noeffect.simulations, which(model=='gca' & !thing=='(Intercept)' & abs(value) > 1.96))])
noeffect.ssa.diffs <- vector(); for(n in unique(noeffect.simulations$.n)){
	if((subset(noeffect.simulations, model=='ssa' & .n==n & thing=='resp.above.both.length')$value > 0) & (subset(noeffect.simulations, model=='ssa' & .n==n & thing=='resp.above.both.length')$value < 47)){
		noeffect.ssa.diffs <- c(noeffect.ssa.diffs, n)
	}
}
length(which(noeffect.gca.diffs %in% noeffect.ssa.diffs))/length(noeffect.ssa.diffs)

yeseffect.gca.diffs <- unique(yeseffect.simulations$.n[with(yeseffect.simulations, which(model=='gca' & !thing=='(Intercept)' & abs(value) > 1.96))])
yeseffect.ssa.diffs <- vector(); for(n in unique(yeseffect.simulations$.n)){
	if((subset(yeseffect.simulations, model=='ssa' & .n==n & thing=='resp.above.both.length')$value > 0) & (subset(yeseffect.simulations, model=='ssa' & .n==n & thing=='resp.above.both.length')$value < 47)){
		yeseffect.ssa.diffs <- c(yeseffect.ssa.diffs, n)
	}
}
length(which(yeseffect.gca.diffs %in% yeseffect.ssa.diffs))/length(noeffect.ssa.diffs)

# Save all plots
