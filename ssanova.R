library(gss)
data(eyetrack)
library(lme4)
library(ggplot2)
library(ggthemes)
library(gridExtra)

# Re-organize data frame so that it has percentages
et <- within(eyetrack, image <- factor(ifelse(color==1 & object==1, 'target', ifelse(color==1, 'color distractor', ifelse(object==1, 'object distractor', 'empty')))))
et$perc <- et$cnt/6

# The original eyetrack dataset has counts, with 0 (for image at time point) being unspecified/implicit.
# Therefore, add 0 percentages for all those
df.0 <- expand.grid(time=unique(et$time), id=unique(et$id), image=levels(et$image), perc=0, KEEP.OUT.ATTRS=FALSE)
et <- et[with(et, order(id,time,image)),]; df.0 <- df.0[with(df.0, order(id,time,image)),]
foo <- merge(et, df.0, by=c('time', 'id', 'image'), all=TRUE)
foo$perc.x[is.na(foo$perc.x)] <- 0
et <- subset(foo, select=-perc.y); names(et)[7] <- 'perc'

# Add orthogonalized polynomial time terms
et <- within(et, time.incr <- match(time, unique(time)))
polynomials3 <- poly(1:max(et$time.incr), 3)
et[,paste('time3', 1:3, sep='')] <- polynomials3[et$time.incr, 1:3]
polynomials7 <- poly(1:max(et$time.incr), 7)
et[,paste('time7', 1:7, sep='')] <- polynomials7[et$time.incr, 1:7]


basic.plot <- ggplot(et, aes(x=time, y=perc, color=image)) +
	stat_summary(fun.data=mean_se, geom='line', alpha=0.5, linetype='solid') +
	theme_tufte(ticks=FALSE, base_size=12) + scale_color_tableau('colorblind10') + 
	ylab('proportion') + coord_cartesian(ylim=c(0,1)) +
	theme(legend.position='top', plot.margin=unit(rep(1,4), 'mm')) + scale_x_continuous(breaks=seq(from=-500,to=1000,by=500)) +
	guides(color=guide_legend(nrow=2, byrow=TRUE))

# Example spline (without SSANOVA): needs to have exactly one y for every x
et.color <- data.frame()
for(t in unique(et$time)){
	et.color <- rbind(et.color, c(t, mean(subset(et, image=='color distractor' & time==t)$perc)))
}; colnames(et.color) <- c('time', 'perc')
color.plot <- ggplot(et.color, aes(x=time, y=perc)) +
	geom_line(linetype='solid') +
	theme_tufte(ticks=FALSE, base_size=12) + 
	ylab('proportion') + coord_cartesian(ylim=c(0,0.35)) +
	theme(plot.margin=unit(rep(1,4), 'mm')) + scale_x_continuous(breaks=seq(from=-500,to=1000,by=500))
et.color$spline5 <- smooth.spline(et.color$time, et.color$perc, nknots=5)$y
et.color$spline10 <- smooth.spline(et.color$time, et.color$perc, nknots=10)$y
et.color$spline20 <- smooth.spline(et.color$time, et.color$perc, nknots=20)$y
color.splines.plot <- ggplot(et.color, aes(x=time, y=perc)) +
	geom_line(linetype='dotted') +
	geom_line(aes(x=time, y=spline5), color='#006BA4') +
	geom_line(aes(x=time, y=spline10), color='#FF800E') +
	geom_line(aes(x=time, y=spline20), color='#ABABAB') +
	theme_tufte(ticks=FALSE, base_size=12) + 
	ylab('proportion') + coord_cartesian(ylim=c(0,0.35)) +
	theme(plot.margin=unit(rep(1,4), 'mm')) + scale_x_continuous(breaks=seq(from=-500,to=1000,by=500))

# Example SSANOVA
et.ssa <- ssanova(perc ~ time * image, random=~1|id, data=et, seed=1)
et$ssa <- predict(et.ssa, newdata=subset(et, select=c('time', 'image', 'id')), se=TRUE)$fit
et$ssa.lower <- et$ssa - (1.96 * predict(et.ssa, newdata=subset(et, select=c('time', 'image', 'id')), se=TRUE)$se.fit)
et$ssa.upper <- et$ssa + (1.96 * predict(et.ssa, newdata=subset(et, select=c('time', 'image', 'id')), se=TRUE)$se.fit)
basic.ssa.plot <- basic.plot +
	geom_line(data=et, aes(y=ssa), stat='summary', fun.y=mean, linetype='dashed', show_guide=FALSE) +
	geom_ribbon(data=et, aes(ymin=ssa.lower, ymax=ssa.upper, fill=image), color=NA, alpha=0.3, show_guide=FALSE) +
	scale_fill_tableau('colorblind10')

#######

# Take a subset of the data, id 61 to id 66
sixty <- subset(et, id %in% seq(from=61, to=66, by=1), select=-c(color, object, cnt))

# Set up a basic plot of that subset
sixty.plot <- ggplot(sixty, aes(x=time, y=perc, color=image)) +
	stat_summary(fun.data=mean_se, geom='line', alpha=0.5, linetype='solid') +
	theme_tufte(ticks=FALSE, base_size=12) + scale_color_tableau('colorblind10') + 
	ylab('proportion') + coord_cartesian(ylim=c(0,1)) +
	theme(legend.position='top', plot.margin=unit(rep(1,4), 'mm')) + scale_x_continuous(breaks=seq(from=-500,to=1000,by=500)) +
	guides(color=guide_legend(nrow=2, byrow=TRUE))


# Fit and plot grand mean as "model"
mean.plot <- sixty.plot + 
	geom_segment(aes(x=min(sixty$time), y=mean(subset(sixty, image=='color distractor')$perc), xend=max(sixty$time), yend=mean(subset(sixty, image=='color distractor')$perc)), color='#006BA4', linetype='dashed') +
	geom_segment(aes(x=min(sixty$time), y=mean(subset(sixty, image=='empty')$perc), xend=max(sixty$time), yend=mean(subset(sixty, image=='empty')$perc)), color='#FF800E', linetype='dashed') +
	geom_segment(aes(x=min(sixty$time), y=mean(subset(sixty, image=='object distractor')$perc), xend=max(sixty$time), yend=mean(subset(sixty, image=='object distractor')$perc)), color='#ABABAB', linetype='dashed') +
	geom_segment(aes(x=min(sixty$time), y=mean(subset(sixty, image=='target')$perc), xend=max(sixty$time), yend=mean(subset(sixty, image=='target')$perc)), color='#595959', linetype='dashed')

# Fit and plot basic linear model
sixty$lm <- predict(lm(perc ~ time * image, data=sixty))
lm.plot <- sixty.plot +
	geom_line(data=sixty, aes(y=lm), linetype='dashed', show_guide=FALSE)

# Fit SSANOVA
sixty.ssa <- ssanova(perc ~ time * image, data=sixty, seed=1, random=~1|id)
sixty$ssa <- predict(sixty.ssa, newdata=subset(sixty, select=c('time', 'image')), se=TRUE)$fit
sixty$ssa.lower <- sixty$ssa - (1.96 * predict(sixty.ssa, newdata=subset(sixty, select=c('time', 'image')), se=TRUE)$se.fit)
sixty$ssa.upper <- sixty$ssa + (1.96 * predict(sixty.ssa, newdata=subset(sixty, select=c('time', 'image')), se=TRUE)$se.fit)
#additivity?
#re-did without "emtpy" level -- slightly different shapes, but that could well be randomness

# Plot SSANOVA
ssa.plot <- sixty.plot +
	geom_line(data=sixty, aes(y=ssa), stat='summary', fun.y=mean, linetype='dashed', show_guide=FALSE) +
	geom_ribbon(data=sixty, aes(ymin=ssa.lower, ymax=ssa.upper, fill=image), color=NA, alpha=0.3, show_guide=FALSE) +
	scale_fill_tableau('colorblind10')
	

# Fit GCAs
sixty.gca3 <- lmer(perc ~ (time31 + time32 + time33) * image +
				  (1 + time31 + time32 + time33 | id) +
				  (1 + time31 + time32 + time33 | id:image),
				  data=sixty, control=lmerControl(optimizer='bobyqa'), REML=FALSE
				 )
sixty$gca3 <- fitted(sixty.gca3)
sixty.gca7 <- lmer(perc ~ (time71 + time72 + time73 + time74 + time75 + time76 + time77) * image +
				  (1 + time71 + time72 + time73 + time74 + time75 + time76 + time77 | id) +
				  (1 + time71 + time72 + time73 + time74 + time75 + time76 + time77 | id:image),
				  data=sixty, control=lmerControl(optimizer='bobyqa'), REML=FALSE
				 )
#throws convergence warnings!
sixty$gca7 <- fitted(sixty.gca7)

# Plot GCAs
gca3.plot <- sixty.plot +
	geom_line(data=sixty, aes(y=gca3), stat='summary', fun.y=mean, linetype='dashed', show_guide=FALSE)
gca7.plot <- sixty.plot +
	geom_line(data=sixty, aes(y=gca7), stat='summary', fun.y=mean, linetype='dashed', show_guide=FALSE)

# Save all plots
