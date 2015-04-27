# this expects an input df with named columns for:
# PID, TID, gaze_time, gaze_x, gaze_y.

setwd('/Users/dmb131/Documents/')
START_TIME <- Sys.time()


# 0. Set parameters.
#"a fixation filter that compiled groups of gaze points occurring within 50 pixels and within at most 200 ms of each other."
RADIUS = 25 #pixels
T_THRESHOLD = 0.2 #seconds


# 1. Load existing data csv files.
gaze_all <- readRDS('g_all.rds')
points_vector = numeric() #this will eventually be row indices. at 60Hz with a 200ms threshold, it should never be longer than 13 items.
gaze_all$mulak <- NA
fixation_no = 1

data <- gaze_all[complete.cases(gaze_all$gaze_x),]

# 2. Do the thing.
for(i in 1:nrow(data)){

	#vector of previous points empty?
	if(length(points_vector > 0)){
		
		#current from from same PID and TID as previous ones?
		current_pid = data$PID[i]
		current_tid = data$TID[i]
		last_pid = data$PID[tail(points_vector,1)]
		last_tid = data$TID[tail(points_vector,1)]
		if( (current_pid==last_pid) & (current_tid==last_tid) ){
			same_pidtid = TRUE
		} else {
			same_pidtid = FALSE
		}
		
		#current row under time threshold (counting from first row)?
		current_time = data$gaze_time[i]
		first_time = data$gaze_time[head(points_vector,1)]
		if(current_time <= (first_time + T_THRESHOLD)){
			under_threshold = TRUE
		} else {
			under_threshold = FALSE
		}
		
		#current and previous rows all within circle?
		x_center = mean(data$gaze_x[c(points_vector,i)])
		y_center = mean(data$gaze_y[c(points_vector,i)])
		leq_vector = logical()
		for(j in c(points_vector,i)){
			if( (data$gaze_x[j] - x_center)^2 + (data$gaze_y[j] - y_center)^2 <= RADIUS^2 ){
				leq_vector <- c(leq_vector, TRUE)
			} else {
				leq_vector <- c(leq_vector, FALSE)
			}
		}
		if(all(leq_vector)){
			in_circle = TRUE
		} else {
			in_circle = FALSE
		}
		
		if(!(same_pidtid & under_threshold & in_circle)){
			if(length(points_vector) > 1){

				#mark previous points as fixation
				data$mulak[points_vector] <- paste('fixation', fixation_no, sep='')
				fixation_no = fixation_no + 1
			}
			
			#empty vector of previous points
			points_vector = numeric()
		}
	}
	
	#add current point to vector
	points_vector <- c(points_vector, i)
	
	if(i%%5000==0){
		print(paste(Sys.time(), 'Finished with row #', i))
	}
}


# 3. Save back to csv.
#write.csv(data, file='/Users/dmb131/Documents/A_gaze_mulak.csv', row.names=FALSE)

gaze_all[rownames(data),] <- data
saveRDS(object=gaze_all, file='g_all_M.rds', compress=TRUE)
FINISH_TIME <- Sys.time()
FINISH_TIME - START_TIME #took <3min
print(fixation_no) #20030
