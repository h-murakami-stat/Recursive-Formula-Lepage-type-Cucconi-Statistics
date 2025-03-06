library(data.table)
library(Rcpp)

Rcpp::sourceCpp("Lepage-type-Statistic.cpp")

### When N is Even
Rec.Lepage.Type.EVEN <- function(m1, m2, alpha) {
start_time <- Sys.time()

	N <- m1 + m2
	w_min <- m1 * (m1 + 1) / 2
	w_max <- N * (N + 1) / 2 - m2 * (m2 + 1) / 2

	if(N %% 2 == 0){
		mo_min <- sum_sequence1(m1)/4
		mo_max <- (sum_sequence1(N) - sum_sequence1(m2)) / 4
	} else if(N %% 2 != 0){
		mo_min <- sum_sequence2(m1)
		mo_max <- sum_sequence2(N) - sum_sequence2(m2)
	}
  
	results <- data.table(Stat = numeric(), Prob = numeric())

	for(i in w_min:w_max){
		for(j in seq(mo_min, mo_max, 2)){
			C.LEPAGE <- Lepage_Type_r(m1, m2, i, j) ## Call Function in C++
			L.STAT <- C.LEPAGE$TT  ## List of Statistic
			L.PROB <- C.LEPAGE$P  ## List of Probability

			if (L.PROB != 0){
				results <- rbindlist(list(results, data.table(Stat = L.STAT, Prob = L.PROB)))
				}
			}
		}
  
	filter2 <- results[, .(Prob = sum(Prob)), by = .(Stat)][order(-Stat)]
	filter2[, CumProb := cumsum(Prob)]
	last_s <- min(filter2$Stat[filter2$CumProb <= alpha], na.rm = TRUE)
	last_p <- max(filter2$CumProb[filter2$CumProb <= alpha], na.rm = TRUE)

end_time <- Sys.time()

	rec_time <- end_time - start_time
  
	return(list(Lepage.Type.Stat = last_s, Cucconi.Stat = last_s/2, Prob = last_p, time = rec_time))
}


### When N is Even
Rec.Lepage.Type.ODD <- function(m1, m2, alpha) {
start_time <- Sys.time()

	N <- m1 + m2
	w_min <- m1 * (m1 + 1) / 2
	w_max <- N * (N + 1) / 2 - m2 * (m2 + 1) / 2
 
 	if(N %% 2 == 0){
		mo_min <- sum_sequence1(m1)/4
		mo_max <- (sum_sequence1(N) - sum_sequence1(m2)) / 4
	} else if(N %% 2 != 0){
		mo_min <- sum_sequence2(m1)
		mo_max <- sum_sequence2(N) - sum_sequence2(m2)
	}

	results <- data.table(Stat = numeric(), Prob = numeric())
  
	for (i in w_min:w_max) {
		for (j in mo_min:mo_max) {
			C.LEPAGE <- Lepage_Type_r(m1, m2, i, j) ## Call Function in C++
			L.STAT <- C.LEPAGE$TT  ## List of Statistic
			L.PROB <- C.LEPAGE$P  ## List of Probability

			if (L.PROB != 0) {
				results <- rbindlist(list(results, data.table(Stat = L.STAT, Prob = L.PROB)))
				}
			}
		}
  
		filter2 <- results[, .(Prob = sum(Prob)), by = .(Stat)][order(-Stat)]
		filter2[, CumProb := cumsum(Prob)]
		last_s <- min(filter2$Stat[filter2$CumProb <= alpha], na.rm = TRUE)
		last_p <- max(filter2$CumProb[filter2$CumProb <= alpha], na.rm = TRUE)
end_time <- Sys.time()

	rec_time <- end_time - start_time
  
	return(list(Lepage.Type.Stat = last_s, Cucconi.Stat = last_s/2, Prob = last_p, time = rec_time))
}

