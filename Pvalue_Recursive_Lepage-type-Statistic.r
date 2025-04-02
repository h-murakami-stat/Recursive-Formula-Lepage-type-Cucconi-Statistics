library(data.table)
library(Rcpp)

Rcpp::sourceCpp("Lepage-type-Statistic.cpp")

### When N is Even
Pvalue.Rec.Lepage.Type.EVEN <- function(x, y) {
  start_time <- Sys.time()
  
  m1 <- length(x)
  m2 <- length(y)
  N <- m1 + m2

  data <- rank(c(x, y))
  
  W_stat <- sum(data[1:m1])
  e_w <- m1 * (m1 + m2 + 1) / 2
  v_w <- m1 * m2 * (m1 + m2 + 1) / 12
  
  MO_stat <- sum((data[1:m1] - ((N+1)/2))^2)
  e_mo <- m1 * (N^2 - 1) / 12
  v_mo <- m1 * m2 * (N + 1) * (N^2 - 4) / 180  

  A <- (W_stat - e_w) / sqrt(v_w)
  B <- (MO_stat - e_mo) / sqrt(v_mo)

  DATA.L.STAT <- A^2 + B^2
  
  w_min <- m1 * (m1 + 1) / 2
  w_max <- N * (N + 1) / 2 - m2 * (m2 + 1) / 2
  if(N %% 2 == 0){
  	mo_min <- sum_sequence1(m1)/4
  	mo_max <- (sum_sequence1(N) - sum_sequence1(m2)) / 4
  }else if(N %% 2 != 0){
  	mo_min <- sum_sequence2(m1)
  	mo_max <- sum_sequence2(N) - sum_sequence2(m2)
  }

  results <- data.table(Stat = numeric(), Prob = numeric())
  
  for (i in w_min:w_max) {
    for (j in mo_min:mo_max, 2) {
      C.LEPAGE <- Lepage_Type_r(m1, m2, i, j) ## Call Function in C++
	  L.STAT <- C.LEPAGE$TT  ## List of Statistic
	  L.PROB <- C.LEPAGE$P  ## List of Probability

	  if (L.PROB != 0) {
			results <- rbindlist(list(results, data.table(Stat = L.STAT, Prob = L.PROB)))
			}
		  }
		}

  p_value <- results[Stat >= DATA.L.STAT, sum(Prob)]

  end_time <- Sys.time()
  rec_time <- end_time - start_time
  
  return(list(
    Lep.Stat = DATA.L.STAT,
    Cucconi.Stat = DATA.L.STAT / 2,
    p_value = p_value,
    time = rec_time
  ))
}

### When N is ODD
Pvalue.Rec.Lepage.Type.ODD <- function(x, y) {
  start_time <- Sys.time()
  
  m1 <- length(x)
  m2 <- length(y)
  N <- m1 + m2

  data <- rank(c(x, y))
  
  W_stat <- sum(data[1:m1])
  e_w <- m1 * (m1 + m2 + 1) / 2
  v_w <- m1 * m2 * (m1 + m2 + 1) / 12
  
  MO_stat <- sum((data[1:m1] - ((N+1)/2))^2)
  e_mo <- m1 * (N^2 - 1) / 12
  v_mo <- m1 * m2 * (N + 1) * (N^2 - 4) / 180  

  A <- (W_stat - e_w) / sqrt(v_w)
  B <- (MO_stat - e_mo) / sqrt(v_mo)

  DATA.L.STAT <- A^2 + B^2
  
  w_min <- m1 * (m1 + 1) / 2
  w_max <- N * (N + 1) / 2 - m2 * (m2 + 1) / 2
  if(N %% 2 == 0){
  	mo_min <- sum_sequence1(m1)/4
  	mo_max <- (sum_sequence1(N) - sum_sequence1(m2)) / 4
  }else if(N %% 2 != 0){
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

  p_value <- results[Stat >= DATA.L.STAT, sum(Prob)]

  end_time <- Sys.time()
  rec_time <- end_time - start_time
  
  return(list(
    Lep.Stat = DATA.L.STAT,
    Cucconi.Stat = DATA.L.STAT / 2,
    p_value = p_value,
    time = rec_time
  ))
}

