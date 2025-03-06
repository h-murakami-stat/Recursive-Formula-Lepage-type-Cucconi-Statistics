### Wilcoxon Rank Sum Test ####
Wilcoxon.Stat = function(x, y){
	m1 = length(x)
	n1 = length(y)
	N = m1 + n1

	x1 = sort(x)
	y1 = sort(y)
	z1 = rank(c(x1, y1))
	zx1 = z1[1:m1]

	W.STAT = sum(zx1)
	W.MU = m1*(N + 1)/2
	W.SD = sqrt(m1*n1*(N + 1)/12)

	return(list(T.W = W.STAT, T.W.EX = W.MU, T.W.SD = W.SD))
	}


### Mood Test ### 
MOOD.Stat = function(x, y){
	m1 = length(x)
	n1 = length(y)
	N = m1 + n1

	x1 = sort(x)
	y1 = sort(y)
	z1 = rank(c(x1, y1))
	zx1 = z1[1:m1]

	M.STAT = sum( (zx1 - (N + 1)/2)^2 )
	M.MU = m1 * (N^2 - 1)/12
	M.SD = sqrt(m1 * n1 * (N + 1) * (N^2 - 4)/180)

	return(list(T.MD = M.STAT, T.MD.EX = M.MU, T.MD.SD = M.SD))
	}


### Lepage-type Test ###
Lepege.Type.Stat <- function(x, y){ 	
	T1 <- (Wilcoxon.Stat(x, y)$T.W - Wilcoxon.Stat(x, y)$T.W.EX)/Wilcoxon.Stat(x, y)$T.W.SD	
	T2 <- (MOOD.Stat(x, y)$T.MD - MOOD.Stat(x, y)$T.MD.EX)/MOOD.Stat(x, y)$T.MD.SD
		
	LEP.STAT <- T1^2 + T2^2
	return(LEP.STAT)
	}

### Exact Critical Value ###
Critical_value <- function(m, n){
	start_time <- Sys.time()
		N <- m + n
		alpha <- 0.05
		z <- c(1:N)
		num.perm <- choose(N,m)
		perm <- combn(c(1:N),m)
		D.Stat = numeric(num.perm)
	
		for(i in 1:num.perm){
			xi <- z[perm[,i]]
			yi <- z[c(-perm[,i])]
			D.Stat[i] <- Lepege.Type.Stat(xi, yi)
			}

	MADDF1 <- sort(D.Stat)
	a_0.90 <- min(MADDF1[MADDF1[num.perm*0.900]<MADDF1])
	a_0.95 <- min(MADDF1[MADDF1[num.perm*0.950]<MADDF1])
	a_0.975 <- min(MADDF1[MADDF1[num.perm*0.975]<MADDF1])
	a_0.99 <- min(MADDF1[MADDF1[num.perm*0.990]<MADDF1])

	end_time <- Sys.time()

	rec_time <- end_time - start_time
	
	return(list(alpha_0.90 = a_0.90, alpha_0.95 = a_0.95, alpha_0.975 = a_0.975, alpha_0.99 = a_0.99, time = rec_time))
	}

