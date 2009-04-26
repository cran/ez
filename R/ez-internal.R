ezANOVA_levene <-
function (y, ...) {
	form <- y
	mf <- model.frame(form, ...)
	if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]])))) stop("Levene's test is not appropriate with quantitative explanatory variables.")
	y <- mf[,1]
	if(dim(mf)[2]==2) group <- mf[,2]
	else {
		if (length(grep("\\+ | \\| | \\^ | \\:",form))>0) stop("Model must be completely crossed formula only.")
		group <- interaction(mf[,2:dim(mf)[2]])
	}
	if (!is.numeric(y)) 
		stop(deparse(substitute(y)), " is not a numeric variable")
	if (!is.factor(group)) {
		warning(deparse(substitute(group)), " coerced to factor.")
		group <- as.factor(group)
	}
	meds <- tapply(y, group, median, na.rm = TRUE)
	resp <- abs(y - meds[group])
	table <- as.data.frame(anova(lm(resp ~ group))[, c(1,2, 4, 5)])
	to_return = data.frame(table$D[1],table$D[2],table$S[1],table$S[2],table$F[1],table$P[1])
	names(to_return)=c("DFn", "DFd", "SSn", "SSd", "F", "p")
	to_return$"p<.05"=ifelse(to_return$p<.05,'*','')
	return(to_return)
}

ezANOVA_summary <-
function(object){
	to_return=list()
	GG <- function(SSPE, P){
		p <- nrow(SSPE)
		if (p < 2) return(NA) 
		lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
		lambda <- lambda[lambda > 0]
		((sum(lambda)/p)^2)/(sum(lambda^2)/p)
	}
	HF <- function(gg, error.df, p){
		((error.df + 1)*p*gg - 2)/(p*(error.df - p*gg))
	}
	mauchly <- function (SSD, P, df) {
		# most of this function borrowed from stats:::mauchly.test.SSD
		if (nrow(SSD) < 2) return(c(NA, NA))
		Tr <- function (X) sum(diag(X))
		p <- nrow(P)
		I <- diag(p)
		Psi <- t(P) %*% I %*% P 
		B <- SSD 
		pp <- nrow(SSD) 
		U <- solve(Psi, B)
		n <- df 
		logW <- log(det(U)) - pp * log(Tr(U/pp))
		rho <- 1 - (2 * pp^2 + pp + 2)/(6 * pp * n)
		w2 <- (pp + 2) * (pp - 1) * (pp - 2) * (2 * pp^3 + 6 * pp^2 + 
				3 * p + 2)/(288 * (n * pp * rho)^2)
		z <- -n * rho * logW
		f <- pp * (pp + 1)/2 - 1
		Pr1 <- pchisq(z, f, lower.tail = FALSE)
		Pr2 <- pchisq(z, f + 4, lower.tail = FALSE)
		pval <- Pr1 + w2 * (Pr2 - Pr1)
		c(statistic = c(W = exp(logW)), p.value = pval)
	}		
	test.statistic <- 1:4
	nterms <- length(object$terms)
	error.df <- object$error.df
	table <- data.frame(matrix(0, nterms, 8))
	table2 <- data.frame(matrix(0, nterms, 7))
	table3 <- data.frame(matrix(0, nterms, 4))
	table3[,1] <- table2[,1] <- table[,1] <- object$terms
	colnames(table) <- c("Effect","DFn", "DFd", "SSn", "SSd", "F", "p", "p<.05")
	colnames(table2) <- c("Effect","GGe", "p[GG]", "p[GG]<.05", "HFe", "p[HF]","p[HF]<.05")
	colnames(table3) <- c("Effect","W", "p", "p<.05")
	for (term in 1:nterms){
		SSP <- object$SSP[[term]]
		SSPE <- object$SSPE[[term]]
		P <- object$P[[term]]
		p <- ncol(P)
		PtPinv <- solve(t(P) %*% P)
		gg <- GG(SSPE, P)
		table[term, "SSn"] <- sum(diag(SSP %*% PtPinv))
		table[term, "SSd"] <- sum(diag(SSPE %*% PtPinv))
		table[term, "DFn"] <- object$df[term] * p
		table[term, "DFd"] <- error.df * p
		table[term, "F"] <-  (table[term, "SSn"]/table[term, "DFn"])/
			(table[term, "SSd"]/table[term, "DFd"])
		table[term, "p"] <- pf(table[term, "F"], table[term, "DFn"],
			table[term, "DFd"], lower.tail=FALSE)
		table[term, "p<.05"] = ifelse(table[term, "p"]<.05,'*','')
		table2[term, "GGe"] <- gg
		table2[term, "HFe"] <- HF(gg, error.df, p)
		table3[term,2:3] <- mauchly(SSPE, P, object$error.df)
		table3[term, "p<.05"] = ifelse(table3[term, "p"]<.05,'*','')		
	}
	to_return$ANOVA=as.data.frame(table)
	table3=table3[!is.na(table3$W),]
	if (nrow(table3) > 0){
		to_return$'Mauchly\'s Test for Sphericity'=table3
		table2[,"p[GG]"] <- pf(table[,"F"], table2[,"GGe"]*table[,"DFn"],table2[,"GGe"]*table[,"DFd"], lower.tail=FALSE)
		table2[, "p[GG]<.05"] = ifelse(table2[, "p[GG]"]<.05,'*','')
		table2[,"p[HF]"] <- pf(table[,"F"], pmin(1, table2[,"HFe"])*table[,"DFn"],	pmin(1, table2[,"HFe"])*table[,"DFd"], lower.tail=FALSE)
		table2[, "p[HF]<.05"] = ifelse(table2[, "p[HF]"]<.05,'*','')
		table2=table2[!is.na(table2$GG),]
		to_return$'Sphericity Corrections'=table2
	}
	return(to_return)
}

ezANOVA_get_wide_lm<-
function(data, dv, sid, within, between){
	to_return = list()
	if(!is.null(within)){
		wide_formula = paste(paste(sid,paste(between,collapse='+'),sep='+'),paste(within,collapse='+'),sep='~')
		wide=cast(data, wide_formula, value = dv)
		to_return$idata=ldply(strsplit(names(wide)[!(names(wide) %in% c(between,sid))],'_'))
		wide_dv=data.matrix(wide[,!(names(wide) %in% c(sid,between))])
		names(to_return$idata)=within
		to_return$idesign_formula = paste('~',paste(within,collapse='*'),sep='')
	}else{
		wide=data
	}
	if(is.null(between)){
		lm_formula=paste('wide_dv~1',sep='')
	}else if(is.null(within)){
		lm_formula=paste(dv,'~',paste(between,collapse='*'),sep='')
	}else{
		lm_formula=paste('wide_dv~',paste(between,collapse='*'),sep='')
	}
	to_return$lm = lm(eval(parse(text=lm_formula)),wide)
	return(to_return)
}

ezANOVA_main <-
function(data, dv, sid, within, between){
	wide_lm = ezANOVA_get_wide_lm(data, dv, sid, within, between)
	if(is.null(within)){
		to_return = list()
		temp = as.data.frame(Anova(wide_lm$lm))
		names(temp) = c('SSn','DFn','F','p')
		temp$DFd = temp$D[length(temp$D)]
		temp$SSd = temp$S[length(temp$S)]
		temp$Effect = row.names(temp)
		row.names(temp) = 1:length(temp[,1])
		temp = temp[1:(length(temp[,1])-1),c(7,2,5,1,6,3,4)]
		temp$'p<.05'=ifelse(temp$p<.05,'*','')
		to_return$ANOVA = temp		
		to_return$'Levene\'s Test for Homogeneity of Variance' = ezANOVA_levene(wide_lm$lm)
		return(to_return)
	}else if(is.null(between)){
		return(suppressWarnings(ezANOVA_summary(Anova(wide_lm$lm,idata=wide_lm$idata,idesign=eval(parse(text=wide_lm$idesign_formula))))))
	}else{
		to_return = ezANOVA_summary(Anova(wide_lm$lm,idata=wide_lm$idata,idesign=eval(parse(text=wide_lm$idesign_formula))))
		no.within = ddply(
			data
			,.variables = structure(as.list(c(sid,between)),class = 'quoted')
			,function(x){
				each(mean)(x[,names(x)==dv])
			}
		)
		levene_formula = paste('mean~',paste(between,collapse='*'),sep='')
		to_return$'Levene\'s Test for Homogeneity of Variance (collapsing cells)' = ezANOVA_levene(eval(parse(text=levene_formula)),data=no.within)
		levene_formula = paste(dv,'~',paste(between,collapse='*'),sep='')
		cells = ddply(
			data
			, within
			,function(x){
				ezANOVA_levene(eval(parse(text=levene_formula)),data=x)
			}
		)
		to_return$'Levene\'s Test for Homogeneity of Variance (within cells)'=cells
		return(to_return)
	}
}


ezPerm_aov <-
function(data, aov_formula){
	this_aov = aov(
		formula(aov_formula)
		,data = data
	)
	f_list=llply(this_aov,function(x){summary(x)[[1]]$F})
	f = NULL
	for(i in f_list){
		f=c(f,i)
	}
	f=f[!is.na(f)]
	return(f)
}
