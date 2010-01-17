ezANOVA_levene <-
function (y) {
	form <- y
	mf <- model.frame(form)
	if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]])))) stop("Levene's test is not appropriate with quantitative explanatory variables.")
	y <- mf[,1]
	if(dim(mf)[2]==2) {
		group <- mf[,2]
	}else {
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
	table <- data.frame(matrix(0, nterms, 9))
	table2 <- data.frame(matrix(0, nterms, 7))
	table3 <- data.frame(matrix(0, nterms, 4))
	table3[,1] <- table2[,1] <- table[,1] <- object$terms
	colnames(table) <- c("Effect","DFn", "DFd", "SSn", "SSd", "F", "p", "p<.05", "pes")
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
		table[term, "pes"] <- table[term, "SSn"] / (table[term, "SSn"] + table[term, "SSd"])
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
	}else{
		to_return = NULL
		try(to_return<-suppressWarnings(ezANOVA_summary(Anova(wide_lm$lm,idata=wide_lm$idata,idesign=eval(parse(text=wide_lm$idesign_formula))))),silent=TRUE)
		if(is.null(to_return)){
			warning('Too few Ss for Anova(), reverting to aov(). See "Warning" section of the help on ezANOVA.',call.=FALSE)
			to_return=list(ANOVA=ezANOVA_aov(data, dv, sid, within, between))
		}
		return(to_return)
	}
}

ezANOVA_aov <-
function(data, dv, sid, within, between){
	aov_formula = paste(
		as.character(dv)
		,'~'
		,paste(as.character(between),collapse = '*')
		,ifelse(is.null(between),'',ifelse(is.null(within),'','*'))
		,paste(as.character(within),collapse = '*')
		,ifelse(
			is.null(within)
			,paste(
				'+Error('
				,as.character(sid)
				,')'
				,sep = ''
			)
			,paste(
				'+Error('
				,as.character(sid)
				,'/('
				,paste(as.character(within),collapse = '*')
				,'))'
				,sep = ''
			)
		)
		,sep = ''
	)	
	this_aov = aov(
		formula(aov_formula)
		,data = data
	)
	ANOVA = NULL
	for(x in summary(this_aov)){
		if(length(x)==1){
			x=x[[1]]
		}
		for(row in 1:length(x[,1])){
			if(!is.na(x$P[row])){
				ANOVA = rbind(
					ANOVA
					, data.frame(
						Effect=strsplit(row.names(x)[row],' ')[[1]][1]
						, DFn=x$D[row]
						, DFd=x$D[length(x$D)]
						, SSn=x$S[row]
						, SSd=x$S[length(x$S)]
						, F=x$F[row]
						, p=x$P[row]
					)
				)
			}
		}
	}
	ANOVA$'p<.05'=ifelse(ANOVA$p<.05,'*','')
	ANOVA$pes = ANOVA$SSn/(ANOVA$SSn+ANOVA$SSd)
	return(ANOVA)
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

ezStats_main <-
function (
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, between_full = NULL
	, collapse_within = FALSE
){
	if(is.null(within) & is.null(between)){
		stop('is.null(within) & is.null(between)\nYou must specify at least one independent variable.')
	}else{
		if(!is.null(within) & !is.null(between)){
			if(!collapse_within){
				warning('Mixed within-and-between-Ss effect requested; FLSD is only appropriate for within-Ss comparisons (see warning in ?ezStats or ?ezPlot).',call.=FALSE)
			}
		}
	}
	if(!is.data.frame(data)){
		stop('"data" must be a data frame.')
	}
	if(!is.numeric(data[,names(data)==dv])){
		stop('"dv" must be numeric.')
	}
	vars = as.character(c(dv,sid,between,within))
	for(var in vars){
		if(!(var %in% names(data))){
			stop(paste('"',var,'" is not a variable in "',data,'".',sep=''))			
		}
	}
	if(!is.factor(data[,names(data)==sid])){
		warning(paste('Converting "',sid,'" to factor for ANOVA.',sep=''),call.=FALSE)
		data[,names(data)==sid]=factor(data[,names(data)==sid])
	}else{
		if(length(unique(data[,names(data)==sid]))!=length(levels(data[,names(data)==sid]))){
			warning(paste('You have removed one or more Ss from the analysis. Refactoring "',sid,'" for ANOVA.',sep=''),call.=FALSE)
			data[,names(data)==sid]=factor(data[,names(data)==sid])
		}
	}
	for(i in within){
		if(!is.factor(data[,names(data)==i])){
			warning(paste('Converting "',i,'" to factor for ANOVA.',sep=''),call.=FALSE)
			data[,names(data)==i]=factor(data[,names(data)==i])
		}
	}
	for(i in between){
		if(!is.factor(data[,names(data)==i])){
			warning(paste('Converting "',i,'" to factor for ANOVA.',sep=''),call.=FALSE)
			data[,names(data)==i]=factor(data[,names(data)==i])
		}
		levs = levels(data[,names(data)==i])
		if(length(levs)==1){
			stop(paste('Grouping variable "',i,'" has only one level.',sep=''))	
		}else{
			for(j in levs){
				if(!any(data[,names(data)==i]==j)){
					if(length(levs)==2){
						stop(paste('Group "',j,'" in "',i,'" has no members and there are only 2 groups specified by "',i,'".',sep=''))			
					}else{
						warning(paste('Group "',j,'" in "',i,'" has no members; removing group "',j,'" from the analysis.',sep=''),call.=FALSE)
						data[,names(data)==i]=factor(data[,names(data)==i])					
					}
				}
			}
		}
	}
	N = ddply(
		cbind(data,dummy = rep(1,length(data[,1])))
		,structure(as.list(c(.(dummy),between)),class = 'quoted')
		,function(x){
			to_return = length(unique(x[,names(x) == as.character(sid)]))
			names(to_return) = 'N'
			return(to_return)
		}
	)
	if(!all(N[,length(N)]==N[1,length(N)])){
		warning('Unbalanced groups. Mean N will be used in computation of FLSD')
		N = mean(N[,length(N)])
	}else{
		N = N[1,length(N)]
	}
	if(is.null(between_full)){
		temp_between = between
	}else{
		temp_between = between_full
	}
	this_ANOVA = ezANOVA(
		data = data
		, within = within
		, between = temp_between
		, sid = sid
		, dv = dv
		, collapse_within = collapse_within
	)$ANOVA
	DFd = this_ANOVA$DFd[length(this_ANOVA$DFd)]
	MSd = this_ANOVA$SSd[length(this_ANOVA$SSd)]/DFd
	Tcrit = qt(0.975,DFd)
	CI = Tcrit * sqrt(MSd/N)
	FLSD = sqrt(2) * CI
	data <- ddply(
		data
		,structure(as.list(c(sid,between,within)),class = 'quoted')
		,function(x){
			to_return = mean(x[,names(x) == as.character(dv)])
			names(to_return) = as.character(dv)
			return(to_return)
		}
	)
	data <- ddply(
		data
		,structure(as.list(c(between,within)),class = 'quoted')
		,function(x){
			N = length(x[,names(x) == as.character(dv)])
			Mean = mean(x[,names(x) == as.character(dv)])
			SD = sd(x[,names(x) == as.character(dv)])
			return(c(N = N, Mean = Mean, SD = SD))
		}
	)
	data$FLSD = FLSD
	return(data)
}

