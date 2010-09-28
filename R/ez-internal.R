ezResample <-
function(
	data
	, dv
	, wid
	, within = NULL
	, between = NULL
	, resample_within = NULL
){
	if(!is.null(between)){
		ids = dlply(
			.data = data
			, .variables = between
			, .fun = function(x){
				to_return = sample(as.character(unique(x[,names(x)==as.character(wid)])),replace=T)
				return(to_return)
			}
		)
		ids = unlist(ids)
		names(ids) = NULL
	}else{
		ids = sample(as.character(unique(data[,names(data)==as.character(wid)])),replace=T)
	}
	id_list = list()
	for(i in 1:length(ids)){
		id_list[[i]] = list(num=i,this_id=ids[i])
	}
	resampled_data = ldply(
		.data = id_list
		, .fun = function(x){
			to_return = data[as.character(data[,names(data)==as.character(wid)])==x$this_id,]
			to_return[,names(to_return)==as.character(wid)] = x$num
			return(to_return)
		}
	)
	resampled_data[,names(resampled_data)==as.character(wid)] = factor(resampled_data[,names(resampled_data)==as.character(wid)])
	if(resample_within){
		to_return = ddply(
			.data = resampled_data
			, .variables = structure(as.list(c(wid,within)),class = 'quoted')
			, .fun = function(x){
	 			to_return = x[sample(1:nrow(x),nrow(x),replace=T),]
				return(to_return)
			}
		)
	}else{
		to_return = resampled_data
	}
	return(to_return)
}

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
	ANOVA = as.data.frame(table)
	to_return$ANOVA=ANOVA
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
function(data, dv, wid, within, between){
	to_return = list()
	if(!is.null(within)){
		for(this_within in within){
			old_levs = levels(data[,names(data)==this_within])
			new_levs = rep(NA,length=length(old_levs))
			temp = strsplit(old_levs,'_')
			for(i in 1:length(old_levs)){
				new_levs[i] = paste(temp[[i]],collapse='.')
			}
			levels(data[,names(data)==this_within]) = new_levs
		}
		wide_formula = paste(paste(wid,paste(between,collapse='+'),sep='+'),paste(within,collapse='+'),sep='~')
		wide=cast(data, wide_formula, value = dv)
		to_return$idata=ldply(strsplit(names(wide)[!(names(wide) %in% c(between,wid))],'_'))
		names(to_return$idata)=within
		for(this_within in within){
			to_return$idata[,names(to_return$idata)==this_within] = factor(to_return$idata[,names(to_return$idata)==this_within])
		}
		wide_dv=data.matrix(wide[,!(names(wide) %in% c(wid,between))])
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
function(data, dv, wid, within, between, observed, diff, reverse_diff){
	vars = as.character(c(dv,wid,between,within))
	for(var in vars){
		if(!(var %in% names(data))){
			stop(paste('"',var,'" is not a variable in the data frame provided.',sep=''))			
		}
	}
	if(is.null(within) & is.null(between)){
		stop('is.null(within) & is.null(between)\nYou must specify at least one independent variable.')
	}
	if(!is.data.frame(data)){
		stop('"data" must be a data frame.')
	}
	if(!is.numeric(data[,names(data)==dv])){
		stop('"dv" must be numeric.')
	}
	if(!is.factor(data[,names(data)==wid])){
		warning(paste('Converting "',wid,'" to factor for ANOVA.',sep=''),call.=FALSE)
		data[,names(data)==wid]=factor(data[,names(data)==wid])
	}else{
		if(length(unique(data[,names(data)==wid]))!=length(levels(data[,names(data)==wid]))){
			warning(paste('You have removed one or more Ss from the analysis. Refactoring "',wid,'" for ANOVA.',sep=''),call.=FALSE)
			data[,names(data)==wid]=factor(data[,names(data)==wid])
		}
	}
	vars = as.character(c(between,within,diff))
	for(var in vars){
		if(!is.factor(data[,names(data)==var])){
			warning(paste('Converting "',var,'" to factor for ANOVA.',sep=''),call.=FALSE)
			data[,names(data)==var]=factor(data[,names(data)==var])
		}
		if(length(unique(data[,names(data)==var]))!=length(levels(data[,names(data)==var]))){
			warning(paste('You have removed one or more levels from variable "',var,'". Refactoring for ANOVA.',sep=''),call.=FALSE)
			data[,names(data)==var]=factor(data[,names(data)==var])
		}
		if(length(levels(data[,names(data)==var]))==1){
			stop(paste('"',var,'" has only one level."',sep=''))			
		}
	}
	if(!is.null(diff)){
		temp <- ddply(
			idata.frame(data)
			,structure(as.list(c(wid,diff)),class = 'quoted')
			,function(x){
				to_return = 0
				return(to_return)
			}
		)
		if(!all(as.data.frame(table(temp[,names(temp) %in% c(wid,within)]))$Freq==2)){
			stop(paste('Variable supplied to "diff" ("',as.character(diff),'") does not appear to be a within variable.',sep=''))
		}
	}
	if(!is.null(diff)){
		data[,names(data)==as.character(diff)] = factor(data[,names(data)==as.character(diff)])
		if(length(unique(data[,names(data)==as.character(diff)]))!=2){
			stop('The column passed as argument "diff" must have precisely 2 levels.')
		}
		if(reverse_diff){
			data[,names(data)==as.character(diff)] = factor(data[,names(data)==as.character(diff)],levels=rev(levels(data[,names(data)==as.character(diff)])))
		}
	}
	temp = idata.frame(cbind(data,ezDV=data[,names(data) == as.character(dv)]))
	data <- ddply(
		temp
		,structure(as.list(c(wid,between,within,diff)),class = 'quoted')
		,function(x){
			to_return = mean(x$ezDV)
			names(to_return) = as.character(dv)
			return(to_return)
		}
	)
	if(any(is.na(data[,names(data)==as.character(dv)]))){
		stop('One or more cells returned NA when aggregated to a mean. Check your data.')
	}
	if(is.null(diff)){
		if(!all(as.data.frame(table(data[,names(data) %in% c(wid,within)]))$Freq==1)){
			stop('One or more cells is missing data.')
		}
	}else{
		if(!all(as.data.frame(table(data[,names(data) %in% c(wid,within,diff)]))$Freq==1)){
			#print(summary(as.data.frame(table(data[,names(data) %in% c(wid,within,diff)]))))
			stop('One or more cells is missing data.')
		}		
	}
	if(!is.null(between)){
		if(any(as.data.frame(table(data[,names(data) %in% c(between)]))$Freq==0)){
			stop('One or more cells is missing data.')
		}
	}
	if(!is.null(diff)){
		warning(paste('Collapsing "',as.character(diff),'" to a difference score ("',levels(data[,names(data)==as.character(diff)])[2],'"-"',levels(data[,names(data)==as.character(diff)])[1],'") prior to computing statistics.',sep=''),call.=FALSE)
		temp = idata.frame(cbind(data,ezDV=data[,names(data) == as.character(dv)]))
		data <- ddply(
			temp
			,structure(as.list(c(wid,within,between)),class = 'quoted')
			,function(x){
				to_return = diff(x$ezDV)
				names(to_return) = as.character(dv)
				return(to_return)
			}
		)
		temp = names(within)
		temp = temp[!(within %in% diff)]
		within = within[!(within %in% diff)]
		names(within) = temp
	}
	wide_lm = ezANOVA_get_wide_lm(data, dv, wid, within, between)
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
	}else{
		to_return = NULL
		try(to_return<-suppressWarnings(ezANOVA_summary(Anova(wide_lm$lm,idata=wide_lm$idata,idesign=eval(parse(text=wide_lm$idesign_formula))))),silent=TRUE)
		if(is.null(to_return)){
			warning('Anova() failed for some reason (maybe there are too few Ss?), trying aov()...',call.=FALSE)
			to_return=list(ANOVA=ezANOVA_aov(data, dv, wid, within, between))
		}
	}
	to_return$data = data
	return(to_return)
}

ezANOVA_aov <-
function(data, dv, wid, within, between){
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
				,as.character(wid)
				,')'
				,sep = ''
			)
			,paste(
				'+Error('
				,as.character(wid)
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

progress_time = function () 
{
	n <- 0
	txt <- NULL
	list(
		init = function(x) {
	    	txt <<- init_progress_time(max = x)
	    	setTxtProgressBar(txt, 0)
		}
		, step = function() {
	    	n <<- n + 1
	    	setTxtProgressBar(txt, n)
		}
		, term = function() close(txt)
	)
}


init_progress_time = function (min = 0, max = 1, initial = 0){
	.start <- proc.time()[3]
    .val <- initial
    .killed <- FALSE
    .nb <- 0L
    .pc <- 0L
    width <- getOption("width")
    width <- width - 10L - 55L
    width <- trunc(width)
    if (max <= min) 
        stop("must have max > min")
    up <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
        .val <<- value
		minutes = TRUE
		if(value>0){
			time_taken = proc.time()[3]-.start
			time_per_value = time_taken/value
			time_left = round((max-value)*time_per_value,0)
			if(time_left>60){
				time_left = round((max-value)*time_per_value/60,0)
			}else{
				minutes = FALSE
			}
		}
        nb <- round(width * (value - min)/(max - min))
        pc <- round(100 * (value - min)/(max - min))
        if (nb == .nb && pc == .pc) 
            return()
        cat(paste(c("\r  |", rep.int(" ", width + 6)), collapse = ""))
		if(minutes){
	        cat(paste(c("\r  |", rep.int('=', nb), rep.int(" ", 
	            (width - nb)), sprintf("| %3d%%", pc),' Approximately ',time_left,' minutes remaining until completion.',rep(' ',4-nchar(as.character(time_left)))), collapse = ""))
		}else{
			cat(paste(c("\r  |", rep.int('=', nb), rep.int(" ", 
	            (width - nb)), sprintf("| %3d%%", pc),' Approximately ',time_left,' seconds remaining until completion.',rep(' ',4-nchar(as.character(time_left)))), collapse = ""))
		}
        flush.console()
        .nb <<- nb
        .pc <<- pc
    }
    getVal <- function() .val
    kill <- function() if (!.killed) {
		cat(paste(c("\r  |", rep.int('=', .nb), rep.int(" ", 
            (width - .nb)), sprintf("| %3d%%", .pc),' Complete.',rep(' ',45)), collapse = ""))
        cat("\n")
        flush.console()
        .killed <<- TRUE
    }
    if (initial > min) 
        up(initial)
    structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}
