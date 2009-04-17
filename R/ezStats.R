ezStats <-
function (
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
){
	if(is.null(within) & is.null(between)){
		stop('is.null(within) & is.null(between)\nYou must specify at least one independent variable.')
	}else{
		if(!is.null(within) & !is.null(between)){
			warning('Mixed within-and-between-Ss effect requested; FLSD is only appropriate for within-Ss comparisons (see warning in ?ezStats or ?ezPlot).',call.=FALSE)
		}
	}
	if(!is.data.frame(data)){
		stop('"data" must be a data frame.')
	}
	if(!is.numeric(data[,names(data)==dv])){
		stop('"dv" must be numeric.')
	}
	if(!is.factor(data[,names(data)==sid])){
		warning(paste('Converting "',sid,'" to factor for ANOVA.',sep=''),call.=FALSE)
		data[,names(data)==sid]=factor(data[,names(data)==sid])
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
	N = mean(N[,length(N)])
	data <- ddply(
		data
		,structure(as.list(c(sid,between,within)),class = 'quoted')
		,function(x){
			mean(x[,names(x) == as.character(dv)])
		}
	)
	this_ANOVA = ezANOVA(
		data = data
		, within = within
		, between = between
		, sid = sid
		, dv = .(V1)
	)$ANOVA
	DFd = this_ANOVA$DFd[length(this_ANOVA$DFd)]
	MSd = this_ANOVA$SSd[length(this_ANOVA$SSd)]/DFd
	Tcrit = qt(0.975,DFd)
	CI = Tcrit * sqrt(MSd/N)
	FLSD = sqrt(2) * CI
	.variables = structure(as.list(c(between,within)),class = 'quoted')
	data <- ddply(
		data
		,.variables
		,function(x){
			N = length(x$V1)
			Mean = mean(x$V1)
			SD = sd(x$V1)
			return(c(N = N, Mean = Mean, SD = SD))
		}
	)
	data$FLSD = FLSD
	return(
		list(
			Descriptives = data
			, ANOVA = this_ANOVA
		)
	)
}

