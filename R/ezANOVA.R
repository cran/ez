ezANOVA <-
function(
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
){
	if(is.null(within) & is.null(between)){
		stop('is.null(within) & is.null(between)\nYou must specify at least one independent variable.')
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
	data=data[,names(data) %in% c(within,between,sid,dv)]
	return(ezANOVA_main(data,dv,sid,within,between))
}

