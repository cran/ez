ezANOVA <-
function(
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, collapse_within = FALSE
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
	data <- ddply(
		data
		,structure(as.list(c(sid,between,within)),class = 'quoted')
		,function(x){
			to_return = mean(x[,names(x) == as.character(dv)])
			names(to_return) = as.character(dv)
			return(to_return)
		}
	)
	if(any(is.na(data[,names(data)==as.character(dv)]))){
		stop('One or more cells returned NA when aggregated to a mean. Check your data.')
	}
	if(!all(as.data.frame(table(data[,names(data) %in% c(sid,within)]))$Freq==1)){
		stop('One or more cells is missing data.')
	}
	if(collapse_within){
		if(is.null(within)){
			stop('When setting collapse_within=TRUE, you must provide a within-Ss variable.')
		}else{
			if(length(unique(data[,names(data)==as.character(within)]))!=2){
				stop('When setting collapse_within=TRUE, you must provide a within-Ss variable with precisely 2 levels.')
			}else{
				warning('Collapsing the within-Ss variable to a difference score prior to computing statistics.',call.=FALSE)
				data <- ddply(
					data
					,structure(as.list(c(sid,between)),class = 'quoted')
					,function(x){
						to_return = diff(x[,names(x) == as.character(dv)])
						names(to_return) = as.character(dv)
						return(to_return)
					}
				)
				within = NULL
			}
		}
	}
	return(ezANOVA_main(data,dv,sid,within,between))
}

