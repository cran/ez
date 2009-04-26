ezPerm <-
function(
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, perms
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
	if(!is.numeric(perms)){
		stop('"perms" must be numeric.')
	}else{
		if(perms<0){
			stop('"perms" must be >= 0.')			
		}else{
			if(perms%%1){
				stop('"perms" must be an integer.')
			}
		}
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
	cat('Permutation test progress:\n')
	start = proc.time()[1]
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
	obs = ezPerm_aov(data,aov_formula)
	sim = matrix(NA,length(obs),perms)
	sim_data=data
	if(!is.null(between)){
		group_info = ddply(
			sim_data
			,.variables = structure(as.list(c(sid,between)),class='quoted')
			,function(x){
				x[1,]
			}
		)
	}
	pb = txtProgressBar(max=perms,style=3)
	for(this_perm in 1:perms){
		for(this_within in within){
			sim_data = ddply(
				sim_data
				,.variables = structure(as.list(c(sid,structure(within[within!=this_within],class='quoted'))),class='quoted')
				,function(x){
					x[,names(x)==this_within] = sample(x[,names(x)==this_within])
					return(x)
				}
			)
		}
		for(this_between in between){
			group_info[,names(group_info)==this_between]=sample(group_info[,names(group_info)==this_between])
			for(this_sid in sim_data[,names(sim_data)==sid]){
				sim_data[sim_data[,names(sim_data)==sid]==this_sid,names(sim_data)==this_between] = group_info[group_info[,names(group_info)==sid]==this_sid,names(group_info)==this_between]
			}
		}
		sim[,this_perm] = ezPerm_aov(sim_data,aov_formula)
		setTxtProgressBar(pb, this_perm)
	}
	Effect = ezANOVA_main(data,dv,sid,within,between)$ANOVA$Effect
	if(is.null(between)){
		Effect = Effect[2:length(Effect)]
	}
	perm_test = data.frame(Effect=Effect)
	perm_test$p = rowMeans(sim>=obs)
	perm_test$'p<.05' = ifelse(perm_test$p<.05,'*','')
	close(pb)
	return(list('Permutation Test'=perm_test,'Test Duration' = as.vector(proc.time()[1] - start)))
}

