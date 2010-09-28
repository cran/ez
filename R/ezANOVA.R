ezANOVA <-
function(
	data
	, dv
	, wid
	, within = NULL
	, between = NULL
	, observed = NULL
	, diff = NULL
	, reverse_diff = FALSE
	, detailed = FALSE
){
	to_return = ezANOVA_main(data,dv,wid,within,between,observed,diff,reverse_diff)

	########
	# Compute effect size
	########
	if(!is.null(observed)){
		obs = rep(F,nrow(to_return$ANOVA))
		for(i in as.character(observed)){
			obs = obs | str_detect(to_return$ANOVA$Effect,i)
		}
		obs_SSn1 = sum(to_return$ANOVA$SSn*obs)
		obs_SSn2 = to_return$ANOVA$SSn*obs
	}else{
		obs_SSn1 = 0
		obs_SSn2 = 0
	}
	to_return$ANOVA$ges = to_return$ANOVA$SSn/(to_return$ANOVA$SSn+sum(unique(to_return$ANOVA$SSd))+obs_SSn1-obs_SSn2)
	
	########
	# Final clean-up
	########

	#remove the data from to_return
	temp = names(to_return)
	temp = temp[temp!='data']
	to_return = to_return[names(to_return)!='data']
	names(to_return) = temp

	#if necessary, remove extra columns and the Intercept row from the anova
	if(!detailed){
		start = ifelse(is.null(within),1,2)
		to_return$ANOVA = to_return$ANOVA[,names(to_return$ANOVA) %in% c('Effect','DFn','DFd','F','p','p<.05','ges')]
		to_return$ANOVA = to_return$ANOVA[start:nrow(to_return$ANOVA),]
	}

	#all done!
	return(to_return)
}

