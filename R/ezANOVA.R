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
){
	to_return = ezANOVA_main(data,dv,wid,within,between,observed,diff,reverse_diff)
	temp = names(to_return)
	temp = temp[temp!='data']
	to_return = to_return[names(to_return)!='data']
	names(to_return) = temp
	return(to_return)
}

