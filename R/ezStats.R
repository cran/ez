ezStats <-
function (
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
){
	from_ezStats_main = ezStats_main(data,dv,sid,within,between)
	return(from_ezStats_main$Descriptives)
}

