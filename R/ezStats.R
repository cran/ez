ezStats <-
function (
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, between_full = NULL
){
	from_ezStats_main = ezStats_main(data,dv,sid,within,between,between_full)
	return(from_ezStats_main$Descriptives)
}

