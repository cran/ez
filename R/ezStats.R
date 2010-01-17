ezStats <-
function (
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, between_full = NULL
	, collapse_within = FALSE
){
	ezStats_main(data,dv,sid,within,between,between_full,collapse_within)
}

