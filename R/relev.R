relev<-
function(
	x
	,levs
){
	for(i in rev(levs)){
		x=relevel(x,ref=i)
	}
	return(x)
}
