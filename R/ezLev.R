ezLev<-
function(
	x
	,new_order
){
	for(i in rev(new_order)){
		x=relevel(x,ref=i)
	}
	return(x)
}
