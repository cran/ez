ezPredict <-
function(
	fit
	, to_predict = NULL
	, numeric_res = 1e3
){
	data = attr(fit,'frame')
	vars = attr(attr(data,'terms'),'variables')
	dv = as.character(vars[2])
	if(is.null(to_predict)){
		fixed = rep(NA,length(vars)-3)
		temp = list()
		j = 1
		for(i in 3:length(vars)){
			fixed[i-2] = as.character(vars[i])
			this_fixed_data = data[,names(data)==fixed[i-2]]
			if(is.numeric(this_fixed_data)){
				temp[[j]] = seq(
					min(this_fixed_data)
					, max(this_fixed_data)
					, length.out=numeric_res
				)
			}else{
				temp[[j]] = unique(this_fixed_data)
			}
			j = j + 1
		}
		to_return = expand.grid(temp)
		names(to_return) = as.character(fixed)
	}else{
		to_return = to_predict
	}
	to_return$ezDV = 0
	names(to_return)[ncol(to_return)] = dv
	mm = model.matrix(terms(fit),to_return)
	to_return = to_return[,1:(ncol(to_return)-1)]
	to_return$value = mm %*% fixef(fit)
	#print(mm)
	vf = vcov(fit)
	#print(vf)
	tc = Matrix::tcrossprod(vf,mm)
	to_return$var = Matrix::diag(mm %*% tc)
	return(to_return)
}