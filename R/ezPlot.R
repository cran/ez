ezPlot <-
function (
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, x
	, do_lines = TRUE
	, do_bars = TRUE
	, bar_width = NULL
	, bar_size = NULL
	, split = NULL
	, row = NULL
	, col = NULL
	, to_numeric = NULL
	, x_lab = NULL
	, y_lab = NULL
	, split_lab = NULL
){
	from_ezStats_main = ezStats_main(data,dv,sid,within,between)
	data = from_ezStats_main$Descriptives
	this_ANOVA = from_ezStats_main$ANOVA
	if(!(x %in% within) & !(x %in% between) ){
		stop('"x" not listed in "within" or "between".')
	}
	if(!is.null(split)){
		if(!(split %in% within) & !(split %in% between) ){
			stop('"split" not listed in "within" or "between".')
		}
	}
	if(!is.null(row)){
		if(!(row %in% within) & !(row %in% between) ){
			stop('"row" not listed in "within" or "between".')
		}
	}
	if(!is.null(col)){
		if(!(col %in% within) & !(col %in% between) ){
			stop('"col" not listed in "within" or "between".')
		}
	}
	if(!is.null(to_numeric)){
		if(!(to_numeric %in% within) & !(to_numeric %in% between) ){
			stop('"to_numeric" not listed in "within" or "between".')
		}
	}
	if(!is.logical(do_lines)){
		stop('"do_lines" must be either TRUE or FALSE.')
	}
	if(!is.logical(do_bars)){
		stop('"do_bars" must be either TRUE or FALSE.')
	}
	if(!is.null(bar_width)){
		if(!is.numeric(bar_width)){
			stop('"bar_width" must be numeric.')
		}else{
			if(bar_width<=0){
				stop('"bar_width" must be > 0.')
			}
		}
	}else{
		if(!is.numeric(data[,names(data)==x])){
			bar_width = .5
		}
	}
	if(!is.null(bar_size)){
		if(!is.numeric(bar_size)){
			stop('"bar_size" must be numeric.')
		}else{
			if(bar_size<=0){
				stop('"bar_size" must be > 0.')
			}
		}
	}
	if(is.null(bar_size)){
		bar_size = data$FLSD
	}
	data$ymin = data$Mean-bar_size/2
	data$ymax = data$Mean+bar_size/2
	for(i in to_numeric){
		data[,names(data) == i] = as.numeric(as.character(data[,names(data) == i]))
	}
	names(data)[names(data) == x] = 'x'
	if(!is.null(split)){
		names(data)[names(data) == split] = 'split'
	}
	if(!is.null(row)){
		names(data)[names(data) == row] = 'row'
	}
	if(!is.null(col)){
		names(data)[names(data) == col] = 'col'
	}
	p = ggplot(
		data = data
		,aes(
			y = Mean
			,x = x
		)
	)
	if(!is.null(split)){
		p = p+geom_point(
			aes(
				colour = split
				,shape = split
			)
		)
		if(!is.null(split_lab)){
			p = p+labs(colour = split_lab,shape = split_lab)
		}
		if(do_lines){
			p = p+geom_line(
				aes(
					colour = split
					,linetype = split
					,x = as.numeric(x)
				)
			)
			if(!is.null(split_lab)){
				p = p+labs(linetype = split_lab)
			}
		}
		if(do_bars){
			p = p+geom_errorbar(
				aes(
					colour = split
					,ymin = ymin
					,ymax = ymax
				)
				,linetype = 1
				,legend = FALSE
				,width = bar_width
			)
		}
	}else{
		p = p+geom_point()
		if(do_lines){
			p = p+geom_line(aes(x = as.numeric(x)))
		}
		if(do_bars){
			p = p+geom_errorbar(
				aes(
					ymin = ymin
					,ymax = ymax
				)
				,linetype = 1
				,legend = FALSE
				,width = bar_width
			)
		}
	}
	if(!is.null(row)){
		if(!is.null(col)){
			p = p+facet_grid(row~col)
		}else{
			p = p+facet_grid(row~.)
		}
	}else{
		if(!is.null(col)){
			p = p+facet_grid(.~col)
		}
	}
	if(!is.null(x_lab)){
		p = p+labs(x = paste('\n',x_lab,sep = ''))
	}
	if(!is.null(y_lab)){
		p = p+labs(y = paste(y_lab,'\n',sep = ''))
	}
	return(p)
}

