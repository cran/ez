ezPlotBoot <-
function(
	from_ezBoot
	, x
	, split = NULL
	, row = NULL
	, col = NULL
	, do_lines = TRUE
	, bar_width = NULL
	, to_numeric = NULL
	, x_lab = NULL
	, y_lab = NULL
	, split_lab = NULL
	, levels = NULL
	, diff = NULL
	, reverse_diff = FALSE
	, row_y_free = FALSE
){
	if(!is.logical(do_lines)){
		stop('"do_lines" must be either TRUE or FALSE.')
	}
	if(!is.null(levels)){
		for(i in 1:length(levels)){
			this_iv = names(levels)[i]
			from_ezBoot$cells[,names(from_ezBoot$cells)==this_iv] = factor(from_ezBoot$cells[,names(from_ezBoot$cells)==this_iv])
			if('new_order' %in% names(levels[[i]])){
				from_ezBoot$cells[,names(from_ezBoot$cells)==this_iv] = factor(from_ezBoot$cells[,names(from_ezBoot$cells)==this_iv],levels=levels[[i]]$new_order)
			}
			if('new_names' %in% names(levels[[i]])){
				levels(from_ezBoot$cells[,names(from_ezBoot$cells)==this_iv]) = levels[[i]]$new_names
			}
			from_ezBoot$boots[,names(from_ezBoot$boots)==this_iv] = factor(from_ezBoot$boots[,names(from_ezBoot$boots)==this_iv])
			if('new_order' %in% names(levels[[i]])){
				from_ezBoot$boots[,names(from_ezBoot$boots)==this_iv] = factor(from_ezBoot$boots[,names(from_ezBoot$boots)==this_iv],levels=levels[[i]]$new_order)
			}
			if('new_names' %in% names(levels[[i]])){
				levels(from_ezBoot$boots[,names(from_ezBoot$boots)==this_iv]) = levels[[i]]$new_names
			}
		}
	}
	cat('ezPlotBoot: Collapsing cells to requested design...')
	cells = ddply(
		.data = idata.frame(from_ezBoot$cells)
		, .variables = structure(as.list(c(x,split,row,col,diff)),class = 'quoted')
		, .fun = function(x){
			to_return = data.frame(
				value = mean(x$value)
			)
			return(to_return)
		}
	)
	cat('\nezPlotBoot: Collapsing boots to requested design...')
	boots = ddply(
		.data = idata.frame(from_ezBoot$boots)
		, .variables = structure(as.list(c(x,split,row,col,diff,expression(iteration))),class = 'quoted')
		, .fun = function(x){
			to_return = data.frame(
				value = mean(x$value)
			)
			return(to_return)
		}
	)
	if(!is.null(diff)){
		if(reverse_diff){
			cells[,names(cells)==as.character(diff)] = factor(
				cells[,names(cells)==as.character(diff)]
				, levels = rev(levels(cells[,names(cells)==as.character(diff)]))
			)
			boots[,names(boots)==as.character(diff)] = factor(
				boots[,names(boots)==as.character(diff)]
				, levels = rev(levels(boots[,names(boots)==as.character(diff)]))
			)
		}
		cat('\nezPlotBoot: Computing requested difference score within cells...')
		cells = ddply(
			.data = cells
			, .variables = structure(as.list(c(x,split,row,col)),class = 'quoted')
			, .fun = function(x){
				to_return = data.frame(
					value = x$value[x[,names(x)==as.character(diff)]==(levels(x[,names(x)==as.character(diff)])[1])] - x$value[x[,names(x)==as.character(diff)]==(levels(x[,names(x)==as.character(diff)])[2])]
				)
				return(to_return)
			}
		)
		cat('\nezPlotBoot: Computing requested difference score within boots...')
		boots = ddply(
			.data = boots
			, .variables = structure(as.list(c(x,split,row,col,expression(iteration))),class = 'quoted')
			, .fun = function(x){
				to_return = data.frame(
					value = x$value[x[,names(x)==as.character(diff)]==(levels(x[,names(x)==as.character(diff)])[1])] - x$value[x[,names(x)==as.character(diff)]==(levels(x[,names(x)==as.character(diff)])[2])]
				)
				return(to_return)
			}
		)	
	}
	cat('\nezPlotBoot: Computing confidence intervals...')
	boot_stats = ddply(
		.data = idata.frame(boots)
		, .variables = structure(as.list(c(x,split,row,col)),class = 'quoted')
		, .fun = function(x){
			to_return = data.frame(
				lo = quantile(x$value,.025)
				, hi = quantile(x$value,.975)
			)
			return(to_return)
		}
	)
	cat('\nezPlotBoot: Building plot...')
	names(cells)[names(cells)==as.character(x)] = 'x'
	names(boot_stats)[names(boot_stats)==as.character(x)] = 'x'
	if(!is.null(split)){
		names(cells)[names(cells)==as.character(split)] = 'split'
		names(boot_stats)[names(boot_stats)==as.character(split)] = 'split'
	}
	if(!is.null(row)){
		names(cells)[names(cells)==as.character(row)] = 'row'
		names(boot_stats)[names(boot_stats)==as.character(row)] = 'row'
	}
	if(!is.null(col)){
		names(cells)[names(cells)==as.character(col)] = 'col'
		names(boot_stats)[names(boot_stats)==as.character(col)] = 'col'
	}
	p = ggplot(
		data = cells
		, mapping = aes(
			x = x
		)
	)
	if(!is.null(split)){
		p = p+geom_point(
			aes(
				colour = split
				, shape = split
				, y = value
			)
			, alpha = .8
		)
		if(!is.null(split_lab)){
			p = p+labs(colour = split_lab,shape = split_lab)
		}
		if(do_lines){
			p = p+geom_line(
				aes(
					colour = split
					, linetype = split
					, x = as.numeric(x)
					, y = value
				)
				, alpha = .8
			)
			if(!is.null(split_lab)){
				p = p+labs(linetype = split_lab)
			}
		}
		p = p+geom_errorbar(
			data = boot_stats
			, mapping = aes(
				colour = split
				, ymin = lo
				, ymax = hi
			)
			, linetype = 1
			, legend = FALSE
			, width = bar_width
			, alpha = .5
		)
	}else{
		p = p+geom_point(
			mapping = aes(
				y = value
			)
		)
		if(do_lines){
			p = p+geom_line(
				mapping = aes(
					x = as.numeric(x)
					, y = value
				)
			)
		}
		p = p+geom_errorbar(
			data = boot_stats
			, mapping = aes(
				, ymin = lo
				, ymax = hi
			)
			, linetype = 1
			, legend = FALSE
			, width = bar_width
			, alpha = .5
		)
	}
	if(!is.null(row)){
		if(!is.null(col)){
			if(row_y_free){
				p = p+facet_grid(row~col,scales='free_y')
			}else{
				p = p+facet_grid(row~col)
			}
		}else{
			if(row_y_free){
				p = p+facet_grid(row~.,scales='free_y')
			}else{
				p = p+facet_grid(row~.)
			}
		}
	}else{
		if(!is.null(col)){
			p = p+facet_grid(.~col)
		}
	}
	if(!is.null(x_lab)){
		p = p+labs(x = x_lab)
	}
	if(!is.null(y_lab)){
		p = p+labs(y = y_lab)
	}
	to_return = list()
	to_return$plot = p
	to_return$cells = cells
	to_return$boots = boots
	to_return$boot_stats = boot_stats
	cat('\nezPlotBoot: Done.')
	return(to_return)
}

