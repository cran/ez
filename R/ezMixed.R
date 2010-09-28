ezMixed <-
function(
	data
	, dv
	, random
	, fixed
	, family = gaussian
	, alarm = TRUE
	, depth = 0
	, return_models = FALSE
	, return_anovas = FALSE
){
	start = proc.time()[3]
	LLR_lmer = function(m0,m1){
		LLRa = AIC(m0) - AIC(m1)
		LLRb = BIC(m0) - BIC(m1)
		LLRa = log(exp(1),base=10)*LLRa
		LLRb = log(exp(1),base=10)*LLRb
		return(c(LLRa,LLRb))
	}
	vars = as.character(c(dv,random,fixed))
	for(var in vars){
		if(!(var %in% names(data))){
			stop(paste('"',var,'" is not a variable in the data frame provided.',sep=''))			
		}
	}
	if(!is.data.frame(data)){
		stop('"data" must be a data frame.')
	}
	if(!is.numeric(data[,names(data)==dv])){
		stop('"dv" must be numeric.')
	}
	fixed = as.character(fixed)
	randoms = NULL
	for(i in as.character(random)){
		randoms = paste(randoms,'+(1|',i,')',sep='')
	}
	full_formula = paste(
		as.character(dv)
		, '~'
		, paste(fixed,collapse='*')
		, randoms
	)
	from_terms = terms(eval(parse(text=full_formula)))
	term_labels = attr(from_terms,'term.labels')
	term_labels = term_labels[!str_detect(term_labels,'1')]
	if(depth>0){
		term_labels = term_labels[laply((strsplit(term_labels,':')),length)<=depth]
	}
	to_return = list()
	to_return$summary = data.frame(
		effect = term_labels
		, p = NA
		, LLRa = NA
		, LLRb = NA
	)
	list_index = 1
	this_baseline_formula = paste(
		as.character(dv)
		, '~(1|'
		, as.character(random)
		, ')'
		, sep = ''
	)
	fit_baseline = lmer(
		formula = eval(parse(text=this_baseline_formula))
		, family = family
		, data = data
		, REML = FALSE
	)
	for(i in 1:length(term_labels)){
		effect_split = strsplit(term_labels[i],':')[[1]]
		levels = length(effect_split)
		if(levels==1){
			this_comparison_formula = paste(
				as.character(dv)
				, '~(1|'
				, as.character(random)
				, ')+'
				, effect_split
				, sep = ''
			)
			fit_comparison = lmer(
				formula = eval(parse(text=this_comparison_formula))
				, family = family
				, data = data
				, REML = FALSE
			)
		}else{
			this_from_terms = terms(
				eval(
					parse(
						text = paste(
							'1~'
							, paste(
								effect_split
								, collapse = '*'
							)
						)
					)
				)
			)
			this_term_labels = attr(this_from_terms,'term.labels')
			this_comparison_formula = paste(
				as.character(dv)
				, '~(1|'
				, as.character(random)
				, ')+'
				, paste(
					this_term_labels
					, collapse = '+'
					, sep = ''
				)
				, sep = ''
			)
			this_baseline_formula = paste(
				as.character(dv)
				, '~(1|'
				, as.character(random)
				, ')+'
				, paste(
					this_term_labels[1:(length(this_term_labels)-1)]
					, collapse = '+'
					, sep = ''
				)
				, sep = ''
			)
			fit_comparison = lmer(
				formula = eval(parse(text=this_comparison_formula))
				, family = family
				, data = data
				, REML = FALSE
			)
			fit_baseline = lmer(
				formula = eval(parse(text=this_baseline_formula))
				, family = family
				, data = data
				, REML = FALSE
			)
		}
		if(return_models|return_anovas){
			to_return[[i+1]] = list()
			names(to_return)[i+1] = term_labels[i]			
		}
		if(return_models){
			to_return[[i+1]]$fit = fit_comparison
		}
		this_anova = anova(fit_baseline,fit_comparison)
		if(return_anovas){
			to_return[[i+1]]$anova = this_anova
			attr(to_return[[i+1]]$anova,'heading') = c("Models:",this_baseline_formula,this_comparison_formula)
		}
		to_return$summary$p[i] = this_anova$P[2]
		LLRs = LLR_lmer(fit_baseline,fit_comparison)
		to_return$summary$LLRa[i] = LLRs[1]
		to_return$summary$LLRb[i] = LLRs[2]
		cat(c(
			term_labels[i],' -> '
			,'p =',signif(to_return$summary$p[i],2)
			,', LLRa =',signif(to_return$summary$LLRa[i],2)
			,', LLRb =',signif(to_return$summary$LLRb[i],2),'\n'
		))
	}
	cat('Time taken for ezBuildME() to complete:',round(proc.time()[3]-start),'seconds\n')
	alarm()
	return(to_return)
}

