# Check argument before each do.call or exec run
CheckArgument = function(fun, arglist){
    avaiable_args    = fun %>% args %>% as.list() %>% names
    not_support_args = setdiff(names(arglist), avaiable_args)
    fun_name = substitute(fun) %>% as.character
    if(length(not_support_args) != 0) stop(str_glue('Found not support argument: "{toString(not_support_args)}" for function "{fun_name}"'))
}