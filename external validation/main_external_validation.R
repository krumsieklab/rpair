rm( list = ls() )
library(survival)
library(magrittr)
library(rpair)
library(glmnet)

# load the data 
load('external validation/data_for_extarnal_validation.Rdata')

# decite traning set, and test sets
# training set is the one which has highest number of samples 
nodes = lapply(nodes, function(nods){
  # cohorts
  i = names(nodes)
  # sample sizes 
  ns = nods %>% lapply(`[[`,'x') %>% sapply(nrow)
  # training set 
  tr = nods[which.max(ns)]
  # test sets 
  ts = nods[-which.max(ns)]
  # return 
  list(tr=tr, ts=ts)
})

# sample sizes for training set, and test sets
ss = 
sapply(nodes, sapply, sapply, function(x) nrow(x$x)) %>% t %>%
  {data.frame(cancer = rownames(.), .)} %>%
  data.table::as.data.table()
ss$k_genes = sapply(nodes, function(x) ncol(x$tr[[1]]$x) )
ss

# add cancer tag
for(i in names(nodes)) nodes[[i]]$cancer = i


# run the models ----------------------------------------------------------
# traing model with traning set and test on test sets for each cancer
res = structure(names(nodes), names = names(nodes)) %>% lapply(function(i){
  print(i)
  # train the model
  tr = nodes[[i]]$tr[[1]]
  
  # fold ids
  set.seed(42)
  fids = rpair:::get_stratified_folds(tr$y, 5)[,4]
  # fit rpair model and glmnet 
  fit_hin = rpair::cv_rpair(x = tr$x, y = tr$y, loss_type = 'sqh', foldid = fids)
  fit_exp = rpair::cv_rpair(x = tr$x, y = tr$y, loss_type = 'exp', foldid = fids)
  fit_log = rpair::cv_rpair(x = tr$x, y = tr$y, loss_type = 'log', foldid = fids)
  fig = glmnet::cv.glmnet(x = tr$x, y = tr$y, family ='cox', pmax = sum(tr$y[,2]), foldid = fids)
  
  # predict function, if not null use lambda.1se, otherwise lambda.min 
  my_predict <- function(a_fit, newx){
    s = paste0( 'lambda.', if(all(as.vector(coef(a_fit))==0)) 'min' else '1se' )
    as.vector(predict(a_fit, newx, s=s))
  }
  
  # collect results for all test cohorts
  re = 
    lapply(nodes[[i]]$ts, function(ts){
      list(
        `rpair[hin]` = concordance(S~yh, data.frame(S= ts$y, yh = my_predict(fit_hin, ts$x)), reverse = T),
        `rpair[log]` = concordance(S~yh, data.frame(S= ts$y, yh = my_predict(fit_log, ts$x)), reverse = T),
        `rpair[exp]` = concordance(S~yh, data.frame(S= ts$y, yh = my_predict(fit_exp, ts$x)), reverse = T),
        cox = concordance(S~yh, data.frame(S= ts$y, yh = my_predict(fig, ts$x)), reverse = T)
      )
    })
   
  # collect results as data.frame
  df = 
    {lapply(re, sapply, `[`, c('concordance', 'var'))} %>% {
      aa = .;
      lapply(names(aa), function(i) 
        data.frame(cohort = i, model = colnames(aa[[i]]), apply(t(aa[[i]]),2,unlist))) 
    } %>% do.call(what = rbind) %>% {.$se = sqrt(.$var);.} 
  df$cancer = i
  df
})
save(file = 'external validation/results_ext.Rdata', res)  

 # df = do.call(rbind, res)
# ggplot(df, aes(x = model, y = concordance, color = model)) + 
#   geom_pointrange(aes(ymin = concordance - se, ymax = concordance + se)) +
#   geom_point(size = 2) + 
#   facet_grid(cancer~cohort, scales = 'free_y') + 
#   theme_bw() + 
#   scale_color_brewer(type = 'qual', palette = 2, direction = -1)


ggs = 
lapply(res, function(df)
  ggplot(df, aes(x = model, y = concordance, color = model)) + 
    geom_pointrange(aes(ymin = concordance - se, ymax = concordance + se)) +
    geom_point(size = 2) + 
    facet_grid(cancer~cohort, scales = 'free_y') + 
    theme_bw() + 
    scale_color_brewer(type = 'qual', palette = 2, direction = -1)
)
ggs$COAD = ggs$COAD+theme(axis.title.x = element_blank())
ggs$OV = ggs$OV+theme(axis.title.x = element_blank()) + 
  geom_hline(yintercept = 0.5, lty =2, size = 0.5)
ggs$PAAD = ggs$PAAD+theme(axis.title.x = element_blank()) + 
  geom_hline(yintercept = 0.5, lty =2, size = 0.5)

library(patchwork)
patchwork::wrap_plots(ggs, ncol = 1, guides = 'collect')
