rm(list = ls())
library(survival)
library(Biobase)
library(SummarizedExperiment)
library(magrittr)

# ovarian cancer cohorts --------------------------------------------------
library(MetaGxOvarian)
esets0 <- MetaGxOvarian::loadOvarianEsets()[[1]]

# decide which data needs to be used 
l = sapply(esets0, function(aa) 
  min( na.omit(aa$vital_status) %>% length, na.omit(aa$days_to_death) %>% length ) )

# keep 4 cohorts which have most sample sizes
l = names(which(l>100)) #
genes = lapply( esets0[l], function(x) rownames( Biobase::assayData(x)$exprs ) )
genes = unlist(genes) %>% table  %>% {names(which(.==max(.)))}
l = l[sapply( esets0[l], function(x) 
  all( genes %in% rownames( Biobase::assayData(x)$exprs ) ))]

# datasets to consider
esets = esets0[l]
esets=
  lapply( esets, function(x){
    j = rownames(Biobase::assayData(x)$exprs) %in% genes
    i = (!is.na(x$days_to_death)) & (!is.na(x$vital_status))
    x = x[j,i]
    print(dim(x))
    x
  })

# prepare datesets for each cohort
nods_ov = lapply(esets, function(a){
  x = scale(t(Biobase::assayData(a)$exprs))
  y = Surv(a$days_to_death+1, a$vital_status=='deceased')
  list(x=x, y=y)
})

# pancreas cancer  --------------------------------------------------------
library(MetaGxPancreas)
esets0 <- MetaGxPancreas::loadPancreasDatasets()[[1]]

# decide which data needs to be used 
l = sapply(esets0, function(aa) min( na.omit(aa$vital_status) %>% length,
                                     na.omit(aa$days_to_death) %>% length ) )
# keep 4 cohorts which have most sample sizes
# l = names(which(l>100)) #
l = names(l[order(l,decreasing = T)[1:4]])

genes = lapply( esets0[l], function(x) rownames( assay(x) ) )
genes = unlist(genes) %>% table  %>% {names(which(.==max(.)))}
l = l[sapply( esets0[l], function(x) all( genes %in% rownames(assay(x) ) ))]

# datasets to consider
esets = esets0[l]

esets=
  lapply( esets, function(x){
    j = rownames(assay(x)) %in% genes
    x$days_to_death = as.numeric(as.character(x$days_to_death)) 
    i = (!is.na(x$days_to_death)) & (!is.na(x$vital_status))
    x = x[j,i]
    print(dim(x))
    x
  })

# prepare data for analysis 
nods = lapply(esets, function(a){
  x = scale(t(assay(a)))
  y = Surv(a$days_to_death+1, a$vital_status=='deceased')
  list(x=x, y=y)
})

# discard constant variance features
inds = rowSums( sapply(nods, function(x) colSums(is.na(x$x))) ) == 0
nods_paad = lapply(nods, function(x){
  x$x = x$x[,inds]
  x
})


# colorectal adenocarcinoma  ----------------------------------------------
library(mcsurvdata)

# get COAD data 
eh <- ExperimentHub()
dat <- query(eh, "mcsurvdata")
nda.coad <- dat[["EH1498"]]

# survival outcome 
S =  Surv(nda.coad$tev+1, nda.coad$evn == 1)
nda.coad = nda.coad[,!is.na(S)] 

# gene expression data
mrna = Biobase::assayData(nda.coad)$exprs %>% t

# survival outcome 
S =  data.frame(t = nda.coad$tev+1, e = nda.coad$evn == 1)

# cohorts 
Z = paste0(nda.coad$dataset, nda.coad$cohort)

# split data to cohorts 
nods_coad <-
  lapply( structure( unique(Z), names = unique(Z) ), function(i){
    print(i)
    inds = (Z==i)
    x = as.matrix(scale(mrna[inds,]))
    y = Surv( S[inds,"t"], S[inds,"e"] )
    list(x=x, y=y)
  })

rm(list = setdiff(ls(), c('nods_coad', 'nods_ov', 'nods_paad')))
names(nods_paad) = gsub(pattern = '_SumExp', replacement = '', names(nods_paad) )

# collect datasets 
nodes = list(COAD = nods_coad, OV = nods_ov, PAAD = nods_paad)

# save datasets for extarnal evaluation
save(file = 'external validation/data_for_extarnal_validation.Rdata', nodes)


