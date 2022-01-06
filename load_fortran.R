# get the fortran procedure
ffortran_compile<- function(x, check = TRUE){ 
  re = system(paste0("R CMD SHLIB ",x))==0
  if(check) print(ffotran_checkerror(x))
  re
}

ffotran_checkerror <- function(x){
  system(paste0("R CMD SHLIB ",x," 2>&1 | tee SomeFile.txt"))
  rl = readLines("SomeFile.txt") 
  ei = grep("Error", rl)
  if(length(ei)>0) rl[(ei.-7):(ei+7)] else ""
} 


ffortran_compile("pnet.f")
dyn.load("pnet.so")

print(
  c( plog = is.loaded("plognet"), 
     pexp = is.loaded("pfishnet") 
   )
)

ffortran_compile("pnet2.f90")
dyn.load("pnet2.so")
print(
  c( phuh = is.loaded("phuhnet"),
     psqh = is.loaded("psqhnet") 
   )
)

