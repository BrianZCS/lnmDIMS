library(devtools)
#create("lnmDIMS")
document("~/Documents/lnmDIMS")

install("lnmDIMS")
library(lnmDIMS)


#build("lnmDIMS")
library(pkgdown)
build_site("~/Documents/lnmDIMS")
#use_vignette("~/lnmDIMS")
#build_vignettes()
