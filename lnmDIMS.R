library(devtools)
#create("lnmDIMS")
document("~/Documents/lnmDIMS")
document("lnmDIMS")


install("lnmDIMS")
library(lnmDIMS)


#build("lnmDIMS")
library(pkgdown)
build_site("~/Documents/lnmDIMS")
build_site("lnmDIMS")

#use_vignette("~/lnmDIMS")
build_vignettes()

## build reference
## Then I can add function to different catogories

## daft for plotting the diagram