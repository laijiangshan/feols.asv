.onAttach <- function(libname, pkgname) {
  cite_info <- utils::citation(pkgname)
  cite_info <- "Jiangshan Lai,Yan He,Mi Hou,Guochun Shen,WenYong Guo,Lingfeng Mao(2026).Quantifying Spatial and Environmental Effects in Spatial Autoregressive Model with the spatialreg.hp R package. Journal of Plant Ecology 19(2):rtaf220"
  packageStartupMessage("Thank you for using this package! If you use this package in your research, please cite the following references \n")
  packageStartupMessage(cite_info)
}
