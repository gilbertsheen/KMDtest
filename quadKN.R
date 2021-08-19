##  Quadratic Kernel
quad_kn <- function(z=NULL, rho=NULL){
	(z %*% t(z) + rho)^2
}
