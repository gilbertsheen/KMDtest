## Gaussian Kernel
gau_kn <- function(z, rho){
	t.rho <- (apply(z,2,sd))^2
	a.rho <- sum(t.rho) * ta
	exp(-(as.matrix((dist(z))^2)/a.rho))
}