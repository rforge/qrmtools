For former implementation in C
if(is.null(tol)) tol <- -1 # for C code
if(is.infinite(maxiter)) maxiter <- -1 # for C code
rearrange <- NULL # to avoid "RA_aux: no visible binding for global variable 'rearrange_'"
.Call("rearrange_", X, method, tol.type, maxiter, tol)

