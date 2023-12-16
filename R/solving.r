
# area_circ_intersect
ACI <- function(r,R,d) {
    r^2 * acos((d^2+r^2-R^2)/(2*d*r)) + R^2 * acos((d^2+R^2-r^2)/(2*d*R)) - 0.5*sqrt( (-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R) )
}

# test - ACI(1,1,2) should equal zero

# p - base rate of true cases
# alpha - P(true|true)
# beta - P(true|false)
solve_euler_diagram <- function(p,alpha,beta,init=c(1,1)) {
    pp <- sqrt(p/(1-p))
    f <- function(x) {
        r3 <- x[1]
        d <- x[2]
        L <- (ACI(1,r3,d)-beta*pi)^2 + (ACI(pp,r3,1+pp-d)-alpha*pi*p/(1-p))^2
        return(L)
    }
    o <- optim( init , f )
    sols <- o$par
    print( c( ACI(1,sols[1],sols[2]) , beta*pi ) )
    print( c( ACI(pp,sols[1],1+pp-sols[2]) , alpha*pi*p/(1-p) ) )
    return(o)
}

# solve_euler_diagram(0.5,0.5,0.5)

