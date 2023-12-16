
col.alpha <- function (acol, alpha = 0.2) 
{
    acol <- col2rgb(acol)
    acol <- rgb(acol[1]/255, acol[2]/255, acol[3]/255, alpha)
    acol
}

draw_circ <- function(x=0,y=0,r=1,pts=200,...) {
    theta = seq(0, 2 * pi, length = pts)
    lines(x = r * cos(theta) + x, y = r * sin(theta) + y,...)
}

pt_in_circ <- function(x,y,r,cx,cy) {
    r^2 > (x-cx)^2 + (y-cy)^2
}

# grid of jittered points in circle of radius r
# first argument is density of points
rptcirc <- function(n=1,r=1,cx=0,cy=0,jitter=0.1,edge_buffer=0.05) {
    i <- 1
    x_range <- seq(-r,r,by=n)+cx
    y_range <- seq(-r,r,by=n)+cy
    m <- length(x_range) * length(y_range)
    out <- data.frame(x=rep(NA,m),y=rep(NA,m))
    for ( x in x_range )
    for ( y in y_range ) {
        # is point in circle?
        xx <- x + rnorm(1,0,jitter)
        yy <- y + rnorm(1,0,jitter)
        do_draw <- pt_in_circ(xx,yy,r-edge_buffer,cx,cy)
        if ( do_draw==TRUE ) {
            # in circle(s), draw it with jitter
            out$x[i] <- xx
            out$y[i] <- yy
            i <- i + 1
        }
    }#xy loop
    return(out[1:(i-1),])
}

#' @export
draw_euler <- function(p,alpha,beta,r3,d,draw_pts=20,jitter=0.1,lwd=3,al=0.6) {
    plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) , xaxt="n" , yaxt="n" , asp=1 , bty="n" , xlab="" , ylab="" )
    pp <- sqrt(p/(1-p))
    set.seed(1)
    draw_circ(1,0,1,col="green",lwd=lwd)
    if ( draw_pts > 0 ) {
        # green
        pts <- rptcirc(draw_pts,1,1,0,jitter=jitter )
        det <- pt_in_circ( pts$x , pts$y , r3 , 1-d , 0 )
        xcols <- ifelse( det==TRUE , "green" , col.alpha("green",al) )
        xcex <- ifelse( det==TRUE , 1.4 , 1 )
        points( pts$x , pts$y , pch=16 , col=xcols , cex=xcex )
    }
    draw_circ(-pp,0,pp,col="red",lwd=lwd)
    if ( draw_pts > 0 ) {
        # red
        pts <- rptcirc(draw_pts,pp,-pp,0, jitter=jitter )
        det <- pt_in_circ( pts$x , pts$y , r3 , 1-d , 0 )
        xcols <- ifelse( det==TRUE , "red" , col.alpha("red",al) )
        xcex <- ifelse( det==TRUE , 1.4 , 1 )
        points( pts$x , pts$y , pch=16 , col=xcols , cex=xcex )
    }
    
    # thicken arcs inside detection circle
    # doesn't work yet - need to focus on included angles I think so it doesn't close the arc with a line segment
    if ( FALSE ) {
        theta = seq(0, 2 * pi, length = 200)
        # green
        x <- 1*cos(theta) + 1
        y <- 1*sin(theta) + 0
        det <- pt_in_circ(x,y,r3,1-d,0)
        lines(x = x[det], y = y[det] , lwd=6 , col="green" )
        # red
        x <- pp*cos(theta) - pp
        y <- pp*sin(theta) + 0
        det <- pt_in_circ(x,y,r3,1-d,0)
        lines(x = x[det], y = y[det] , lwd=6 , col="red" )
    }

    # detection circle on top
    draw_circ(1-d,0,r3,col="black",lwd=lwd*2,lty=1)
}

# draw_euler(0.5,0.5,0.5,1.1587212,1)

#' Solve and draw an Euler screening diagram
#'
#' This function take the base rate and detection probabilities. It then solves
#' and draws the corresponding Euler diagram for the implied screening scenario.
#'
#' @param p Base rate of positive cases
#' @param alpha Probability of true positive, p(T|T)
#' @param beta Probability of false positive, p(T|F)
#' @return Solution parameters
#' @export
oiler_diagram <- function(p,alpha,beta,init=c(sqrt(p/(1-p)),1),draw_pts=0,jitter=0.1) {
    suppressWarnings( theta <- solve_euler_diagram(p,alpha,beta,init=init)$par )
    draw_euler( p , alpha , beta , theta[1] , theta[2] , draw_pts=draw_pts , jitter=jitter )
    print(theta) # radius and offset (from center of green) of black/detection circle
}

# oiler_diagram( p=0.3 , alpha=0.8, beta=0.2 , draw_pts=0.12 , jitter=0.02 )

