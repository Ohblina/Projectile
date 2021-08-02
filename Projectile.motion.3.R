###Installing libraries
install.packages("dlm")
library(dlm)

#if we assume that the projectile's motion is influenced
#by factors such as drag and wind gusts, then the equations of
#motions are (in Cartesian coordinates):

#for the x-direction:
#position = x0 + v0*cos(alpha)*t+noise
#velocity = v0*cos(alpha)+noise
#acceleration = 0

#for the y-direction:
#position = y0 + v0*sin(alpha)*t - g*t+noise
#velocity = v0*sin(alpha) - .5*g*t^2+noise
#acceleration = -g

#where x0 is the x-position at t=0, y0 the y-position at t=0,
#v0 the initial (launching) velocity, alpha the angle with the x-axis,
#and g the acceleration in the y direction (which is 9.81)


# we assume that we measure at each time step
#- the position in the x and y direction
#- the velocity in the x and y direction
set.seed(1)

##TRUE STATE

#x-position at t=0
x0 <- 0

#y-position at t=0
y0 <- 0

#velocity at t=0
v0 <- 16

#acceleration in y-direction
g <- 9.81

#launching angle (in degrees)
alphaDeg <- 75
#convert degrees to radians (by multiplying with pi/180)
alpha <- alphaDeg*(pi/180)

#flight time
2*v0/g*sin(alpha)

#observation times
t <- seq(0.1, 3, .1)

#compute x-position at the observation times
x <- x0 + v0*cos(alpha)*t+rnorm(length(t), 0, 0.1)

#compute y-position at the observation times
y <- y0 + v0*sin(alpha)*t - .5*g*t^2+rnorm(length(t), 0, 0.1)

#trajectory of the projectile
plot(x, y, type="l", xlab="x-position", ylab="y-position", main = "True state")

#compute velocity in the x-direction at the observation times
vx <- v0*cos(alpha)+rnorm(length(t), 0, 0.1)

#compute velocity in the y-direction at the observation times
vy <- v0*sin(alpha) - g*t+rnorm(length(t), 0, 0.1)


##Measurements

#assume that the measurement noise is normally distributed
#with mean=0 and sd=0.5

#x-position
xm <- x0 + v0*cos(alpha)*t + rnorm(length(t), 0, 0.5)

#y-position
ym <- y0 + v0*sin(alpha)*t - .5*g*t^2 + rnorm(length(t), 0, 0.5)

#velocity in x-direction
vxm <- rep(v0*cos(alpha), length(t)) + rnorm(length(t), 0, 0.5)

#velocity in y-direction
vym <- v0*sin(alpha) - g*t + rnorm(length(t), 0, 0.5)

#measured trajectory of the projectile
plot(x, y, type="l", xlab="x-position", ylab="y-position", ylim=c(-2, 14), main = "True state vs Measurements")
lines(xm, ym, type="o", lty=2, col="red")
legend("topleft", pch=c(NA, 1), lty=c(1, 2),
       col=c("black", "red"), legend=c("True State", "Measurements"),
       bty="n", y.intersp=1.2)

#store measurement data
dataEx1 <- cbind(xm, vxm, ym, vym)

##KALMAN FILTER
dt <- 0.1 #time step
#dynamic linear model
#specifying 5 states, namely [x, vx, y, vy, uy]
#where x is the x-position, vx the velocity in the x direction,
#y the y-position, vy the velocity in the y direction, and uy
#the acceleration in the y direction
ex1 <- dlm(m0=c(0, v0*cos(alpha), 0, v0*sin(alpha), -g), #initial state estimates
           #error covariances of the initial state estimates:
           #this vector reflects the uncertainty in our initial state estimates
           #you may change the values in this vector and see how they influence
           #the behavior of the Kalman filter
           C0=diag(c(rep(0.1,5))),
           #observation matrix
           FF=matrix(c(1,0,0,0,0,
                       0,1,0,0,0,
                       0,0,1,0,0,
                       0,0,0,1,0), nrow=4, byrow=TRUE),
           #measurement noise
           V=diag(c(rep(0.5,4))),
           #state transition matrix
           GG=matrix(c(1,dt,0,0,0,
                       0,1,0,0,0,
                       0,0,1,dt,.5*dt^2,
                       0,0,0,1,dt,
                       0,0,0,0,1), nrow=5, byrow=TRUE),
           #process noise
           W=diag(rep(0.1,5)))

#note that all covariances in the process noise matrix (W) are set not to zero
#this makes sense since, for instance, the change in x-position at each time
#step is not fully explained by our physical model describing the projectile's motion
#furthermore, we assumed that there were factors (such as random wind gusts)
#that could possibly influence the x-position at each time step in some
#unpredictable (random) way

##compute the filtered (a posteriori) state estimates
filteredState <- dlmFilter(dataEx1, ex1)
fullset<-cbind(x,y,xm,ym,filteredState$m[2:31,1],filteredState$m[2:31,3])
colnames(fullset)<-c("x","y","xm","ym","x.EKF","y.EKF")
fullset<-as.data.frame(fullset)
xy<-cbind(x,y)

#plot the filtered state estimates 
plot(x, y, type="o", xlab="x-position", ylab="y-position",
      col=gray(level=.7),main = "True State, Measurements, EKF State")
lines(xm, ym, type="o", lty=2, col="red", cex=.7)
lines(filteredState$m[,1], filteredState$m[,3],type="o",lty=2, col="blue", cex=.7)

legend("topleft", pch=c(1, 1, 1), lty=c(1, 2, 2), lwd=c(1, 1, 1),
       col=c(gray(level=.7), "red", "blue"),
       legend=c("True State", "Measurements", "EKF State"),
       bty="n", y.intersp=1.2, cex=.7)

EKF.xy<-cbind(filteredState$m[2:31,1], filteredState$m[2:31,3])
Measure.xy<-cbind(xm,ym)

##Finding Mean Squer Error and verify that EKF has lower error than measurements
MSE_MS<-mean((xy-Measure.xy)^2)
MSE_EKF<-mean((xy-EKF.xy)^2)

MSE_MS
MSE_EKF


#### Exercise EKF
#What if Vo = 21, launching angle = 45
#Q=Process noise=0.05, 
#R=Measurement noise=1
#What are MSE for Measurements and for EKF?
set.seed(1)

##TRUE STATE

#x-position at t=0
x0 <- 0

#y-position at t=0
y0 <- 0

#velocity at t=0
v0 <- 21

#acceleration in y-direction
g <- 9.81

#launching angle (in degrees)
alphaDeg <- 45
#convert degrees to radians (by multiplying with pi/180)
alpha <- alphaDeg*(pi/180)

#flight time
2*v0/g*sin(alpha)

#observation times
t <- seq(0.1, 3, .1)

#compute x-position at the observation times
x <- x0 + v0*cos(alpha)*t+rnorm(length(t), 0, 0.1)

#compute y-position at the observation times
y <- y0 + v0*sin(alpha)*t - .5*g*t^2+rnorm(length(t), 0, 0.1)

#trajectory of the projectile
plot(x, y, type="l", xlab="x-position", ylab="y-position", main = "True state")

#compute velocity in the x-direction at the observation times
vx <- v0*cos(alpha)+rnorm(length(t), 0, 0.1)

#compute velocity in the y-direction at the observation times
vy <- v0*sin(alpha) - g*t+rnorm(length(t), 0, 0.1)


##Measurements

#assume that the measurement noise is normally distributed
#with mean=0 and sd=0.5

#x-position
xm <- x0 + v0*cos(alpha)*t + rnorm(length(t), 0, 0.5)

#y-position
ym <- y0 + v0*sin(alpha)*t - .5*g*t^2 + rnorm(length(t), 0, 0.5)

#velocity in x-direction
vxm <- rep(v0*cos(alpha), length(t)) + rnorm(length(t), 0, 0.5)

#velocity in y-direction
vym <- v0*sin(alpha) - g*t + rnorm(length(t), 0, 0.5)

#measured trajectory of the projectile
plot(x, y, type="l", xlab="x-position", ylab="y-position", ylim=c(-2, 14), main = "True state vs Measurements")
lines(xm, ym, type="o", lty=2, col="red")
legend("topleft", pch=c(NA, 1), lty=c(1, 2),
       col=c("black", "red"), legend=c("True State", "Measurements"),
       bty="n", y.intersp=1.2)

#store measurement data
dataEx1 <- cbind(xm, vxm, ym, vym)

##KALMAN FILTER
dt <- 0.1 #time step
#dynamic linear model
#specifying 5 states, namely [x, vx, y, vy, uy]
#where x is the x-position, vx the velocity in the x direction,
#y the y-position, vy the velocity in the y direction, and uy
#the acceleration in the y direction
ex1 <- dlm(m0=c(0, v0*cos(alpha), 0, v0*sin(alpha), -g), #initial state estimates
           #error covariances of the initial state estimates:
           #this vector reflects the uncertainty in our initial state estimates
           #you may change the values in this vector and see how they influence
           #the behavior of the Kalman filter
           C0=diag(c(rep(0.1,5))),
           #observation matrix
           FF=matrix(c(1,0,0,0,0,
                       0,1,0,0,0,
                       0,0,1,0,0,
                       0,0,0,1,0), nrow=4, byrow=TRUE),
           #measurement noise=R
           V=diag(c(rep(1,4))),
           #state transition matrix
           GG=matrix(c(1,dt,0,0,0,
                       0,1,0,0,0,
                       0,0,1,dt,.5*dt^2,
                       0,0,0,1,dt,
                       0,0,0,0,1), nrow=5, byrow=TRUE),
           #process noise=Q
           W=diag(rep(0.01,5)))

#note that all covariances in the process noise matrix (W) are set not to zero
#this makes sense since, for instance, the change in x-position at each time
#step is not fully explained by our physical model describing the projectile's motion
#furthermore, we assumed that there were factors (such as random wind gusts)
#that could possibly influence the x-position at each time step in some
#unpredictable (random) way

##compute the filtered (a posteriori) state estimates
filteredState <- dlmFilter(dataEx1, ex1)
fullset<-cbind(x,y,xm,ym,filteredState$m[2:31,1],filteredState$m[2:31,3])
colnames(fullset)<-c("x","y","xm","ym","x.EKF","y.EKF")
fullset<-as.data.frame(fullset)
xy<-cbind(x,y)

#plot the filtered state estimates 
plot(x, y, type="o", xlab="x-position", ylab="y-position",
     col=gray(level=.7),main = "True State, Measurements, EKF State")
lines(xm, ym, type="o", lty=2, col="red", cex=.7)
lines(filteredState$m[,1], filteredState$m[,3],type="o",lty=2, col="blue", cex=.7)

legend("topleft", pch=c(1, 1, 1), lty=c(1, 2, 2), lwd=c(1, 1, 1),
       col=c(gray(level=.7), "red", "blue"),
       legend=c("True State", "Measurements", "EKF State"),
       bty="n", y.intersp=1.2, cex=.7)

EKF.xy<-cbind(filteredState$m[2:31,1], filteredState$m[2:31,3])
Measure.xy<-cbind(xm,ym)

##Finding Mean Squer Error and verify that EKF has lower error than measurements
MSE_MS<-mean((xy-Measure.xy)^2)
MSE_EKF<-mean((xy-EKF.xy)^2)

MSE_MS
MSE_EKF

