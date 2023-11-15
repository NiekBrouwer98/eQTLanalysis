x = rlnorm(500,0,0.01)
grid = seq(0.5,1.5,0.01)

plot(dlnorm(grid,0,0.01),type="l",xlab="x",ylab="f(x)", col=1)
par(new=TRUE)
plot(dlnorm(grid,0,0.1),type="l",xlab="x",ylab="f(x)", col=2)
par(new=TRUE)
plot(dlnorm(grid,0,0.25),type="l",xlab="x",ylab="f(x)", col=3)
par(new=TRUE)
plot(dlnorm(grid,0,0.5),type="l",xlab="x",ylab="f(x)", col=4)
par(new=TRUE)
plot(dlnorm(grid,0,1),type="l",xlab="x",ylab="f(x)", col=5)
# lines(density(x),col="red")


legend("topright",c("std = 0.01", "std = 0.1","std = 0.25", "std = 0.5", "std = 1"),lty=1,col=1:5, cex= 0.5)
