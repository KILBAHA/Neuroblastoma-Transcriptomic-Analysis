x = c(2:10)
y = c(1456759.08087539, 1412827.20942556, 1359164.51352206, 1365709.1021994, 1276135.79335051, 1276067.53978284, 1233304.68863817, 1237096.79937369, 1215747.61413016)

plot(x,y, type ='p', main = "Elbow Plot", xlab = "K", ylab = "WCSS")

model <- lm(y ~ x + I(x^2) + I(x^3))

# I can get the features of this model :
#summary(model)
#model$coefficients
#summary(model)$adj.r.squared

# For each value of x, I can get the value of y estimated by the model, and add it to the current plot !
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], col=1, lwd=2 )  
