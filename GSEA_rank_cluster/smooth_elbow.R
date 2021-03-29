y = read.csv("elbow_para.txt", sep = "\t")
y = y$WCSS
x = c(2:10)
#y = c(955486.888094832,889707.495703698,843749.680404026,815371.073186371,788220.554284853,783514.149646629,763837.507496071,758439.367853814,743315.599335635)

plot(x,y, type ='p', main = "Clustering Performance with \n Varying Cluster Centroid Numbers", xlab = "Number of Cluster Centroids (K)", ylab = "WCSS")

model <- lm(y ~ x + I(x^2) + I(x^3))

# I can get the features of this model :
#summary(model)
#model$coefficients
#summary(model)$adj.r.squared

# For each value of x, I can get the value of y estimated by the model, and add it to the current plot !
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], col=1, lwd=2 )  
