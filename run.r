# setup
# clear the environment
# rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
    if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
    dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
    result <- paste(x1,x2,sep="")
    return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
    if (x %in% rownames(installed.packages())) { 
        print(concat("package already installed: ", x))
    }
    else { 
        install.packages(x, repos="http://mirror.las.iastate.edu/CRAN/") 
    }
    library(x, character.only=TRUE)
}

load_package("MVA")
load_package("car")
load_package("psych")
load_package("ggplot2")
load_package("data.table")
load_package("mvtnorm")

#######################################################
# measure.csv                                         #
#######################################################
# chest
# waist
# hips
# gender

# import the file
data <- read.csv(concat(DATA_DIR,'/measure.csv'), header=T)
# describe the data
str(data)
head(data)
write.table(describe(data), file=concat(OUTPUT_DIR,'/measure - descriptions.csv'), sep=",")
# create variable vectors
chest  <- data$chest
waist  <- data$waist
hips   <- data$hips
gender <- data$gender

# function to plot histograms
plot_histogram <- function (d, v, name) {
    # v is the vector we are plotting
    # name is the name of the vector
    y_max_offset <- function(y, offset) {
        y - (0.1 * y * offset)
    }
    m <- mean(v)
    y_max <- max(v)
    ggplot(d, aes(x=v), environment = environment()) +
        geom_histogram(binwidth=1, col="black", aes(fill=..count..)) +
        geom_vline(data=d, aes(xintercept=m),
                   linetype="dashed", size=1, colour="black") +
        scale_fill_gradient("Count", low = "green", high = "red") +
        annotate("text", x = m+0.5*m, y = y_max, label = concat("mean: ", round(m,4)), hjust = 0) +
        annotate("text", x = m+0.5*m, y = y_max_offset(y_max, 1), label = concat("skew: ", round(skew(v),4)), hjust = 0) +
        annotate("text", x = m+0.5*m, y = y_max_offset(y_max, 2), label = concat("kurtosis: ", round(kurtosi(v),4)), hjust = 0) +
        labs(title=concat("Histogram for ", name)) +
        labs(x=name, y="Count")
}

# histogram of chest
plot_histogram(data, chest, "chest")

# covariance matrix - have to exclude gender b/c it is a factor variable
sigma <- cov(data[,1:3])
cov(data[gender == "male",1:3])
cov(data[gender == "female",1:3])

# correlation matrix
cor(data[,1:3])
cor(data[gender == "male",1:3])
cor(data[gender == "female",1:3])

# test for normality
x <- data[, 1:3]
cm <- colMeans(x)
S <- cov(x)
d <- apply(x, MARGIN=1, function(x)
    t(x - cm) %*% solve(S) %*% (x - cm))
qqnorm(data[,"chest"], main="chest"); qqline(data[,"chest"])
qqnorm(data[,"waist"], main="waist"); qqline(data[,"waist"])
qqnorm(data[,"hips"], main="hips"); qqline(data[,"hips"])
plot(qc <- qchisq((1:nrow(x) - 1/2) / nrow(x), df=3), sd <- sort(d),,
     xlab=expression(paste(chi[3]^2, " Quantile")),
     ylab="Ordered distances")
oups <- which(rank(abs(qc - sd), ties="random") > nrow(x) - 3)
text(qc[oups], sd[oups] - 1.5, names(oups))
abline(a=0,b=1)


###########################
# US Air pollution data    #
###########################
# test for normality with qqplots and chi-squared plot
data("USairpollution")
layout(matrix(1:8, nc=2))
sapply(colnames(USairpollution), function(x) {
    qqnorm(USairpollution[[x]], main=x)
    qqline(USairpollution[[x]])
})
layout(matrix(1:1,nc=1))
x <- USairpollution
cm <- colMeans(x)
S <- cov(x)
d <- apply(x, MARGIN=1, function(x)
    t(x - cm) %*% solve(S) %*% (x - cm))
p <- plot(qc <- qchisq((1:nrow(x) - 1/2) / nrow(x), df=6), 
          sd <- sort(d),
          xlab=expression(paste(chi[3]^2, " Quantile")),
          ylab="Ordered distances",
          xlim=range(qc) * c(1, 1.1))
oups <- which(rank(abs(qc - sd), ties="random") > nrow(x) - 3)
text(qc[oups], sd[oups] - 1.5, names(oups))
abline(a=0,b=1)

#############################
# cars data from 1979       #
#############################
data("mtcars")
str(mtcars)
layout(matrix(1:10, nc=2))
d <- mtcars[,2:10]
sapply(colnames(d), function(x) {
    qqnorm(d[[x]], main=x)
    qqline(d[[x]])
})
layout(matrix(1:1,nc=1))
x <- d
cm <- colMeans(x)
S <- cov(x)
d <- apply(x, MARGIN=1, function(x)
    t(x - cm) %*% solve(S) %*% (x - cm))
p <- plot(qc <- qchisq((1:nrow(x) - 1/2) / nrow(x), df=9), 
          sd <- sort(d),
          xlab=expression(paste(chi[3]^2, " Quantile")),
          ylab="Ordered distances",
          xlim=range(qc) * c(1, 1.1))
oups <- which(rank(abs(qc - sd), ties="random") > nrow(x) - 3)
text(qc[oups], sd[oups] - 1.5, names(oups))
abline(a=0,b=1)

############################
# student risk-taking      #
############################
data("students", package = "HSAUR2")
layout(matrix(1:2, ncol = 2))
boxplot(low ~ treatment, data = students, ylab = "low")
boxplot(high ~ treatment, data = students, ylab = "high")
str(students)

#################################
# pottery data                  #
#################################
data("pottery", package = "HSAUR2")
scatterplotMatrix(pottery)
str(pottery)
df <- pottery[,1:9]
png(concat(IMAGES_DIR,'/pottery - qqplots.png'), 
    width = 2048, height = 2048)
layout(matrix(1:10, nc=2))
sapply(colnames(df), function(x) {
    qqnorm(df[[x]], main=x)
    qqline(df[[x]])
})
dev.off()
layout(matrix(1:1,nc=1))
x <- df
cm <- colMeans(x)
S <- cov(x)
d <- apply(x, MARGIN=1, function(x)
    t(x - cm) %*% solve(S) %*% (x - cm))
p <- plot(qc <- qchisq((1:nrow(x) - 1/2) / nrow(x), df=8), 
          sd <- sort(d),
          xlab=expression(paste(chi[3]^2, " Quantile")),
          ylab="Ordered distances",
          xlim=range(qc) * c(1, 1.1))
oups <- which(rank(abs(qc - sd), ties="random") > nrow(x) - 3)
text(qc[oups], sd[oups] - 1.5, names(oups))
abline(a=0,b=1)


####################################
####################################
##                                ##
## CHAPTER 2                      ##
##                                ##
####################################
####################################
# setup
# clear the environment
# rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
    if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
    dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
    result <- paste(x1,x2,sep="")
    return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
    if (x %in% rownames(installed.packages())) { 
        print(concat("package already installed: ", x))
    }
    else { 
        install.packages(x, repos="http://mirror.las.iastate.edu/CRAN/") 
    }
    library(x, character.only=TRUE)
}

load_package("MVA")
load_package("car")
load_package("psych")
load_package("ggplot2")
data("USairpollution")
mlab <- "Manufacturing enterprises with 20 or more workers"
plab <- "Population size (1970 census) in thousands"
# scatter plot of manu and popul
plot(popul ~ manu, data=USairpollution, xlab=mlab, ylab=plab)
# divide the device into three plotting areas
rug(USairpollution$manu, side=1)
rug(USairpollution$popul, side=2)
layout(matrix(c(2,0,1,3), nrow=2, byrow=TRUE),
       widths=c(2,1), heights=c(1,2), respect=TRUE)
xlim <- with(USairpollution, range(manu)) * 1.1
plot(popul ~ manu, data=USairpollution, cex.lab=0.9,
     xlab=mlab,ylab=plab, type="n",xlim=xlim)
with(USairpollution, text(manu, popul, cex=0.6,labels=abbreviate(row.names(USairpollution))))
with(USairpollution, hist(manu, main="",xlim=xlim))
with(USairpollution, boxplot(popul))
# show outliers
outcity <- match(lab <- c("Chicago", "Detroit", "Cleveland", "Philadelphia"),
                 rownames(USairpollution))
x <- USairpollution[, c("manu", "popul")]
bvbox(x, mtitle="",xlab=mlab, ylab=plab)
text(x$manu[outcity], x$popul[outcity], labels=lab,
     cex=0.7, pos=c(2,2,4,2,2))
# find the correlations between manu and popul
with(USairpollution, cor(manu, popul))
outcity <- match(c("Chicago", "Detroit", "Cleveland", "Philadelphia"),
                 rownames(USairpollution))
with(USairpollution, cor(manu[-outcity],popul[-outcity]))
# use the convex hull to eliminate outliers
(hull <- with(USairpollution, chull(manu,popul)))
with(USairpollution,
     plot(manu,popul, pch=1,xlab=mlab,ylab=plab))
with(USairpollution,
     polygon(manu[hull], popul[hull], density=15, angle=30))
with(USairpollution, cor(manu[-hull], popul[-hull]))
# analyize the independence using a chiplot
layout(matrix(c(2,2,1,1), nrow=2, byrow=TRUE),
       widths=c(2,2), heights=c(2,2), respect=TRUE)
layout.show(2)
with(USairpollution,
     plot(manu,popul, pch=1,xlab=mlab,ylab=plab),cex.lab=0.9)
with(USairpollution, chiplot(manu,popul))
# bubble plot for 3 variables
ylim <- with(USairpollution, range(wind)) * c(0.95,1)
layout(matrix(1))
layout.show(1)
plot(wind~temp, data=USairpollution,
     xlab="Average Annual Temperature (Fahrenheit)",
     ylab="Average Annual Wind Speed (mph)", pch=10,
     ylim=ylim)
with(USairpollution, symbols(temp, wind, circles=SO2,
                             inches=0.5, add=TRUE))
# represent all 7 variables with a star plot
stars(USairpollution, cex=0.55)
plot(USairpollution)
scatterplotMatrix(USairpollution)
round(cor(USairpollution), 4)

# 3-D scatterplots
epa <- function(x,y) ((x^2 + y^2) < 1) * 2/pi * (1-x^2-y^2)
x <- seq(from=-1.1, to=1.1, by=0.05)
epavals <- sapply(x, function(a) epa(a,x))
persp(x=x,y=x,z=epavals,xlab="x",ylab="y",
      zlab=expression(K(x,y)), theta=-35, axes=TRUE,
      box=TRUE)
library("KernSmooth")
CYGOB1d <- bkde2D(CYGOB1, bandwidth=sapply(CYGOB1, dpik))
plot(CYGOB1, xlab="log surface temperature",
     ylab="log light intensity")
contour(x=CYGOB1d$x1, y=CYGOB1d$x2,
        z=CYGOB1d$fhat, add=TRUE)
persp(x=CYGOB1d$x1, y=CYGOB1d$x2, z=CYGOB1d$fhat,
      xlab="log surface tempurature",
      ylab="log light intensity",
      zlab="density")
load_package("scatterplot3d")
data <- read.csv(concat(DATA_DIR,'/measure.csv'), header=T)
with(data, scatterplot3d(chest,waist,hips,
                            pch=(1:2)[gender],type="h",angle=55))

with(USairpollution, scatterplot3d(temp, wind, SO2, type="h", angle=55))
stalac(USairpollution)

####################################
####################################
##                                ##
## CHAPTER 3 - PCA                ##
##                                ##
####################################
####################################
# headsize data
library("boot")
data("frets")
str(frets)
summary(frets)
head(frets)
head_dat <- frets[,c("l1","l2")]
colMeans(head_dat)
cov(head_dat)
cor(head_dat)
(head_pca <- princomp(x=head_dat))
summary(head_pca, loadings=TRUE, scale=TRUE)
describe(head_dat)
a1 <- 183.84-0.721*185.72/0.693
b1 <- 0.721/0.693
a2 <- 183.84-(-0.693*185.72/0.721)
b2 <- -0.693/0.721
plot(head_dat, xlab="First son's head length (mm)",
     ylab="Second son's head length (mm)")
abline(a1, b1)
abline(a2,b2, lty=2)
xlim <- range(head_pca$scores[,1])
plot(head_pca$scores, xlim=xlim, ylim=xlim)
abline(0,0)

# heptathlon data from 1988
data("heptathlon")
# score all seven events in the same direction so that large values indicate better performance
heptathlon$hurdles <- with(heptathlon, max(hurdles)-hurdles)
heptathlon$run200m <- with(heptathlon, max(run200m)-run200m)
heptathlon$run800m <- with(heptathlon, max(run800m)-run800m)
score <- which(colnames(heptathlon)=="score")
round(cor(heptathlon[,-score]), 2)
scatterplotMatrix(heptathlon[-score])
# papua new guinea participant throws off correlations
heptathlon <- heptathlon[-grep("PNG", rownames(heptathlon)),]
score <- which(colnames(heptathlon) == "score")
round(cor(heptathlon[,-score]),2)
scatterplotMatrix(heptathlon[-score])
heptathlon_pca <- prcomp(heptathlon[,-score], scale=TRUE)
heptathlon_pca
summary(heptathlon_pca)
a1 <- heptathlon_pca$rotation[,1]
a1
center <- heptathlon_pca$center
scale <- heptathlon_pca$scale
hm <- as.matrix(heptathlon[,-score])
drop(scale(hm, center=center, scale=scale) %*%
         heptathlon_pca$rotation[,1])
predict(heptathlon_pca)[,1]
plot(heptathlon_pca)
cor(heptathlon$score, heptathlon_pca$x[,1])
plot(heptathlon$score, heptathlon_pca$x[,1])
PCbiplot <- function(PC, rownames, x="PC1", y="PC2") {
    # code is modified but mostly borrowed from:
    #    http://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
    #    posted by http://stackoverflow.com/users/577462/crayola
    # PC being a prcomp object
    data <- data.frame(obsnames=rownames, PC$x)
    plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
    plot <- plot + geom_hline(alpha=0.4, size=.2, yintercept=0) + geom_vline(alpha=0.4, size=.2, xintercept=0)
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
    )
    datapc <- transform(datapc,
                        v1 = .7 * mult * (get(x)),
                        v2 = .7 * mult * (get(y))
    )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
    plot
}
PCbiplot(heptathlon_pca, row.names(heptathlon))
