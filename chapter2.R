# Plots and Examples for Chapter 2
# Original code developed by Mr Jake Kraksa (Monash University)
# Built on R 3.6.3 with R Studio 1.2.5033

### load libraries

library(irtDemo)
library(mirt)

### iccs for different dichomotous parameter models

irtDemo("dich")

iccplot <- function(irtmodel = NULL, df = NULL, title = NULL, colour = 24) {
  par(lab=c(7,3,3))
  theta <- seq(-3, 3, .1)
  a <- df[,1]
  b <- df[,2]
  c <- df[,3]
  d <- df[,4]
  for (i in 1:length(b)) {
    if (irtmodel == "rasch") {
      P <- 1 / (1 + exp(-1 * (theta - b[i])))
    } else if (irtmodel == "twopl") {
      P <- 1 / (1 + exp(-a[i] * (theta - b[i])))
    } else if (irtmodel == "threepl") {
      P <- c[i] + (1 - c[i]) * (1 / (1 + exp(-a[i] * (theta - b[i]))))
    } else if (irtmodel == "fourpl") {
      P <- c[i] + (d[i] - c[i]) * (1 / (1 + exp(-a[i] * (theta - b[i]))))
    }
    plot(theta, 
         P, 
         type = "l", 
         xlim = c(-3,3), 
         ylim = c(0, 1), 
         xlab = "Ability/Theta", 
         ylab = "Probability of Correct Response",
         main = title,
         col = colour[i])
    text(x = b[i], y = (d[i]+c[i])/2, labels = paste0("Item ", i, "\n a = ", a[i], "\n b = ", b[i], "\n c = ", c[i], "\n d = ", d[i]), col = colour[i])
    par(new = TRUE)
  }
  par(mar = c(0, 0, 0, 0))
}

linecols <- c("blue", "red", "black", "green", "orange")

raschitems <- data.frame("a" = 1, "b" = seq(-2,2, length.out = 5), "c" = 0, "d" = 1)
iccplot(irtmodel = "rasch", df = raschitems, colour = linecols)
dev.off()

twoplitems <- data.frame("a" = round(seq(0,2.8, length.out = 5),2), "b" = seq(-2,2, length.out = 5), "c" = 0, "d" = 1)
iccplot(irtmodel = "twopl", df = twoplitems, colour = linecols)
dev.off()

negdiscrimination <- data.frame("a" = round(seq(1, -2.8, length.out = 3),2), "b" = seq(-2,2, length.out = 3), "c" = 0, "d" = 1)
iccplot(irtmodel = "twopl", df = negdiscrimination, colour = linecols)
dev.off()

threeplitems <- data.frame("a" = 1, "b" = seq(-2,2, length.out = 5), "c" = round(seq(0,.35, length.out = 5), 2), "d" = 1)
iccplot(irtmodel = "threepl", df = threeplitems, colour = linecols)
dev.off()

fourplitems <- data.frame("a" = 1, "b" = seq(-2, 1, length.out = 3), "c" = 0, "d" = seq(0.5, 1,length.out = 3))
iccplot(irtmodel = "fourpl", df = fourplitems, colour = linecols)
dev.off()

### polytomous category response curves

irtDemo("grm")

### information curves

dat <- expand.table(LSAT6)
mod <- mirt(dat, 1)
M2(mod)
itemfit(mod)
one <- itemplot(mod, 1, 'info')
two <- itemplot(mod, 2, 'info')
four <- itemplot(mod, 4, 'info')
five <- itemplot(mod, 5, 'info')

print(one + two + four + five)

plot(mod, type = 'info', facet_items=FALSE, shiny = TRUE)
