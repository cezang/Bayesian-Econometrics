library(tidyverse)
setwd("C:\\Users\\Cezary\\Desktop\\szko³a\\ekonometria bayesowska\\projekt zaliczeniowy\\college-basketball-dataset")
bb19 <- read.csv("cbb19.csv", header = TRUE, stringsAsFactors = FALSE)
bb18 <- read.csv("cbb18.csv", header = TRUE, stringsAsFactors = FALSE)


bb19 <- bb19[,-c(1,2,22,23)]
sum(is.na(bb19))
bb19 <- bb19[,-20]

model <- lm(W~.,data = bb19)
summary(model)

BIC.lm = step(model, k = log(nrow(bb19)), trace = 0)
summary(BIC.lm)

columns <- append(names(BIC.lm$coefficients)[-1], "W")

bb19.2 <- bb19[,columns]
model <- lm(W~., data= bb19.2)
summary(model)

library(corrplot)
cormat <- cor(bb19.2)
corrplot(cormat,method = "circle", type="full")

"https://statisticalhorizons.com/multicollinearity"
library(car)
vif<- vif(model)
# <10 to regula kciuka
col_vif <- vif[vif<5]
columns_name_vif <- append(names(col_vif), "W")


bb19.3 <- bb19.2[,columns_name_vif]

model <- lm(W~., data=bb19.3)
summary(model)



#Bayesowskie podejscie
library(manipulate)
library(mvtnorm)
#OLS
y <- as.matrix(bb19.3[, c('W')])
N.data <- length(y)
X <- cbind(as.matrix(rep(1, N.data)), 
           as.matrix(bb19.3[, c("G","EFG_O","TOR","TORD","ORB","DRB")]))

Beta.ols.data <- model$coefficients
v.data <- model$df.residual
XTX.data <- t(X) %*% X
s2.data <- sum((model$residuals) ^ 2) / v.data



#informacje a priori z poprzedniego modelu
bb18 <- bb18[,columns_name_vif]
sum(is.na(bb18))

model2 <- lm(W~., data=bb18)
summary(model2)

Beta.prior <- as.vector(model2$coefficients)
sm2.prior <- solve(var(model2$residuals))
U.prior <- vcov(model2)*as.vector(sm2.prior)
v.prior <- length(bb18$W)
vs2.prior <- v.prior / sm2.prior

#informacje a posteriori
Beta.posterior <- solve(solve(U.prior) + XTX.data) %*% (solve(U.prior) %*% Beta.prior + XTX.data %*% Beta.ols.data)
U.posterior <- solve(solve(U.prior) + XTX.data)
v.posterior <- v.prior + N.data
vs2.posterior <- v.prior / sm2.prior + v.data * s2.data + t(Beta.ols.data - Beta.prior) %*% solve(U.prior + solve(XTX.data)) %*% (Beta.ols.data - Beta.prior)
sm2.posterior <- 1 / (vs2.posterior / v.posterior)


#testy normalnosci shapiro wilka
library(nortest)
shapiro.test(model$residuals)
shapiro.test(model2$residuals)



#wykresy


#ZMIENNE PRZEDZIALY BETA SPACE

#Definicje kolorów
blue <- rgb(61, 88, 184, 200, names = NULL, maxColorValue = 255)
orange <- rgb(255, 128, 64, 200, names = NULL, maxColorValue = 255)

#Z mrof lub bez
#par(mfrow = c(3, 3))


#rozklady brzegowe
for(ii in 2:length(Beta.posterior)) {
  beta.space <- seq(from = -1.5, to = 1.5, by = 0.01)
  if (ii==1){
    beta.space <- seq(from = -1.5, to = 1.5, by = 0.01)
  } else if (ii==2){
    beta.space <- seq(from = .4, to = 1.2, by = 0.005)
  } else if (ii==3){
    beta.space <- seq(from = .6, to = 1.2, by = 0.005)
  } else if (ii==4){
    beta.space <- seq(from = -1, to = -0.2, by = 0.005)
  } else if (ii==5){
    beta.space <- seq(from = .4, to = 1.20, by = 0.005)
  } else if (ii==6){
    beta.space <- seq(from = 0, to = 0.8, by = 0.005)
  } else if (ii==7){
    beta.space <- seq(from = -0.8, to = 0, by = 0.005)
  }
  
  n_eval_points <- length(beta.space)
  
  prior.marg.dens.beta <- matrix(NA, nrow = 1, ncol = n_eval_points)
  posterior.marg.dens.beta <- matrix(NA, nrow = 1, ncol = n_eval_points)
  
  prior.marg.dens.beta <- apply(as.matrix(beta.space), 1, dmvt,
                                delta = Beta.prior[ii], 
                                sigma = as.matrix(U.prior[ii, ii] / sm2.prior),
                                df = v.prior, log = FALSE)
  posterior.marg.dens.beta <-apply(as.matrix(beta.space), 1, dmvt,
                                   delta = Beta.posterior[ii],
                                   sigma = as.matrix(U.posterior[ii, ii] / sm2.posterior), 
                                   df = v.posterior, log = FALSE)
  
  plot(beta.space, prior.marg.dens.beta, las = 1, lwd = 2, bty = "n", col = blue,
       ylim = c(0, max(c(max(prior.marg.dens.beta),max(posterior.marg.dens.beta))) + 1), type = "l", ylab = "gêstoœæ",xlab="", main = colnames(bb19.3)[ii-1])

  abline(v = Beta.prior[ii], col = blue, lwd = 3)

  abline(v = Beta.ols.data[ii], col = rgb(0, 0, 0, 1), lwd = 3)
  
  lines(beta.space, posterior.marg.dens.beta, lwd = 2, col = orange)

  abline(v = Beta.posterior[ii], col = orange, lwd = 3)
  
  prior_text<- paste0("E(beta) a priori = ", round(Beta.prior[ii], 4))
  ols_text<- paste0("parametr OLS = ", round(Beta.ols.data[ii], 4))
  post_text <-paste0("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4))
  legend("topright", legend=c(prior_text,ols_text,post_text),
         col=c(blue,rgb(0, 0, 0, 1),orange),lty=c(1,1,1),cex=0.75, ncol=1, lwd=c(2,2,2), text.font = 2)
}


####################
#przedzialy ufnosci hpdi
#par(mfrow = c(3, 3))
for(ii in 2:length(Beta.posterior)) {
  beta.space <- seq(from = -1.5, to = 1.5, by = 0.01)
  if (ii==1){
    beta.space <- seq(from = -1.5, to = 1.5, by = 0.01)
  } else if (ii==2){
    beta.space <- seq(from = .4, to = 1.2, by = 0.005)
  } else if (ii==3){
    beta.space <- seq(from = .6, to = 1.2, by = 0.005)
  } else if (ii==4){
    beta.space <- seq(from = -1, to = -0.2, by = 0.005)
  } else if (ii==5){
    beta.space <- seq(from = .4, to = 1.20, by = 0.005)
  } else if (ii==6){
    beta.space <- seq(from = 0, to = 0.8, by = 0.005)
  } else if (ii==7){
    beta.space <- seq(from = -0.8, to = 0, by = 0.005)
  }
  
  n_eval_points <- length(beta.space)
  
  prior.marg.dens.beta <- matrix(NA, nrow = 1, ncol = n_eval_points)
  posterior.marg.dens.beta <- matrix(NA, nrow = 1, ncol = n_eval_points)
  
  prior.marg.dens.beta <- apply(as.matrix(beta.space), 1, dmvt,
                                delta = Beta.prior[ii], 
                                sigma = as.matrix(U.prior[ii, ii] / sm2.prior),
                                df = v.prior, log = FALSE)
  posterior.marg.dens.beta <-apply(as.matrix(beta.space), 1, dmvt,
                                   delta = Beta.posterior[ii],
                                   sigma = as.matrix(U.posterior[ii, ii] / sm2.posterior), 
                                   df = v.posterior, log = FALSE)
  #poziom ufnosci
  conf_level = 95
  highest <- data.frame(1:n_eval_points, posterior.marg.dens.beta, stringsAsFactors = FALSE)
  highest <- highest[order(highest[, 2], decreasing = TRUE), ]
  highest <- data.frame(highest, cumsum(highest[, 2]) < conf_level, stringsAsFactors = FALSE)
  highest <- highest[order(highest[, 1], decreasing = FALSE), ]
  credible_set_indicator <- as.vector(as.integer(highest[, 3]))
  credible_set_begin <- match(1, credible_set_indicator)
  credible_set_end <- length(credible_set_indicator) - match(1, rev(credible_set_indicator))
  #Lewy i prawy brzeg HPDI
  x1 <- beta.space[credible_set_begin]
  x2 <- beta.space[credible_set_end]
  
  posterior.cs <- posterior.marg.dens.beta * credible_set_indicator
  #Poziom ufnoœci
  HPDI_probab <- conf_level * 0.01
  #Wykres gêstoœci a posteriori
  plot(beta.space, posterior.marg.dens.beta, las = 1, lwd = 2, bty = "n", col = orange,
       ylim = c(0, max(posterior.marg.dens.beta + 1)), type = "l", ylab = "gêstoœæ",xlab="", main = colnames(bb19.3)[ii - 1])
  polygon(c(beta.space, rev(beta.space)), 
          c(posterior.marg.dens.beta, rep(0, length(beta.space))), 
          col = orange, border = NA)
  #text(Beta.posterior[ii], max(posterior.marg.dens.beta) + 0.6, paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = orange)
  abline(v = Beta.posterior[ii], col = orange, lwd = 3)
  #Pole oznaczaj¹ce gêstoœæ a posteriori w przedziale ufnoœci HPDI
  polygon(c(beta.space, rev(beta.space)), 
          c(posterior.cs, rep(0, length(beta.space))), 
          col = blue, border = NA)
  
  #Wyœwietl poziom ufnoœci i granice przedzia³u
  #text(Beta.posterior[ii], max(posterior.marg.dens.beta) + 0.2, paste(round(HPDI_probab * 100, digits = 1), "% przedzia³ HPDI: (", round(x1, digits = 2), " , ", round(x2, digits = 2), ")"), col = blue)
  hpdi_text <- paste0(round(HPDI_probab * 100, digits = 1), "% przedzia³ HPDI: (", round(x1, digits = 2), " , ", round(x2, digits = 2), ")")
  post_text <-paste0("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4))
  legend("topright", legend=c(hpdi_text,post_text),
         col=c(blue,orange),
         density = c(NA,1),
         fill=c(blue,orange),
         lty=c(NA,1),
         lwd =c(NA,2),
         border=c(NA,NA),
         cex=0.7, ncol=1,
         x.intersp=c(0.5,2),text.font = 2)
}




