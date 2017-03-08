###Part 1
### Get data
data_o<-read.csv("HUDM6026.csv")
data<-data_o[,!names(data_o) %in% c("X...CHILDID","filter_.")]

#explore the structure of data frame
str(data)
head(data)

#build full model with all variables
ful_mod <- lm(C7R2SSCL ~ ., data = data)
summary(ful_mod)

# build null model containing only response variable
null_mod <- lm(C7R2SSCL ~ 1, data = data)

## use Stepwise Selection to select best fit model

library(MASS)
### AIC model selection
bestAIC <- stepAIC(ful_mod,direction = "backward") #start with full model, backward stepwise selection
summary(bestAIC)


bestAIC2 <- stepAIC(null_mod, 
                    scope = list(upper=ful_mod, lower=null_mod),
                    direction = "forward") # start with null model, forward stepwise selection
summary(bestAIC2)

bestAIC3 <- stepAIC(null_mod,scope = list(upper=ful_mod, lower=null_mod),
                    direction = "both") # from both direction, the result is the 
summary(bestAIC3)                       # forward or backward conditioning on initial model 

# BIC model selection
n <- nrow(data)
bestBIC  <- stepAIC(ful_mod, direction = "both",  k = log(n))
summary(bestBIC)

###########################################
### Part 2

#####################################################################
################    Ols Model    ####################################
#####################################################################

### Fit ols by glm
library(boot)
glm1 <- glm(C7R2SSCL ~ ., data = data, family = gaussian)  # Full model
glm2 <- glm(C7R2SSCL ~ 1, data = data, family = gaussian)  # Null model
glm3 <- glm(C7R2SSCL~GENDER+WKWHITE+WKSESL+W8INCCAT+C7ANGRY+C7LIKRD+C7LONLY+C7FLGOOD+S7REGSKL 
+C7R4RSCL+C7R4MSCL+S2AFTSCH+S2CMNITY+S2TRFFIC+S2TNSION,data = data, family = gaussian) # Stepwise selection 
summary(glm1)
### The default cost function is mean squared error, which works for continuous outcome
set.seed(1)
cv.err <- cv.glm(data = data, glmfit = glm1, K = 10)
cv.err1 <- cv.glm(data = data, glmfit = glm2, K = 10)
cv.err2 <- cv.glm(data = data, glmfit = glm3, K = 10)

### Get 10-fold CV error estimate
cv.err$delta[1] 
cv.err1$delta[1] 
cv.err2$delta[1] 



#####################################################################
################ Lasso and Ridge ####################################
#####################################################################

### Scale data before using lasso or ridge
data<-data[,!names(data) %in% c("C7R2SSCL")]
sds <- apply (data, 2, sd)
matsds <- matrix(rep(sds, times = nrow(data)), nrow(data), byrow = TRUE)
data <- data/matsds
data<-cbind(data, C7R2SSCL = data_o$C7R2SSCL)

### Load relative packages
library(Matrix)
library(foreach)
library(glmnet)

### X and Y matrices
x=model.matrix(C7R2SSCL~.,data)[,-1]  # For use glmnet, construct x design matrix
y=data$C7R2SSCL                           # it can transform qualitative variables into dummy variables

### Ridge
grid =10^ seq (-2,10, length =100)    # Conclude null and full model
ridge.mod = glmnet(x,y,alpha =0,lambda =grid )

dim(coef(ridge.mod)) #29*100 coef for 100 models
l2.norm<-sqrt(apply(coef(ridge.mod)[-1,]^2,2,sum))
plot(log(ridge.mod$lambda),l2.norm,type="l")
plot(log(ridge.mod$lambda),coef(ridge.mod)[2,],type="l",ylim=range(coef(ridge.mod)[-1,]))
for (i in 1:27){points(log(ridge.mod$lambda),coef(ridge.mod)[i+1,],type="l")}
## there is almost no change after log(lambda)=10 which is 
## in magnitude of 10e+4-10e+5
grid = 10^ seq (-7,3, length =100) 
ridge.mod = cv.glmnet(x,y,alpha =0,lambda =grid)
min(ridge.mod$cvm)
ridge.mod$lambda.min
plot(ridge.mod)
cbind("ridge"=coef(ridge.mod, s = "lambda.min"), "full model"=coef(glm1))

### Lasso
grid =10^ seq (-2,10, length =100) 
lasso.mod = glmnet(x,y,alpha =1,lambda =grid)    ## to find the lambda range containing 
l1.norm<-apply(coef(lasso.mod)[-1,],2,sum)         ## both null and full model
plot(log(lasso.mod$lambda),l1.norm,type="l",xlab = "log(lambda)")
plot(log(lasso.mod$lambda),coef(lasso.mod)[2,],type="l",ylim=range(coef(ridge.mod)[-1,]), xlab = "log(lambda)", ylab = "coefficients")
for (i in 1:27){points(log(lasso.mod$lambda),coef(lasso.mod)[i+1,],type="l")}

grid =10^ seq (-4,3, length =100) 
lasso.mod = cv.glmnet(x,y,alpha =1,lambda =grid)
min(lasso.mod$cvm)
lasso.mod$lambda.min
plot(lasso.mod)




###################################### 
# ANOTHER TEST CODED BY OURSELVES ####
######################################

### Generate index for test and train
kfcv.sizes = function(n, k=10) {
  # generate sample sizes for k-fold cross validation on a data set of
  # size n
  
  # author: Matthias C. M. Troffaes
  # date: 22 Nov 2010
  # license: GPLv3
  
  sizes = c()
  for (i in 1:k) {
    first = 1 + (((i - 1) * n) %/% k)
    last = ((i * n) %/% k)
    sizes = append(sizes, last - first + 1)
  }
  sizes
}

kfcv.testing = function(n, k=10) {
  # generate testing sample indices for k-fold cross validation on a
  # data set of size n
  
  # author: Matthias C. M. Troffaes
  # date: 22 Nov 2010
  # license: GPLv3
  
  indices = list()
  sizes = kfcv.sizes(n, k=k)
  values = 1:n
  for (i in 1:k) {
    # take a random sample of given size
    s = sample(values, sizes[i])
    # append random sample to list of indices
    indices[[i]] = s
    # remove sample from values
    values = setdiff(values, s)
  }
  indices
}

### Generate train and test sets
set.seed(1)
test_index<-kfcv.testing(4551)
testset = list()
trainset = list()
for (i in 1:10){
  testset[[i]]<-data[test_index[[i]],]
  train_index <- setdiff(seq(nrow(data)),test_index[[i]])
  trainset[[i]]<-data[train_index,]
}

### Calculate error rate for OLS model
full.mse<-vector()
null.mse<-vector()
for (i in 1:10){
  
  # fit models with trainset
  olsfull.mod<-lm(C7R2SSCL~.,data = trainset[[i]])
  olsnull.mod<-lm(C7R2SSCL~1,data = trainset[[i]])
  
  # predict models use testset
  new = testset [[i]]
  new.x = new[,-29]
  full.pred = predict(olsfull.mod,newdata=new.x)
  null.pred = predict(olsnull.mod,newdata=new.x)
  
  # error
  full<-mean((full.pred-new[,29])^2)
  full.mse<-c(full.mse,full)
  null<-mean((null.pred-new[,29])^2)
  null.mse<-c(null.mse,null)
}
cbind("Full MSE"=mean(full.mse),"Null MSE"=mean(null.mse))

###Part 3: Random Forest
install.packages("randomForest")
library(randomForest)
set.seed(1314)
rf1 <- randomForest(C7R2SSCL ~ ., data = data)
imp <- data.frame(importance(rf1))
impnames <- rownames(imp)
impvals <- as.numeric(unlist(imp))
impnames2 <- impnames[order(impvals, decreasing = TRUE)]
impvals2 <- impvals[order(impvals, decreasing = TRUE)]
barplot(impvals2, names.arg = impnames2, las = 3) 
print(rf1)

### Estimate estimate OOB error rate for the RF fit
### First get the predicted values based on the OOB units
preds <- rf1$predicted
### Then calculate the OOB MSPE
Y <- data$C7R2SSCL
(MSPE <- mean((preds - Y)^2))
