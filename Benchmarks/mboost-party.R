load("benchmark.data.sets.RData")

library(mboost)
dlbcl.y <- Surv(dlbcl.surv, dlbcl.event)
dlbcl.x <- data.frame(t(dlbcl.covar))

unix.time(
cf1 <- cforest(dlbcl.y ~ ., data = dlbcl.x, control = cforest_classical(ntree = 10))
          ) ## takes forever; more than a 1 day

unix.time(
          cf2 <- cforest(dlbcl.y ~ ., data = dlbcl.x[, 1:200], control = cforest_classical(ntree = 500))
          ) ## 2.421


unix.time(
          cf3 <- cforest(dlbcl.y ~ ., data = dlbcl.x[, 1:400], control = cforest_classical(ntree = 500))
          ) ## 7 seconds

unix.time(
          cf4 <- cforest(dlbcl.y ~ ., data = dlbcl.x[, 1:1000], control = cforest_classical(ntree = 1000))
          ) ## 72 seconds




x2 <- t(dlbcl.covar - rowMeans(dlbcl.covar))
unix.time(
          gb1 <- glmboost(dlbcl.y ~ ., data = dlbcl.x, family = CoxPH())
##          gb1 <- glmboost(y = dlbcl.y, x = dlbcl.x, family = CoxPH())
          )
unix.time(
          gb2 <- glmboost(x2, dlbcl.y, family = CoxPH())
          )
sum(abs(coef(gb2)) > 0)


gb3 <- glmboost(x2, dlbcl.y, family = CoxPH(), control = boost_control(mstop = 500))
mstop(aic <- AIC(gb3)) ## nope

folds <- cbind(c(rep(1, 120), rep(0, 40)), c(rep(0, 40), rep(1, 120)))

cv2 <- cvrisk(gb2, folds = folds)


bb1 <- blackboost(x2, dlbcl.y, family = CoxPH())
cvbb1 <- cvrisk(bb1, folds = folds)


m1 <- coxph(dlbcl.y ~ ., data = x2[, c(1:4)])
m0 <- coxph(dlbcl.y ~ x2[, c(1)])
