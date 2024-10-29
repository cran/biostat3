## Date: 2015-03-03
## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist
## Modified: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 7
###############################################################################

## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("muhaz")
## install.packages("car")

## @knitr loadDependencies
library(biostat3)
library(dplyr)    # for data manipulation
library(car)      # for car::linearHypothesis -> biostat3::lincom

## @knitr loadPreprocess
## Read melanoma data

## Create a new dataset with only localised cancer
melanoma.l <- filter(biostat3::melanoma, stage=="Localised")
head( melanoma.l )
summary(melanoma.l)

## @knitr 7.a.i

## Plot Kaplan-Meier curve using survfit
## Create a new event indicator
melanoma.l <- mutate(melanoma.l,
                   death_cancer = as.numeric(status=="Dead: cancer") )

## Create a fitted object for our subcohort
## using survfit
sfit7a1 <- survfit(Surv(surv_mm, event=death_cancer) ~ year8594,
                   data = melanoma.l )

## Have a look at the fitted object
str(sfit7a1, 1)

## Plot the survival curve (with some bells and whistles)
plot(sfit7a1,
     ## No automatic labelling of the curve (we do that ourselves)
     mark.time=F,
     ## Time is measured in months,  but we want to see it in years
     xscale=12,
     ## Make the plot prettier
     xlab="Years since diagnosis",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))
## Add legend too
legend("bottomleft",legend=levels(melanoma.l$year8594),col=c("blue","red"),lty=c("solid","dashed"), bty="n")

### TRY IF YOU WANT ###
if (FALSE) {
    library(survMisc)
    ## Note: `autoplot(sfit7a1)` was broken; I have submitted a pull request to fix this
    ## autoplot(sfit7a1)
    ## alternatively:
    autoplot(sfit7a1, timeTicks = "custom", times= seq(0, 20*12, 5*12))
}

## @knitr 7.a.ii

## To plot smoothed hazards, we use the muhaz package (using the muhaz2 wrapper)

plot(muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594, data=melanoma.l),
     xlab="Years since diagnosis", col=c("blue","red"), lty=1:2)


## @knitr 7.a.iii
## Compare with Kaplan-Meier plot
par(mfrow=1:2)
plot(sfit7a1,
     ## No automatic labelling of the curve (we do that ourselves)
     mark.time=FALSE,
     ## Time is measured in months,  but we want to see it in years
     xscale=12,
     ## Make the plot prettier
     xlab="Years since diagnosis",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))

plot(muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594, data=melanoma.l),
     xlab="Years since diagnosis", col=c("blue","red"),lty=c("solid","dashed"))

## @knitr 7.b

survRate(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma.l)


## @knitr 7.c.i

## Calculate the incidence rate by time of diagnosis
## but with new variables
melanoma.l2 <- mutate(melanoma.l,
                      ## Update the death indicator (only count deaths within 120 months)
                      ## death_cancer = death_cancer * as.numeric(surv_mm<=120),
                      death_cancer = ifelse(surv_mm<=120 & status == "Dead: cancer",1,0),
                      ## Create a new time variable
                      ## surv_mm = pmin(surv_mm, 120)
                      surv_mm = ifelse(surv_mm<=120, surv_mm, 120)
                      )

## Calculate the rates on the truncated data
rates_by_diag_yr2 <- survRate(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma.l2)
rates_by_diag_yr2

## @knitr 7.c.ii

## MRR full data
rates_by_diag_yr2[2, "rate"] / rates_by_diag_yr2[1, "rate"]
with(rates_by_diag_yr2[2:1,], poisson.test(event, tstop))

## @knitr 7.c.iii
## Use glm to estimate the rate ratios
## we scale the offset term to 1000 person-years
poisson7c <- glm( death_cancer ~ year8594 + offset( log( surv_mm/12/1000 ) ), family=poisson, data=melanoma.l2 )
summary( poisson7c )

## also for collapsed data
summary(glm( event ~ year8594 + offset( log( tstop/12/1000 ) ), family=poisson, data=rates_by_diag_yr2))


## IRR
eform(poisson7c)


## Note that the scaling of the offset term only has an impact on the intercept
summary( glm( death_cancer ~ year8594 + offset( log( surv_mm ) ),
             family=poisson, data=melanoma.l2 ) )

## @knitr 7.d

## Add a new variable for year
melanoma.l2 <- mutate( melanoma.l2, surv_yy1 = surv_mm/12)

## Split follow up by year
melanoma.spl <- survSplit(melanoma.l2, cut=0:9, end="surv_yy1", start="start",
                           event="death_cancer")

## Calculate persontime and
## recode start time as a factor
melanoma.spl <- transform(melanoma.spl,
                          pt = surv_yy1 - start,
                          fu = as.factor(start) )

## @knitr 7.e

## Calculate the incidence rate by observation year
yearly_rates <- survRate(Surv(pt/1000,death_cancer)~fu, data=melanoma.spl) |>
    transform(start=as.numeric(levels(fu))[fu]) |>
    transform(mid=start+0.5)

yearly_rates2 <- survRate(Surv(pt/1000,death_cancer)~fu+sex, data=melanoma.spl) |>
    transform(start=as.numeric(levels(fu))[fu]) |>
    transform(mid=start+0.5)

## Plot by year
with(yearly_rates,
     matplot(fu,
             cbind(rate, lower, upper),
             lty=c("solid","dashed","dashed"),
             col=c("black","gray","gray"),
             type="l",
             main="Cancer deaths by years since diagnosis",
             ylab="Incidence rate per 1000 person-years",
             xlab="Years since diagnosis") )

ci_polygon = function(x, lower, upper, ...)
    polygon(c(x,rev(x)), c(lower, rev(upper)), ...)
with(yearly_rates, {
    matplot(start, cbind(lower, upper), type="n",
            main="Cancer deaths by years since diagnosis",
            ylab="Incidence rate per 1000 person-years",
            xlab="Years since diagnosis")
    ci_polygon(start,lower,upper,
               col="lightgrey",
               border="transparent")
    lines(start,rate)
})

library(tinyplot)
with(yearly_rates, {
    plt(rate~start, ymin=lower, ymax=upper, type="ribbon",
        main="Cancer deaths by years since diagnosis",
        ylab="Incidence rate per 1000 person-years",
        xlab="Years since diagnosis")
})
with(yearly_rates2, {
    plt(rate~start|sex, ymin=lower, ymax=upper, type="ribbon",
        col = c(blue = "#0072B2", orange = "#E69F00"),
        main="Cancer deaths by years since diagnosis",
        ylab="Incidence rate per 1000 person-years",
        xlab="Years since diagnosis")
})


library(scales) # alpha() for transparency
ci_polygon = function(x, lower, upper, ...)
    polygon(c(x,rev(x)), c(lower, rev(upper)), ...)
with(yearly_rates2, {
    col = c(blue = "#0072B2", orange = "#E69F00")
    group = levels(sex)
    matplot(start, cbind(lower, upper), type="n",
            main="Cancer deaths by years since diagnosis",
            ylab="Incidence rate per 1000 person-years",
            xlab="Years since diagnosis")
    for (i in 1:length(group)) {
        index = sex == group[i]
        ci_polygon(start[index],lower[index],upper[index],
                   col=alpha(col[i],0.4),
                   border="transparent")
        lines(start[index],rate[index],col=col[i])
    }
    legend("topright", legend=group, lty=1, col=col, bty="n")
})

library(Hmisc)
xYplot(Cbind (rate, lower, upper) ~ start,
       groups=sex, 
       method="filled bands",
       data=yearly_rates2,
       type="l",
       keys='lines',
       col.fill=alpha(trellis.par.get("superpose.line")$col,0.4),
       main="Cancer deaths by years since diagnosis",
       ylab="Incidence rate per 1000 person-years",
       xlab="Years since diagnosis")

library(tactile)
with(yearly_rates2,
     xyplot(rate~start,
            lower=lower,
            upper=upper,
            type="l",
            groups=sex,
            auto.key=TRUE,
            panel=function(...) {tactile::panel.ci(...); panel.xyplot(...)},
            main="Cancer deaths by years since diagnosis",
            ylab="Incidence rate per 1000 person-years",
            xlab="Years since diagnosis"))

#' @description Superpose plots for grouped data
#' @param x vector for the x axis
#' @param y vector for the estimated values
#' @param INDEX vector for the index
#' @param panel a function with signature function(x, y, i, subscripts, value), where x and y are subsets of the full x and y, where i is a counter for the index set, subscripts is an indicator vector for the subset and value is the ith INDEX value
#' @param ... other parameters passed to `plot`

make_panel = function(FUN,col=1:7,lty=1:7,pch=1:7,...)
    function(x,y,i,subscripts,value)
        FUN(x,y,col=col[i],lty=lty[i],pch=pch[i],...)
plot.superpose = function(x, y, INDEX, panel=make_panel(points), index=NULL,
                          legend=unique(INDEX), legend.args=list(), add=FALSE, ...) {
    if (!add) plot(x,y,...,type="n")
    uids = unique(INDEX)
    for(i in 1:length(uids)) {
        subscripts = INDEX == uids[i]
        panel(x[subscripts], y[subscripts], i=i,
              subscripts=subscripts, value=uids[i])
    }
    if (!add)
        do.call("legend", modifyList(list(x="topright",legend=legend), legend.args))
}
ci_polygon = function(x, lower, upper, border="transparent", ...)
    polygon(c(x,rev(x)), c(lower, rev(upper)), border=border, ...)
with(yearly_rates2, {
    col = trellis.par.get("superpose.line")$col
    plot.superpose(start,rate,sex,
                   legend.args=list(col=col,lty=1,bty="n"),
                   ylim=c(0,80),
                   panel=function(x,y,i,subscripts,...) {
                       ci_polygon(x, lower[subscripts],upper[subscripts],
                                  col=alpha(col[i],0.4))
                       lines(x,y,col=col[i])
                   })
})

## matrix y
make_panel = function(FUN,col=1:7,lty=1:7,pch=1:7,...)
    function(x,y,i,value)
        FUN(x,y,col=col[i],lty=lty[i],pch=pch[i],...)
plot.superpose = function(x, y, INDEX, panel=make_panel(points), index=NULL,
                          legend=unique(INDEX), legend.args=list(), add=FALSE, ...) {
    if (!add) matplot(x,y,...,type="n")
    uids = unique(INDEX)
    for(i in 1:length(uids)) {
        subscripts = INDEX == uids[i]
        panel(x[subscripts], y[subscripts,], i=i, value=uids[i])
    }
    if (!add)
        do.call("legend", modifyList(list(x="topright",legend=legend), legend.args))
}
ci_polygon = function(x, lower, upper, border="transparent", ...)
    polygon(c(x,rev(x)), c(lower, rev(upper)), border=border, ...)
with(yearly_rates2, {
    col = trellis.par.get("superpose.line")$col
    plot.superpose(start, cbind(rate,lower,upper), sex,
                   legend.args=list(col=col,lty=1,bty="n"),
                   panel=function(x,y,i,...) {
                       ci_polygon(x, y[,2],y[,3],
                                  col=alpha(col[i],0.4))
                       lines(x,y[,1],col=col[i])
                   })
})

model.matrix(cbind(rate,lower,upper)~start+sex+1, data=yearly_rates2)[,-1]

## matrix y (simplified)
matplot2 = function(x, y, group=NULL, panel=matpoints,
                    xlab=NULL, ylab=NULL, ...) {
    stopifnot(length(x)==NROW(y),
              is.null(group) || length(x)==length(group),
              is.matrix(y) || is.vector(y))
    if (is.null(group)) group=1 # arbitrary value
    if (is.null(xlab)) xlab=deparse(substitute(x))
    if (is.null(ylab)) ylab=deparse(substitute(y))
    if (is.vector(y)) y=matrix(y)
    group_levels = if (is.factor(group)) levels(group) else unique(group)
    matplot(x, y, ..., type="n", xlab=xlab, ylab=ylab)
    for(i in 1:length(group_levels)) {
        subscripts = group == group_levels[i]
        panel(x[subscripts], y[subscripts,], i=i, group=group_levels[i])
    }
}
ci_polygon = function(x, lower, upper, border="transparent", ...)
    polygon(c(x,rev(x)), c(lower, rev(upper)), border=border, ...)
with(yearly_rates, {
    matplot2(start, cbind(rate,lower,upper), 
             panel=function(x,y,...) {
                 ci_polygon(x, y[,2], y[,3],
                            col="lightgrey")
                 lines(x,y[,1],col=1)
             },
             main="Cancer deaths by\nyears since diagnosis",
             ylab="Incidence rate per 1000 person-years",
             xlab="Years since diagnosis")
})
with(yearly_rates2, {
    col = trellis.par.get("superpose.line")$col
    matplot2(start, cbind(rate,lower,upper), sex,
             panel=function(x,y,i,...) {
                 ci_polygon(x, y[,2], y[,3],
                            col=alpha(col[i],0.4))
                 lines(x,y[,1],col=col[i])
             },
             main="Cancer deaths by\nyears since diagnosis",
             ylab="Incidence rate per 1000 person-years",
             xlab="Years since diagnosis")
    legend("topright", legend=levels(sex), lty=1, col=col, bty="n")
})

with(yearly_rates2, {
    matplot2(start, cbind(rate,lower,upper), sex,
             panel=function(x,y,i,...) matlines(x,y,col=i,lty=c(1,2,2)))
    legend("topright", legend=levels(sex), lty=1, col=1:2, bty="n")
})

with(yearly_rates, {
    par(mfrow=1:2)
    matplot(start, cbind(rate,lower,upper), type="l", main="matplot",
            lty=c(1,2,2), col=1)
    matplot2(start, cbind(rate,lower,upper),
             panel=function(x,y,...) matlines(x,y,lty=c(1,2,2),col=1),
             main="matplot2")
})


## Example using coplot()
panel.ci = function(x, y, lower, upper, subscripts, col, pch,
                    border="transparent", fill.col="grey",
                    line.col="black", ...) {
    polygon(c(x,rev(x)),
            c(lower[subscripts], rev(upper[subscripts])),
            border=border, col=fill.col)
    lines(x,y,col=line.col)
}


## @knitr 7.f
# Plot smoothed hazards

par(mfrow=1:2)
with(yearly_rates, matplot(as.numeric(as.character(fu))+0.5,
                           cbind(rate, lower,
                                 upper),
                           ylim=c(0,max(upper)),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths by time since diagnosis",
                           ylab="Mortality rate per 1000 person-years",
                           xlab="Years since diagnosis") )

library(bshazard)
hazfit7f <- bshazard(Surv(surv_mm/12, status == "Dead: cancer") ~ 1, data = melanoma.l)
## scale hazard by 1000
plot(hazfit7f, xlab="Years since diagnosis",col="blue",lty="solid")


## @knitr 7.g
## Run Poisson regression
summary(poisson7g <- glm( death_cancer ~ fu + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
eform(poisson7g)

## @knitr 7.h
summary(poisson7h <- glm( death_cancer ~ fu + year8594 + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
eform(poisson7h)

# Add interaction term
summary(poisson7h2 <- glm( death_cancer ~ fu*year8594 + offset( log(pt) ), family=poisson, data=melanoma.spl ))
## IRR
eform(poisson7h2)

## @knitr 7.i

summary(poisson7i <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
eform(poisson7i)

## Test if the effect of age is significant using a likelihood ratio test
drop1(poisson7i, ~agegrp, test="Chisq")
## For this we can also use the car package and a Wald test
linearHypothesis(poisson7i,c("agegrp45-59 = 0","agegrp60-74 = 0","agegrp75+ = 0"))
## ADVANCED:
## Alternative approach for the likelihood ratio test
# poisson7i_2 <- update(poisson7i,. ~ . - agegrp)
# anova(poisson7i_2,poisson7i,test="Chisq")

## @knitr 7.j

summary(poisson7j <- glm( death_cancer ~ fu + agegrp + year8594*sex + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
eform(poisson7j)

## @knitr 7.k.i
# hand calculations
hz7k <- exp(coef(poisson7j))
hz7k["sexFemale"]
hz7k["sexFemale"]*hz7k["year8594Diagnosed 85-94:sexFemale"]

## @knitr 7.k.ii
## You will need the "car" package to use lincom. If it is not already installed:
## install.packages("car")
lincom(poisson7j,c("sexFemale + year8594Diagnosed 85-94:sexFemale"),eform=TRUE)


## @knitr 7.k.iii
## Create dummies and Poisson regression
melanoma.spl <- melanoma.spl %>%
    ## Add confidence intervals for the rates
    mutate(femaleEarly = sex=="Female" & year8594=="Diagnosed 75-84",
           femaleLate = sex=="Female" & year8594=="Diagnosed 85-94")

summary(poisson7k <- glm( death_cancer ~ fu + agegrp + year8594 + femaleEarly +
                         femaleLate + offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))

## IRR
eform(poisson7k)

## @knitr 7.k.iv
## Add interaction term
summary(poisson7k2 <- glm( death_cancer ~ fu + agegrp + year8594 + year8594:sex +
                         offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))
eform(poisson7k2)


## @knitr 7.l

summary( poisson7l.early <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 75-84" ) )
eform(poisson7l.early)

summary( poisson7l.late <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 85-94" ) )
eform(poisson7l.late)

# compare with results in i
eform(poisson7i)

# compare with results in j
eform(poisson7j)


# Poisson-regression with effects specific for diagnose period
summary(poisson7l2 <- glm( death_cancer ~ fu + fu:year8594 + agegrp + agegrp:year8594
                          + sex*year8594 + offset( log(pt) ),
                          family=poisson, data=melanoma.spl ))
eform(poisson7l2)

## @knitr 7.m
## Split follow up by month
library(splines)
time.cut <- seq(0,10,by=1/12)
nrow(biostat3::melanoma)
melanoma.spl <- survSplit(Surv(surv_mm/12,status=="Dead: cancer")~., data=biostat3::melanoma,
                          cut=time.cut,
                          subset=stage=="Localised")
nrow(melanoma.spl)
melanoma.spl <- transform(melanoma.spl, mid=(tstop+tstart)/2, risk_time=tstop-tstart)
poisson7m <- glm(event ~ ns(mid,df=6) + agegrp + year8594 +
                     offset(log(risk_time)),
                 family=poisson,
                 data=melanoma.spl)
df <- data.frame(agegrp="0-44", year8594="Diagnosed 75-84",
                 mid=time.cut[-1], risk_time=1)
## plot the rate at the baseline values
plot(df$mid, predict(poisson7m, newdata=df, type="response"),
     type="l", ylab="Rate", xlab="Time since diagnosis (years)",
     ylim=c(0,0.05))

add_ci = function(df, level=0.95, exponentiate=FALSE,
                  itrans=c("I","log"), ...) {
    stopifnot(is.list(df),
              (".fitted" %in% names(df) && ".se.fit" %in% names(df)) ||
              ("fit" %in% names(df) && "se.fit" %in% names(df)),
              is.numeric(level),
              length(level) %in% c(1,nrow(df)),
              level>0,
              level<1,
              is.logical(exponentiate))
    itrans <- match.arg(itrans)
    if (inherits(df,"list")) df = as.data.frame(df)
    a <- (1 - level)/2
    a <- cbind(a, 1 - a)
    fac <- qnorm(a)
    if (".fitted" %in% names(df) && ".se.fit" %in% names(df))
        within(df, {
            if (itrans=="log") {
                .se.fit <- .se.fit/.fitted
                .fitted <- log(.fitted)
            }
            .upper <- .fitted+fac[,2]*.se.fit
            .lower <- .fitted+fac[,1]*.se.fit
            transf <- function() {
                .upper <<- exp(.upper)
                .lower <<- exp(.lower)
                .fitted <<- exp(.fitted)
                .se.fit <<- .se.fit*.fitted
            }
            if (itrans=="log") transf()
            if (exponentiate) transf()
            transf <- NULL
        })
    else
        within(df, {
            if (itrans=="log") {
                se.fit <- se.fit/fit
                fit <- log(fit)
            }
            conf.high <- fit+fac[,2]*se.fit
            conf.low <- fit+fac[,1]*se.fit
            transf <- function() {
                conf.high <<- exp(conf.high)
                conf.low <<- exp(conf.low)
                fit <<- exp(fit)
                se.fit <<- se.fit*fit
            }
            if (itrans=="log") transf()
            if (exponentiate) transf()
            transf <- NULL
        })
}
predict(poisson7m, newdata=df,se.fit=TRUE) |> add_ci() |> head()
predict(poisson7m, newdata=df,se.fit=TRUE) |> add_ci(exp=TRUE) |> head()
predict(poisson7m, newdata=df,se.fit=TRUE) |> add_ci(exp=TRUE) |> add_ci(itrans="log") |> head()
augment(poisson7m, newdata=df, se_fit=TRUE) |> add_ci()
augment(poisson7m, newdata=df, se_fit=TRUE) |> add_ci(exp=TRUE)
augment(poisson7m, newdata=df, se_fit=TRUE) |> add_ci(exp=TRUE) |> add_ci(itrans="log")


as.data.frame.by <- function(x, row.names=NULL, option=FALSE, ...) {
    out <- array2DF(x, simplify=FALSE) |>
        subset(!sapply(Value,is.null))
    cbind(out[-length(out)], do.call(rbind,out$Value))
}
melanoma.spl4 <-
    by(melanoma.spl2,
       with(melanoma.spl2, data.frame(sex, agegrp, stage, year8594, timeband)),
       function(data)
           with(data,
                data.frame(risk_times=sum(risk_time),
                           events=sum(event),
                           mid=sum(risk_time*mid)/sum(risk_time)))) |>
    as.data.frame()
head(melanoma.spl4)

## @knitr 7.n
## using melanoma.spl and df from previous chunk
poisson7n <- glm(event ~ ns(mid,df=4) + agegrp + year8594 +
                     ifelse(year8594=="Diagnosed 85-94",1,0):ns(mid,df=3) +
                     offset(log(risk_time)),
                 family=poisson,
                 data=melanoma.spl)
library(rstpm2)
library(ggplot2)
## get log(RR) confidence interval using predictnl (delta method)
pred <- predictnl(poisson7n, function(object)
    log(predict(object, newdata=transform(df, year8594="Diagnosed 85-94"),
                type="response") /
        predict(object, newdata=df, type="response")))
pred2 <- transform(pred, time = df$mid, rr = exp(fit), ci = exp(confint(pred)))
ggplot(pred2, aes(x=time, y=rr, ymin=ci.2.5.., ymax=ci.97.5..)) +
    ggplot2::geom_line() + ggplot2::geom_ribbon(alpha=0.6) +
    xlab("Time since diagnosis (years)") +
    ylab("Rate ratio")
    
## Calculate the rate difference
band <- function(x,yy,col="grey")
    polygon(c(x,rev(x)), c(yy[,1], rev(yy[,2])), col=col, border=col)
pred <- predictnl(poisson7n, function(object)
    predict(object, newdata=transform(df, year8594="Diagnosed 85-94"),
                type="response") -
    predict(object, newdata=df, type="response"))
rd <- pred$fit
ci <- confint(pred)
matplot(df$mid,
        ci,
        type="n", # blank plot area
        xlab="Time since diagnosis (years)",
        ylab="Rate difference") 
band(df$mid,ci) # add confidence band
lines(df$mid, rd) # add rate difference

## @knitr 7.o
## Calculate survival from a smooth Poisson regression model 
if (requireNamespace("deSolve")) {
    twoState <- function(object, ...) {
        out <- as.data.frame(markov_msm(list(object),matrix(c(NA,1,NA,NA),2,byrow=TRUE), ...))
        transform(subset(out, state==1),
                  S=P,
                  S.lower=P.lower,
                  S.upper=P.upper)
    }
    df2 <- expand.grid(agegrp=levels(biostat3::melanoma$agegrp),
                       year8594=levels(biostat3::melanoma$year8594))
    df2 <- transform(df2,risk_time=1)
    pred <- twoState(poisson7n, t=c(0,df$mid), newdata = df2, tmvar = "mid")
    ggplot(pred, aes(x=time,y=S,ymin=S.lower,ymax=S.upper,fill=year8594)) +
        ggplot2::geom_line() + ggplot2::geom_ribbon(alpha=0.6) +
        facet_grid(~agegrp) +
        xlab("Time since diagnosis (years)") +
        ylab("Survival")
} else cat("To run this example, please install the deSolve package\n")
    


## Additional random code

## q2.R
hazards <- by(melanoma, melanoma$stage,
              function(data) with(bshazard(Surv(surv_mm/12, death_cancer)~1,
                                           data=data, verbose=FALSE),
                                  data.frame(time,hazard,lower.ci,upper.ci))) |>
    array2DF()
with(hazards, plt(hazard~time|`melanoma$stage`, 
                  ymin=lower.ci, ymax=upper.ci, type="ribbon",
                  ylim=c(0,1), # bug -- not used
                  xlab="Time since diagnosis (years)", ylab="Hazard"))
