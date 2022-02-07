# -Analysis-of-variance-ANOVA-
The aim of the analysis is usually to investigate whether the mean of the numeric #variable differs significantly between the groups and to quantify the divergence.
with(InsectSprays, tapply(count, spray, length)) #12 obs. per group
with(InsectSprays, tapply(count, spray, summary)) #compare distr.
with(InsectSprays, tapply(count, spray, sd)) #standard deviations
library(ggplot2)
ggplot(InsectSprays, aes(spray, count)) +
   geom_point(shape = 1, position = position_jitter(width = 0.2 +
     height = 2))
	 
ggplot(InsectSprays, aes(spray, count))+
  geom_point(shape = 1, position = position_jitter(width = 0.2))

library(lattice)
stripplot(count~spray, data=InsectSprays, jitter=0.2)
plot(count ~ spray, data=InsectSprays, las=1)

install.packages(sciplot)
library(sciplot)
bargraph.CI(spray, count, col=(gray(0.88)), data=InsectSprays,
               xlab="spray", ylab="count", ylim=c(0,20))
lineplot.CI(spray, count, type="p", data = InsectSprays,
               xlab="spray", ylab="count", ylim=c(0,20))

## Parametric models: one-way ANOVA 
ins.lm <- lm(count ~ spray, data=InsectSprays)
anova(ins.lm)

library(car)
qqPlot(resid(ins.lm))
qqPlot(rnorm(72))

shapiro.test(resid(ins.lm)) 

plot(jitter(fitted(ins.lm),1), resid(ins.lm), las=1,
      xlab="Jittered fitted values", ylab="Residuals")
     + abline(h=0)

bartlett.test(count ~ spray, data=InsectSprays)
leveneTest(count ~ spray, data=InsectSprays)
oneway.test(count ~ spray, data=InsectSprays)

par(mfrow=c(2,2))
plot(ins.lm)
par(mfrow=c(1,1))
with(InsectSprays, pairwise.wilcox.test(count, spray, "holm"))
a <- aov(count~spray, data=InsectSprays)
TukeyHSD(a, "spray")
library(agricolae)
HSD.test(ins.lm, "spray", group=TRUE, console=TRUE)
library(multcomp)
summary(glht(ins.lm, mcp(spray="Dunnett")))

InsectSprays$b.ref <- relevel(InsectSprays$spray, ref="B")
summary(glht(lm(count~b.ref, data=InsectSprays), mcp(b.ref="Dunnett")))

InsectSprays$b.ref <- relevel(InsectSprays$spray, ref="B")
summary(glht(lm(count~b.ref, data=InsectSprays), mcp(b.ref="Sequen")))

#RMSE and R2
anova(ins.lm)
summary(ins.lm)$sigma
summary(ins.lm)$r.squared

# the reference level comparisons
summary(ins.lm)

# look at the sample group means again:
with(InsectSprays, tapply(count,spray,mean))

## to remove the intercept and get  the
# estimate mean for every category, use the term - 1
summary(lm(count ~ spray - 1 , data = InsectSprays))$coefficients

kruskal.test(count ~ spray, data=InsectSprays)
library(NSM3)
with(InsectSprays, pSDCFlig(x=count, g=as.numeric(spray), method=NA))


## 1.5 Robust method 
library(asbio)
with(InsectSprays, BDM.test(count, spray))
with(InsectSprays, trim.test(count, spray))
groupmeans <- c(120, 130, 140, 150)
power.anova.test(groups = length(groupmeans)
                   + between.var = var(groupmeans)
                   + within.var = 500, power = .80, sig.level = 0.05)

##The two-way layout, The turnip yield data
library(agridat)
turnip <- mcconway.turnip
turnip$density <- as.factor(turnip$density)
turnip

ggplot(turnip, aes(x=density, y=yield, fill=gen)) + geom_boxplot()
turnip[c(1:4,61:64),]

qplot(density, yield, data=turnip, facets=~block, shape=gen, size=I(2))

turnip$density <- as.factor(turnip$density)
with(turnip, interaction.plot(density, gen, yield))

turnip.full <- lm(yield ~ gen*density, data=turnip)
coef(turnip.full)

with(turnip, tapply(yield, list(gen, density), mean))

coef(lm(yield ~ gen*density, data=turnip))
 
(mean(turnip$yield))
#Ftest
summary(turnip.full)

anova(turnip.full)

summary(turnip.full)

library(car)
qqPlot(resid(turnip.full), xlab="Normal quantiles", ylab="Residuals")

shapiro.test(resid(turnip.full))	#null hypothesis rejecetd
 
turnip.sqrt <- lm(sqrt(yield) ~ gen*density, data=turnip)
qqPlot(resid(turnip.sqrt), xlab="Normal quantiles", ylab="Residuals")

shapiro.test(resid(turnip.sqrt))$p.value

plot(fitted(turnip.sqrt), resid(turnip.sqrt), xlab="Fitted values",
      ylab="Residuals"); abline(h=0)
bartlett.test(sqrt(yield) ~ interaction(gen,density), data=turnip)
leveneTest(sqrt(yield) ~ interaction(gen,density), data=turnip)

library(car)
Anova(turnip.sqrt, white.adjust = "hc3")

par(mfrow=c(2,2))
plot(turnip.sqrt)
par(mfrow=c(1,1))

coef(turnip.sqrt)
dat <- expand.grid(gen = c("Marco", "Barkant"),
                    density = as.factor(c(1,2,4,8)))
dat
pred <- predict(turnip.sqrt, dat)
pred
cbind(dat, pred.yield = pred^2)

anova(lm(Y1 ~ Loc + Var, immer))

anova(lm(Y1 ~ Var, immer))
# no significan effect
