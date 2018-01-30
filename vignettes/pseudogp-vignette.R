## ----setup, include = FALSE----------------------------------------------
library(ggplot2)
theme_set(theme_bw())
knitr::opts_chunk$set(echo = TRUE, cache = TRUE, 
                      message = FALSE, warning = FALSE,
                      fig.center = TRUE, fig.width = 6, fig.height = 4)

## ----simple-pca, warning = FALSE, message = FALSE, fig.width = 5, fig.height=3----
library(pseudogp)
library(ggplot2)
data(monocle_le)

ggplot(data.frame(monocle_le)) + geom_point(aes(x = X1, y = X2)) + theme_bw()

## ----simple-fit, warning = FALSE, message = FALSE, cache=TRUE, eval=FALSE----
#  le_fit <- fitPseudotime(monocle_le, smoothing_alpha = 30, smoothing_beta = 6, iter = 1000, chains = 1)

## ----load-le-fit---------------------------------------------------------
data(le_fit)

## ----poscurve, warning = FALSE, fig.width = 5, fig.height = 3------------
posteriorCurvePlot(monocle_le, le_fit)

## ----ind-trace, fig.width = 5, fig.height = 3----------------------------
posteriorCurvePlot(monocle_le, le_fit, posterior_mean = FALSE)

## ----simple-boxplot, fig.width = 6, fig.height = 4-----------------------
posteriorBoxplot(le_fit)

## ----rstan-plotting, fig.width=6, fig.height=4---------------------------
rstan::plot(le_fit, pars="lambda")
rstan::traceplot(le_fit, pars="sigma")

## ----extract-pst---------------------------------------------------------
pst <- rstan::extract(le_fit, pars="t")$t
print(str(pst))

## ----posplot, fig.width=7, fig.height=3, message=FALSE-------------------
set.seed(1L)
to_sample <- sample(nrow(monocle_le), 4)
traces <- data.frame(pst[,to_sample])
names(traces) <- paste0("Cell",1:4)
traces_melted <- reshape2::melt(traces, variable.name="Cell", value.name="Pseudotime")
pstplt <- ggplot(traces_melted) + geom_density(aes(x = Pseudotime, fill = Cell), alpha = 0.5) +
  theme_bw()
pstplt

## ----posmode, fig.width=7, fig.height=3----------------------------------
library(MCMCglmm)
library(coda)
tmcmc <- mcmc(traces)
tmap <- posterior.mode(tmcmc)
for(i in 1:length(tmap)) {
  pstplt <- pstplt + geom_vline(xintercept = tmap[i], linetype = 2)
}
pstplt

## ----curve-examples, fig.width=6, fig.height=3---------------------------
x <- runif(100, -1, 1)
y_structured <- rnorm(100, x^2, sd = 0.1)
y_unstructured <- rnorm(100, x^2, sd = 1)
dfplt <- data.frame(x, y_structured, y_unstructured)
dfmelt <- reshape2::melt(dfplt, id.vars = "x", value.name = "y")
ggplot(dfmelt) + geom_point(aes(x=x, y=y)) + facet_wrap(~variable) + theme_bw()

## ----firstfit, message = FALSE, warning = FALSE, fig.width=7, fig.height=4----
## first let's load the data
data(monocle_le)
data(monocle_pca)
data(monocle_tsne)

X <- data.frame(rbind(monocle_le, monocle_pca, monocle_tsne))
names(X) <- c("x1", "x2")
X$cell <- rep(1:nrow(monocle_le), times = 3)
X$representation <- rep(c("LE", "PCA", "tSNE"), each = nrow(monocle_le))
ggplot(X) + geom_point(aes(x = x1, y = x2)) + facet_wrap( ~ representation) + theme_bw()

## ----stanfit, cache=TRUE, eval=FALSE-------------------------------------
#  data <- list(LE = monocle_le, PCA = monocle_pca, tSNE = monocle_tsne)
#  fit <- fitPseudotime(data, chains = 2, iter = 1000, smoothing_alpha = 12, smoothing_beta = 3)

## ----sessioninfo---------------------------------------------------------
devtools::session_info()

