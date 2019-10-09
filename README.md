# Clustering-oriented selection of the amount of smoothing in kernel density estimation

Density-based clustering relies on the idea of linking groups to some specific features of the probability distribution underlying the data.
The reference to a true, yet unknown, population structure allows framing
the clustering problem in a standard inferential setting, where the concept of ideal population clustering is defined as the partition induced by the true density function. The nonparametric formulation of this approach,
known as modal clustering, draws a correspondence between the groups and
the domains of attraction of the density modes. Operationally, a nonparametric density estimate is required and a proper selection of the amount of
smoothing, governing the shape of the density and hence possibly the modal
structure, is crucial to identify the final partition. 

We consider kernel density estimation as a tool for the final purpose of modal clustering, addressing the issue from an asymptotic perspective. Two different unimodal modal clustering-oriented bandwidth selectors are derived as the minimizers of the asymptotic version of the Expected Distance in Measure (EDM).

This repository is associated with the article: [Casa, A., Chacón, J.E. & Menardi, G. (2019). *Modal clustering asymptotics with applications to bandwidth selection*](https://arxiv.org/pdf/1901.07300.pdf)

The R package can be installed as follows.

```R
devtools::install_github("AlessandroCasa/BsMc")
```
