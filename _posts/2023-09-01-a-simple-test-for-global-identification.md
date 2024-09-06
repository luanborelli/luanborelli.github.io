---
layout: post
title: A Simple Test for Global Identification
category: statistics
---

Consider an extremum (point) estimator $\hat\theta$, i.e., the argmax of a sample criterion function $Q_n(\theta)$. Examples include maximum likelihood (ML) and generalized method of moments (GMM) estimators. It is well known that the consistency and asymptotic normality of extremum estimators rely on a set of regularity conditions, which are often difficult to verify or test. One key condition is (global) identification.

A point $\theta_0 \in \Theta \subset \mathbb{R}^p$ is said to be _globally identified_ if the distribution of the model's data at $\theta_0$ differs from that at any other possible $\theta$. From an extremum perspective, $\theta_0$ is globally identified only if

$$\arg\max_{\theta \in \Theta} Q(\theta) = \{\theta_0\},$$

where $Q$ is the population criterion; i.e. $Q_T \xrightarrow[]{p} Q$.

Consider the case of a generalized method of moments estimator. In this case, the criterion function takes the form 

$$
\begin{align}
\label{GMM}
Q(\theta) = - \mathbb{E}[g(Z_t, \theta)]'W\mathbb{E}[g(Z_t, \theta)],
\end{align}
$$

for some *moment functions* $g(Z_t, \theta)$. When the moment functions are linear in $\theta$, the (necessary and sufficient) global identification condition simplifies to an easily verifiable rank condition. However, when the moment functions are nonlinear in $\theta$, specifying primitive conditions for identification becomes much more challenging, as identification then requires conditions for the uniqueness of solutions to nonlinear equations. 

What, then, can be done in practice in the nonlinear case? Surprisingly, a common approach is simply to assume identification without knowing whether it actually holds. As Newey and McFadden (1994) emphasize in their classic handbook chapter:

> A practical "solution" to the problem of global GMM identification, that has often been adopted, is to simply assume identification.

Another approach is to attempt a full characterization of the set of (globally) identified structures and then assess whether, after imposing all reasonable theory-based restrictions, there is reason to believe that the "true structure" lies within this set. Even better, one could statistically test for it.

Bravo, Escanciano, and Otsu (2010) propose a test for the global identification of conditional moment restriction models. These models satisfy moment restrictions of the form

$$\mathbb{E}[h(Z_t, \theta) \mid X_t] = 0,$$

which, by defining $g(Z_t, \theta) \equiv a(X_t) h(Z_t, \theta)$ for some given vector of _instruments_ $a(X_t)$, imply unconditional moment restrictions that underlie estimators of the form $(\ref{GMM})$. The testing problem is formulated as

$$
\begin{align}
    &H_0: P_0 \in \mathcal{P}_u \\ 
    &H_1: P_0 \in \mathcal{P}_c \backslash \mathcal{P}_u,
\end{align}
$$

where $P_0$ is the true unknown probability measure of $Z_t \equiv (Y_t', X_t')'$ — the *true structure* — and $\mathcal{M}$ is the set of all possible measures for $Z_t$ — the *statistical model*. Further,

$$\mathcal{P}_{c} = \left\{P \in \mathcal{M} : \int h(y, \theta_0) ~ dP_{Y|X=x} = 0 \text{ at some unique } \theta_0 \in \Theta \right\},$$

and

$$\mathcal{P}_u = \left\{ P \in \mathcal{P}_c : \int a(x)h(y, \theta) dP = 0 \text{ only at } \theta = \theta_0 \right\}.$$

Under the null hypothesis, $\Theta_I = {\theta_0}$, meaning the parameter $\theta$ is globally identified, and $\hat\theta$ is consistent for $\theta_0$. Under the alternative, $\Theta_I \neq {\theta_0}$, meaning the parameter $\theta$ is _not_ globally identified, and $\hat\theta$ is not consistent.

Bravo, Escanciano, and Otsu (2010) offer the insight that there are estimators consistent for $\theta_0$ _both_ under $H_0$ and $H_1$. Let $\hat\theta_{c}$ be such an estimator. Since $\hat\theta$ is consistent and asymptotically normal only under the null hypothesis, this suggests constructing a Hausman-type test statistic for $H_0$ based on the [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance#:~:text=The%20Hausdorff%20distance%20is%20the,point%20in%20the%20other%20set.) between ${\hat\theta_{c}}$ and $\hat{\Theta}$.

$$T_n \equiv \max_{\theta \in \hat{\Theta}} n(\hat{\theta}_{c} - \theta)' \Sigma_n^{-1} (\hat{\theta}_c - \theta),$$

where $\hat{\Theta}$ is some consistent (set) estimator for the identified set, and $\Sigma_{n}$ is a consistent estimator for the asymptotic variance matrix $\Sigma$ of $\sqrt{T}(\hat{\theta}_{c} - \hat{\theta})$.

For the estimator $\hat{\theta}_c$, Bravo, Escanciano and Otsu (2010) adopt Dominguez and Lobato's (2004) estimator, based on the following characterization of the conditional moments by an infinite number of unconditional moments: 

$$\mathbb{E}[h(Y_t, \theta_0)|X_t] = 0 \iff \mathbb{E}[h(Y_t, \theta_0) I(X_t \leq x)] = 0 \quad \forall x \in \mathbb{R}^{d_x}.$$


Thus, letting $H(\theta, x) \equiv h(Y_t, \theta) I(X_t \leq x)$, $\theta_0$ satisfies 

$$\theta_0 = \arg\min_{\theta \in \Theta} \int H(\theta, x)^2 ~ dP_{X_t}(x).$$

This suggests 

$$\hat{\theta}_c = \arg \min_{\theta \in \Theta} \frac{1}{n^3} \sum_{l = 1}^n \left( \sum_{t=1}^n h(Y_t, \theta) I(X_t \leq X_l) \right)^2.$$


Finally, for the set estimator $\hat{\Theta}$, Bravo, Escanciano, and Otsu (2010) consider a straightforward contour set estimator similar to the one proposed by Chernozhukov, Hong, and Tamer (2007) and introduce a simple numerical algorithm for its computation. The procedure involves generating $m \in \mathbb{N}$ random, independent initial guesses, using numerical optimization methods to compute the $m$ minima $Q_n(\hat{\theta}^{(j)})$ for $j = 1, \dots, m$, and defining

$$\hat{\Theta} = \{\hat{\theta}^{(j)} ~ : ~ Q_n(\hat{\theta}^{(j)}) \leq \min_{1 \leq j \leq m} Q_n(\hat{\theta}^{(j)}) + a_n\},$$


where $a_n = 1/(n \log(n))$.

Bravo, Escanciano, and Otsu (2010) show that, under mild regularity conditions, the test statistic $T_n$ has a limiting chi-squared distribution with $p$ degrees of freedom under the null hypothesis; i.e., $T_n \xrightarrow[]{p} \chi_p^2$. While they do not fully characterize the asymptotic power properties of the test for general alternatives, their Monte Carlo experiments suggest that the proposed test performs well in finite samples, exhibiting strong size and power properties under both the alternatives of lack of identification and weak identification.


---
{: data-content=" references"}

1. Bravo, F., Escanciano, J.C. and Otsu, T., 2012. A simple test for identification in GMM under conditional moment restrictions. In _Essays in Honor of Jerry Hausman_ (Vol. 29, pp. 455-477). Emerald Group Publishing Limited.
2. Chernozhukov, V., Hong, H. and Tamer, E., 2007. Estimation and confidence regions for parameter sets in econometric models 1. _Econometrica_, _75_(5), pp.1243-1284.
3. Newey, W.K. and McFadden, D., 1994. Large sample estimation and hypothesis testing. _Handbook of econometrics_, _4_, pp.2111-2245.


---
{: data-content="//"}
