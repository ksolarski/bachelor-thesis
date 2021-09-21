# Bachelor Thesis

This repo contains my bachelor thesis that was written as a part of my undergraduate degree in Econometrics at Maastricht University. Both paper and the code written in R can be found here. The abstract is provided below.

# Abstract

I investigate the most common problems encountered in binary dependent variable models and proposed solutions. In particular, I examine the bias and inconsistency of the following estimators: linear probability model, non-linear fixed effects, and non-linear random effects. I discuss the bias correction for the non-linear fixed effects model proposed by Fernández-Val (2009) and a solution to the initial conditions problem proposed by Wooldridge (2005). Consequently, I conduct a Monte Carlo simulation to compare the performance of various estimation techniques. The simulation indicates that for the estimations of the coefficients one should use the random effects estimator for small values of T and be indifferent between the random effects estimator and the bias-corrected fixed effects estimator for higher values of T. When marginal effects are the quantity of interest, one should use either the fixed effects estimator or the bias-corrected fixed effects estimator.

# Citations

<details>
    <summary>
      <em>"Simple solutions to the initial conditions problem in dynamic, nonlinear panel data models with unobserved heterogeneity." (Wooldrdge, 2005)</em>
    </summary>
    <br/>
    <pre>
@article{wooldridge2005simple,
  title={Simple solutions to the initial conditions problem in dynamic, nonlinear panel data models with unobserved heterogeneity},
  author={Wooldridge, Jeffrey M},
  journal={Journal of applied econometrics},
  volume={20},
  number={1},
  pages={39--54},
  year={2005},
  publisher={Wiley Online Library}
}
</pre>
  </details>
  
<details>
    <summary>
      <em>"Fixed effects estimation of structural parameters and marginal effects in panel probit models" (Fernández-Val, 2009)</em>
    </summary>
    <br/>
    <pre>
@article{fernandez2009fixed,
  title={Fixed effects estimation of structural parameters and marginal effects in panel probit models},
  author={Fern{\'a}ndez-Val, Iv{\'a}n},
  journal={Journal of Econometrics},
  volume={150},
  number={1},
  pages={71--85},
  year={2009},
  publisher={Elsevier}
}
</pre>
  </details>
