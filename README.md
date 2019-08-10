# mbQuantile

Multiplier Bootstrap for Quantile Regression

## Description

In this project, we provide simulation code for the paper "Multiplier bootstrap for quantile regression: non-asymptotic theory under random design". See [this link](https://www.math.ucsd.edu/~xip024/QR_Boot.pdf) for the paper and [this link](https://www.math.ucsd.edu/~xip024/QR_Supp.pdf) for the supplement. The code can reproduce numerical results in Section 3 and Appendix B.

Specifically, to duplicate results in Section 3.1 and Appendix B.1, run the file **mb_ci.R**, to replicate results in Section 3.2 and Appendix B.2, run the file **mb_ht.R**. In both files, we allow the following various settings, with details stated in the paper:

* Model type: homoscedastic model / heteroscedastic model
* Error distribution: student's t / normal mixture type I / normal mixture type II
* Covariates design: independent / weakly correlated / equally correlated

## Authors

Xiaoou Pan <xip024@ucsd.edu>, Wen-Xin Zhou <wez243@ucsd.edu> 

## Main reference

Chen, K., Ying, Z., Zhang, H. and Zhao, L. (2008). Analysis of least absolute deviation. Biometrika 95 107–122. [Paper](https://academic.oup.com/biomet/article-abstract/95/1/107/219099)

Feng, X., He, X. and Hu, J. (2011). Wild bootstrap for quantile regression. Biometrika 98 995–999. [Paper](https://academic.oup.com/biomet/article-abstract/98/4/995/234840)

Koenker, R. (2005). Quantile Regression. Cambridge Univ. Press, Cambridge. [Book](https://www.cambridge.org/core/books/quantile-regression/C18AE7BCF3EC43C16937390D44A328B1)

Koenker, R. (2018). Package ‘quantreg’. [CRAN](https://cran.r-project.org/web/packages/quantreg/index.html)

Koenker, R. and Bassett, G. (1978). Regression quantiles. Econometrica 46 33-50. [Paper](https://www.jstor.org/stable/1913643?seq=1#metadata_info_tab_contents)

Pan, X. and Zhou, W.-X. (2019). Multiplier bootstrap for quantile regression: non-asymptotic theory under random design. [Paper](https://www.math.ucsd.edu/~xip024/QR_Boot.pdf) | [Supplement](https://www.math.ucsd.edu/~xip024/QR_Supp.pdf)

Parzen, M. I., Wei, L. J. and Ying, Z. (1994). A resampling method based on pivotal estimating functions. Biometrika 81 341–350. [Paper](https://academic.oup.com/biomet/article-abstract/81/2/341/468184)

