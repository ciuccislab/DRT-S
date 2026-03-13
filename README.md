# DRT-S
This repository contains some of the source code used for the paper titled 'Entropy-based regularized regression for advanced distribution of relaxation times deconvolution' in the Journal of Power Sources, 644 (2025) 236910. https://doi.org/10.1016/j.jpowsour.2025.236910. The article is available online and in the docs folder.

Electrochemical impedance spectroscopy (EIS) is an experimental technique often leveraged to study the electrochemical behavior of batteries and fuel cells to cite only two [1-2]. The success of this method is related to its simplicity, broad applicability, and non-invasiveness, yet the analysis of EIS spectra still remains a cumbersome quest [3-4]. In this context, the distribution of relaxation times (DRT) has become a well-known approach for circumventing the limitations of physical models and equivalent circuits [1-2]. To solve the ill-posed inverse problem of reconstructing the DRT from the measured impedances, Tikhonov regularization or ridge regression is often used to obtain the DRT by minimizing the penalized distance between the experimental and DRT-based impedances [5-8]. In this article, we delve into an alternative regularization approach, namely the entropy, to demonstrate its superiority compared to ridge regression as the former prevents the appearance of meaningless DRT peaks amongst other attractive properties [9].

Dependencies
numpy

scipy

matplotlib

time

cvxopt

cxpy

bayes_opt

pymc

Tutorials
norm-DRT_distant2ZARC.ipynb: this notebook shows how to deconvolve with different norms the DRT of impedance data artificially generated using two ZARCs in series;
norm-DRT_SOFC.ipynb : this notebook shows how to recover with different norms the DRT of real impedance EIS data measured by our group on a solid oxide fuel cell.

Citation
@article{py2025entropy,
  title={Entropy-based regularized regression for advanced distribution of relaxation times deconvolution},
  author={Py, Baptiste and Wang, Zilong and Wang, Yuhao and Ciucci, Francesco},
  journal={Journal of Power Sources},
  volume={644},
  pages={236910},
  year={2025},
  publisher={Elsevier}
}

References
[1] Z. Wang, Y. Wang, F. Ciucci, Distribution of relaxation times: Foundations, methods, diagnostics, and prognosis for electrochemical systems, Curr. Opin. Electrochem. (2025) 101789. https://doi.org/10.1016/j.coelec.2025.101789.

[2] A. Maradesa, B. Py, F. Ciucci et al., Advancing electrochemical impedance analysis through innovations in the distribution of relaxation times method, Joule. 8 (2024) 1958-1981. https://www.cell.com/joule/fulltext/S2542-4351(24)00236-8.

[3] B. Py, A. Maradesa, F. Ciucci, Gaussian processes for the analysis of electrochemical impedance spectroscopy data: Prediction, filtering, and active learning. Electrochim. Acta. 439 (2023) 141688. https://doi.org/10.1016/j.electacta.2022.141688.

[4] B. Py, C. Zhao, F. Ciucci, Q. Meyer, Gaussian processes for fast and accurate measurements of the polarization resistance of hydrogen fuel cells from impedance spectroscopy, J. Electrochem. Society. 172 (2025) 074502.https://iopscience.iop.org/article/10.1149/1945-7111/ade82c/meta.

[5] B. Py, Z. Wang, Y. Wang, F. Ciucci, Entropy-based regularized regression for advanced distribution of relaxation times deconvolution, J. Power Sources. 644 (2025) 236910. https://doi.org/10.1016/j.jpowsour.2025.236910.

[6] B. Py, F. Ciucci, Beyond ridge regression: Enhancing distribution of relaxation times deconvolution, J. Electrochem. Society. 171 (2024) 060529. https://iopscience.iop.org/article/10.1149/1945-7111/ad576a/meta.

[7] A. Maradesa, B. Py, T.H. Wan, M.B. Effat, F. Ciucci, Selecting the regularization parameter in the distribution of relaxation times, J. Electrochem. Society. 170-3 (2023) 030502. https://doi.org/10.1149/1945-7111/acbca4.

[8] M. Saccoccio, T.H. Wan, C. Chen, F. Ciucci, Optimal regularization in distribution of relaxation times applied to electrochemical impedance spectroscopy: ridge and lasso regression methods-a theoretical and experimental study, Electrochim. Acta. 147 (2014) 470-482. https://doi.org/10.1016/j.electacta.2014.09.058.

[9] B. Py, Z. Wang, Y. Wang, F. Ciucci, Entropy-based regularized regression for advanced distribution of relaxation times deconvolution, J. Power Sources. 644 (2025) 236910. https://doi.org/10.1016/j.jpowsour.2025.236910.
