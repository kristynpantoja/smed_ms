# smed_ms

Model Selection using sequential minimum energy design of Joseph, Dasgupta, Tuo, Wu

Relevant Files:

-smed_ms_functions.R - has helper functions and main functions for selecting designs to compare (for now, just linear) models given by two hypotheses. Some functions correspond to the One-at-a-Time Algorithm proposed by Joseph et. al. 2015, and some functions correspond to the Fast/Efficient Algorithm in Joseph et. al. 2018.

-bayes_linear_regression.R - testing/visualizing the One-at-a-Time Algorithm functions (sources smed_ms_functions.R)

-efficient_bayes_lm.R - testing/visualizing the Fast/Efficient Algorithm functions (sources smed_ms_functions.R)

-gp_testcode.R - Gaussian Process examples from the book, and first attempts to write Fast/Efficient Algorithm for distinguishing between two gaussian process models.

