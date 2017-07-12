`sum_of_log_intensity(M, D, ped)`

\[
\ln L(M,D) = \sum_{d\in D}\ln(I(M, d)) - \text{ped.}
\]

\[
\frac{d\ln L(M,D)}{d\lambda} = \sum_{d\in D} \frac{dI(M, d)}{d\lambda} / I(M, d) - \frac{d\text{ped.}}{d\lambda}
\]

`intensity(M, d)` -> `intensity(M.components(), d)`

\[
I(M, d) = \sum_{c\in M} I(c, d)
\]

\[
\frac{dI(M,d)}{d\lambda} = \sum_{c\in M} \frac{dI(c,d)}{d\lambda}
\]

`intensity(c, d)` -> $\beta_c \times$ `intensity(c.decayTrees(), d)`

\[
I(c, d) = \beta_c I(
\]
