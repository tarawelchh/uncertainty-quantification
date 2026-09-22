
# Emulation, History Matching and Optimisation of the Friedman Function

Complex physical models are often too computationally expensive for thorough analysis, sometimes taking days or weeks to run. To address this, Bayes Linear emulators can be used as statistical surrogate functions, as they are much faster to evaluate. 

The following project explores the Friedman function [1]: 
  $f(x) = 10 \sin (\pi x_1 x_2) + 20(x_3-0.5)^2 + 10x_4 + 5x_5$ with $x \in [0,1]^5$.

To replicate the constraints of a simulator, the true function is modelled as being computationally expensive.
Therefore the experimental design is restricted to a minimal number of runs.

## 1D Emulation
* We conduct 1D emulation of the variable $x_3$, fixing all other inputs at 0.5. This variable was chosen to explore
the non-linearity, and also because it has significant effect on the output.
* By performing runs on the boundary points x3 = 0,x3 = 1, we eliminate extrapolation so that emulation is only
performed between runs. This keeps variance bounded.
* The effects of the hyperparameters $\theta$ (correlation length) and %$\sigma^2$ (prior variance) are explored and plotted in figures 2 and 3.
<img width="956" height="745" alt="Screenshot 2026-09-22 at 14 56 41" src="https://github.com/user-attachments/assets/27bf24d9-755a-4514-b85f-1c433f629801" />
## 2D Emulation
* We conduct 2D emulation of $x_1$ and $x_2$ as they have an interaction term.
* A Latin Hypercube Design is used to allocate 20 well-spaced input points.

## History Matching 
* We perform 2D history matching on the pairs ($x_1$,$x_2$) and then to ($x_1$,$x_3$), each with two waves, again with other inputs
set to 0.5 and $z=14.5$.
* Emulator adjusted expectation and variance are plotted, as well as plots of non-implausible regions.
* These are compared with the true $f(x)$, though this would not be possible for an unknown function.


<img width="836" height="594" alt="Screenshot 2026-09-22 at 14 56 06" src="https://github.com/user-attachments/assets/0ddc350e-f900-42b9-ac43-44c5234b7b82" />
 <img width="621" height="505" alt="Screenshot 2026-09-22 at 14 55 46" src="https://github.com/user-attachments/assets/a61b27a9-a34c-4411-a325-f30704a9d151" />

  ## Optimisation
* We can approach optimisation with a history matching style strategy to find a set of inputs that maximise our function
$f(x)$. For this, to obtain interesting plots, we again focus on the variables $(x_1,x_2)$, since the quadratic and linear terms
are simply maximised at their boundaries. To investigate the global max, we set $x_3 = x_4 = x_5 = 1$.
<img width="957" height="378" alt="Screenshot 2026-09-22 at 14 54 50" src="https://github.com/user-attachments/assets/91d99d8b-bd4d-44bc-a4a5-691f2fcfcced" />
<img width="794" height="316" alt="Screenshot 2026-09-22 at 14 55 28" src="https://github.com/user-attachments/assets/2237f596-ffa1-497e-91ec-49de13e0d4ac" />


## Overview
This project explores the use of simple Bayes Linear emulators on the non-linear Friedman function.
* Tuning of sigma and theta was performed, highlighting its influence on ensuring the emulator expectation follows
the shape of the true function, and to obtain sufficiently narrow prediction intervals without making overly confident
predictions.
* History matching of additive terms ($x_1$ vs $x_3$) may only require a single wave for accurate results, whereas for
interaction terms ($x_1$ vs $x_2$), 2 waves were required to emulate the region successfully.
* In all cases of history matching, the emulator’s non-implausible region successfully encapsulated the true model
contours. This is important as it ensures there are no false negatives.
* It is important to acknowledge that the Friedman function is smooth, and so emulator behaviour may represent a
best case scenario. A more rigorous stress-test could include non-smooth, discontinuous functions, possibly with
multiple interaction terms.

Overall this report demonstrated techniques for exploring complex, non-linear models, using an emulator as a robust,
low-cost surrogate. The BL emulators were able to characterise surfaces of parameter interaction as well as identify global
optima, using iterative design to ensure low computational cost.


  [1] J.H. Friedman. Multivariate adaptive regression splines. The Annals of Statistics, 1991.
