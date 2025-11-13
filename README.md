# minuit1

- expFit.cpp: example using the more general fitting interface in minuit
- expFit.ipynb: equivalent example using lmfit
- SimultaneousExps(lm).ipynb: generation of histograms with correlated signals for simultaneous fit exercise
- rootExample.cpp: just another example of using ROOT classes in a C++ program
- *.root files: various input histograms for fitting exercises

-----

Team member names and computing IDs: hxd5bc (just me)

I made the three pdfs with comments for the exercises in overleaf but it makes the plots small. All of the plots are in my repo so you can view them outside of the pdfs. I also copied my comments into the readme (hopefully i got all of them)

I added two cpp scripts for the first exercise. The first is expFit2.cpp which is a root macro to do the chi2 stuff. the example script we were given is run with a makefile and i don't like that as much. The result is just test.png because we didn't necessarily need to show the plot. The second script is FitDistros.cpp which did the double gauss and gumbel fits. Besides ex1.pdf, the plot can be found in fitDistros_results.png
-----

Exercise 2 comments:
--
from the pdf that we were supposed to make (maybe idk)

The output from my script is:

Mean  = 74.6126 ± 0.452681
Sigma = 4.53804 ± 0.416155
chi2/ndf = 76.3408 / 93 = 0.820869
Chi2 probability = 0.894986

The generated signal has mean 75 and sigma 4.5 so this fit does a very good job extracting that information. The reduced chi2 is below 1 so errors may have been overestimated. The fit looks good to me because of the high probability and close-to-actual mean and sigma. 

The plot can also be viewed as simulFit_results.png

-----

Exercise 3 comments:
--

here's the script output: 

=== Fit Parameter Results ===
A = 33575.6 ± 275.146
mu1 = 3.55834 ± 0.00587431
mu2 = 1.56372 ± 0.0207186
sigma1 = 1 ± 7.70741e-05
sigma2 = 1.589 ± 0.0170731
B = 19736.5 ± 248.462

Chi2 = 6765.34   NDF = 3594   Reduced Chi2 = 1.8824

A should be the number of signals so there were 33575.6 ± 275.146 signal events. B covers the background which has 19736.5 ± 248.462 events. Visually my plots look right. The background is very heavy in the one corner, and my script successfully removes background events from that corner. The model being a 2D gaussian with that wavy background is visually what my fit result looks like, and the reduced chi2 is 1.88, meaning the model does well to describe the data. 

The plot is also found as fit2D_results.png
-----
