---
title: 'Unfold.jl: Regression models for Cognitive Neuroscience'
tags:
  - Julia
  - EEG
  - neuroimaging
authors:
  - name: Benedikt V. Ehinger
    orcid:  0000-0002-6276-3332 
    equal-contrib: false
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  
affiliations:
 - name: Stuttgart Center for Simulation Science, University of Stuttgart, Germany
   index: 1
 - name: Institute for Visualization and Interactive Systems, University of Stuttgart, Germany
   index: 2
date: 06 November 2023
bibliography: paper.bib
---

# Summary


Overlapping event responses, (non-)linear effects and confounds, outliers, item effects - all pretty common features limiting the application of neuroimaging methods. Here, we focus on EEG without loss of generality. Especially event-related potentials in naturalistic paradigms (e.g. ET+EEG, VR+EEG, MoBi) but also classical paradigms (evidence Integration, language, even oddball tasks) require such techniques. 

Unfold.jl implements the model fitting based on the popular wilkinson~formula+syntax. It is feature-par with matlab-unfold. In addition it offers speed-improvements (up to 100x on GPUs) and extensions for outlier-robust fits, back2back regression, and mixed models.

While implemented in Julia, it is straight forward to also call from Python and tutorials are available.



# Statement of need



# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
