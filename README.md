# TRINITY <img src="https://raw.githubusercontent.com/andlab-um/trinity/main/demo.png" align="right" width="561px">

[![GitHub repo size](https://img.shields.io/github/languages/code-size/andlab-um/trinity?color=brightgreen&label=repo%20size&logo=github)](https://github.com/andlab-um/trinity)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2022.04.11.487870-blue)](https://doi.org/10.1101/2022.04.11.487870)<br />
[![Twitter URL](https://img.shields.io/twitter/url?label=%40lizhn7&style=social&url=https%3A%2F%2Ftwitter.com%2Flizhn7)](https://twitter.com/lizhn7)
[![Twitter URL](https://img.shields.io/twitter/url?label=%40ANDlab3&style=social&url=https%3A%2F%2Ftwitter.com%2Flizhn7)](https://twitter.com/ANDlab3)

**Code and data for: <br />**
**Li, Z., Dong, Q., Hu, B., & Wu, H. (2022). Every individual makes a difference: A trinity derived from linking individual brain morphometry, connectivity and mentalising ability.** *bioRxiv*. <br />
[DOI: 10.1101/2022.04.11.487870](https://doi.org/10.1101/2022.04.11.487870).
___

## Highlights
- A trinity existed in idiosyncratic patterns of brain MMS, rs-FC and IMQ scores.
- A region-related mentalising specificity emerged from this trinity.
- Rs-FC gates the MMS predicted similarity in mentalising ability.
- IS-RSA can be used to link brain morphometry, connectivity and mentalising ability.
- CPP-SD robustly improved the relatedness between MMS and the other two modalities.
___

## Outreach

- A poster for the Social & Affective Neuroscience Society (SANS) Annual Meeting 2022 is available on [ResearchGate](https://www.researchgate.net/publication/360262009_Every_individual_makes_a_difference_A_trinity_derived_from_linking_individual_brain_morphometry_connectivity_and_mentalising_ability).
- A 4.6-min video for the SANS 2022 is available on [YouTube](https://youtu.be/kmTiUy0SowA). The related slides are available on [ResearchGate](https://www.researchgate.net/publication/360262895_Every_individual_makes_a_difference_A_trinity_derived_from_linking_individual_brain_morphometry_connectivity_and_mentalising_ability).
- [IMQ: Interactive mentalising questionnaire](https://github.com/andlab-um/IMQ).
- [Previous resting-network study](https://github.com/andlab-um/restDishonesty).
___

## Structure

**This repository contains:**
```
root
 ????????? data               # Preprocessed regression data & psychometric data & morphometry data & rs-FC data
 ???    ????????? dra 
 ???    ????????? imq 
 ???    ????????? mms
 ???    ????????? rsfc
 ????????? bs_log             # Bootstrapping reuslts
 ???    ????????? IMQ-MMS
 ???    ????????? IMQ-rsFC
 ???    ????????? rsFC-MMS
 ????????? util.py            # To provide various basic function
 ????????? bootstrap.py       # To conduct subject-wise bootstrapping 
 ????????? eval.py            # To lay the foundation for later analyses
 ????????? evaluation.ipynb   # To test Hypothesis 1 and Hypothesis 2
 ????????? dra_data_pre.ipynb # To pave the way for dyadic regression analysis (dra)
 ????????? dra.ipynb          # To test Hypothesis 3 and plot figure 7
 ????????? pipeline.ipynb     # To analyse the effects of different pipelines and plot figure 8
 ????????? figure_2.ipynb     # To plot figure 2
 ????????? figure_3.ipynb     # To plot figure 3
 ????????? figure_9.ipynb     # To plot figure 9
 ????????? LICENSE
 ????????? README.md
```
**Note 1**: to properly run all scripts, you may need to set the root of this repository as your working directory. <br />
**Note 2**: to properly run dyadic regression analysis, you may need to install [Pymer4](https://eshinjolly.com/pymer4/). <br />
**Note 3**: to reproduce figure 3, you may need version 3.7 of Python. <br />
___

## Citation

    @article {Li2022.04.11.487870,
      author = {Li, Zhaoning and Dong, Qunxi and Hu, Bin and Wu, Haiyan},
      title = {Every individual makes a difference: A trinity derived from linking individual brain morphometry, connectivity and mentalising ability},
      elocation-id = {2022.04.11.487870},
      year = {2022},
      doi = {10.1101/2022.04.11.487870},
      publisher = {Cold Spring Harbor Laboratory},
      URL = {https://www.biorxiv.org/content/early/2022/04/24/2022.04.11.487870},
      eprint = {https://www.biorxiv.org/content/early/2022/04/24/2022.04.11.487870.full.pdf},
      journal = {bioRxiv}
    }
    
___

For bug reports, please contact Zhaoning Li ([yc17319@umac.mo](mailto:yc17319@umac.mo), or [@lizhn7](https://twitter.com/lizhn7)).

Thanks to [shields.io](https://shields.io/) and [Dr. Lei Zhang](https://github.com/lei-zhang) for the [good example](https://github.com/lei-zhang/SIT).
___

## LICENSE

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>, which gives you the right to re-use and adapt, as long as you note any changes you made, and provide a link to the original source.
