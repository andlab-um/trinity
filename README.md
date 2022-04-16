# TRINITY

[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2022.04.11.487870-blue)](https://doi.org/10.1101/2022.04.11.487870)
![Twitter URL](https://img.shields.io/twitter/url?label=%40lizhn7&style=social&url=https%3A%2F%2Ftwitter.com%2Flizhn7)

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

## Structure

**This repository contains:**
```
root
 ├── data   # Preprocessed regression data & psychometric data & morphometry data & rs-FC data
 │    ├── dra 
 │    ├── imq 
 │    ├── mms
 │    ├── rsfc
 ├── bs_log # Bootstrapping reuslts
 │    ├── IMQ-MMS
 │    ├── IMQ-rsFC
 │    ├── rsFC-MMS
 ├── util.py            # To provide various basic function
 │── bootstrap.py       # To conduct subject-wise bootstrapping 
 │── eval.py            # To lay the foundation for later analyses
 │── evaluation.ipynb   # To test Hypothesis 1 and Hypothesis 2
 │── dra_data_pre.ipynb # To pave the way for dyadic regression analysis (dra)
 │── dra.ipynb          # To test Hypothesis 3 and plot figure 7
 │── pipeline.ipynb     # To analyse the effects of different pipelines and plot figure 8
 │── figure_2.ipynb     # To plot figure 2
 │── figure_3.ipynb     # To plot figure 3
 │── figure_9.ipynb     # To plot figure 9
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
	     URL = {https://www.biorxiv.org/content/early/2022/04/14/2022.04.11.487870},
	     eprint = {https://www.biorxiv.org/content/early/2022/04/14/2022.04.11.487870.full.pdf},
	     journal = {bioRxiv}
    }
    
___

For bug reports, please contact Zhaoning Li ([yc17319@umac.mo](mailto:yc17319@umac.mo), or [@lizhn7](https://twitter.com/lizhn7)).
___

## LICENSE

This license (CC BY-NC 4.0) gives you the right to re-use and adapt, as long as you note any changes you made, and provide a link to the original source. Read [here](https://creativecommons.org/licenses/by-nc/4.0/) for more details. 

![](https://upload.wikimedia.org/wikipedia/commons/9/99/Cc-by-nc_icon.svg)
