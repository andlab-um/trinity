# TRINITY
**Every individual makes a difference:** A trinity derived from linking individual brain morphometry, connectivity and mentalising ability

___
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

### LICENSE

This license (CC BY-NC 4.0) gives you the right to re-use and adapt, as long as you note any changes you made, and provide a link to the original source. Read [here](https://creativecommons.org/licenses/by-nc/4.0/) for more details. 

![](https://upload.wikimedia.org/wikipedia/commons/9/99/Cc-by-nc_icon.svg)
