# Code for A translucent box; interpretable machine learning in ecology

This is the full code used in the paper "Code for A translucent box; interpretable machine learning in ecology".
The main script is `code/translucent_analysis.R`. 
This script no longer downloads all the data needed and instead you need to get the pantheria data from [here](https://wiley.figshare.com/articles/dataset/Full_Archive/3531875).
The script runs the full analysis and makes all the figures.
It is formatted as an rmarkdown document so running `rmarkdown::render('translucent_analysis.R')` will run the entire analysis and combine the results into a pdf.
However, the full analysis will take a day or more to run on a typical desktop.
`code/translucent_analysis.pdf` is the compiled document.

`code/helpers.R` includes a few small helper functions for plotting outputs from caret models.

In the paper I introduce two new (as far as I know) methods for combining machine learning methods with phylogenetic data.
These might be of particular interest to some.
They start on lines 730 and 814 of the script.
Please see the full paper for details.

The paper is available open access. [https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1422](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1422)
