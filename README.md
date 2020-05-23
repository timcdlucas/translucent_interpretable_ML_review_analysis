# Code for A translucent box; interpretable machine learning in ecology

This is the full code used in the paper "Code for A translucent box; interpretable machine learning in ecology".
The main script is `code/translucent_analysis.R`. 
This script downloads all the data needed and performs the full analysis and makes all the figures.
It is formatted as an rmarkdown document so running `rmarkdown::render('translucent_analysis.R')` will run the entire analysis and combine the results into a pdf.
However, the full analysis will take a day or more to run on a typical desktop.
`code/translucent_analysis.pdf` is the compiled document.

`code/helpers.R` includes a few small helper functions for plotting outputs from caret models.

