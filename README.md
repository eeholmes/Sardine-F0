# Sardine-F0

Repository for Holmes et al. 2021. Improving landings forecasts using environmental covariates: a case study on the Indian oil sardine (*Sardinella longiceps*). Fisheries Oceanography, in press April 2021.

## Paper

Knit `Manuscript_w_Appendices.Rmd`. This will call the files in the Rmd folder that are needed to build the paper. The paper is formatted to knit to html will look weird in other outputs.

Knit `Supplements.Rmd` to make the Supplements. The Supplements are formatted to be knitted to PDF. They might not look good if knitted to html.

## Data

The raw data are in the SardineForecast package, included as a tar.gz in the data folder. csv of the data are also in the data folder but not used in the scripts

`setup.Rmd` is the key file that takes the monthly data in SardineForecast and prepares the data frames for various monthly averaged covariate data and seasonal catch data used in the paper. This Rmd is called by all the Rmds.

The paper used SardineForecast package verison 1.11 which is included as a tar.gz file in the data folder in the repository. To install from a tar.gz file, use
```
install.packages("data/SardineForecast-1.11-paper.tar.gz", repos = NULL)
```
You will need to be able to build packages from source in order for this to work. If you install packages from GitHub, you already do this.

Once the SardineForecast package is loaded, you can use `?` to get the background information on each data set. Type `?SardineForecast` and then Index at the bottom to find the help files.

If for some reason loading the package does not work, you will find csv files of the data used in the paper in the `data` folder.

## Other files

`AppendixFiles` are Rmds called by `Appendices.Rmd` to make each table. The folder also includes functions that standardize the table building.

`SupplementFiles` has the individual Rmd files and data used to make the tables and figures in the Supplements. Some of the appendix functions in `AppendixFiles` are also used.

