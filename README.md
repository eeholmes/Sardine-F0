# Sardine-F0

## Data

The data for the paper are in the SardineForecast package verison 1.11 which is included as a tar.gz file in the data folder in the repository. To install from a tar.gz file, use
```
install.packages("data/SardineForecast-1.11-paper.tar.gz", repos = NULL)
```
You will need to be able to build packages from source in order for this to work. If you install packages from GitHub, you already do this.

Once the SardineForecast package is loaded, you can use `?` to get the background information on each data set. Type `?SardineForecast` and then Index at the bottom to find the help files.

If for some reason loading the package does not work, you will find csv files of the data used in the paper in the `data` folder.

## Paper

Knit `Manuscript_w_Appendices.Rmd`. This will call the files in the Rmd folder that are needed to build the paper.

Knit `Supplements.Rmd` to make the Supplements.

## Other files

`AppendixFiles` are Rmds called by `Appendices.Rmd` to make each table. The folder also includes functions that standardize the table building.

`SupplementFiles` has the individual Rmd files and data used to make the tables and figures in the Supplements. Some of the appendix functions in `AppendixFiles` are also used.
