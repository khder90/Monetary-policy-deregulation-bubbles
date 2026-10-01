# Monetary-policy-deregulation-bubbles
# Replication Code

This repository contains the R code used for the empirical analysis in:

"The Impact of Monetary Policy and Central Bank Information Shocks on U.S. Stock Market Bubble Dynamics: The Role of Deregulation"

The script implements:

- Linear local projections
- Interaction local projections
- Smooth varying-coefficient local projections
- Threshold estimation
- Smooth-transition local projections
- Newey-West HAC inference
- Cumulative-response tests
- Block-bootstrap procedures
- Simultaneous confidence bands
- Alternative shock decompositions
- Lag-order robustness
- Outlier robustness
- European comparison exercises

## Main file

`analysis.R`

## Data

The data are not distributed in this repository.

The script expects a workbook named:

`data.xlsx`

with the following worksheets:

- `US_MP_CBI`
- `US_U1U4`
- `EU_Panel`

The underlying data sources are documented in the manuscript.

## Running the analysis

Place `data.xlsx` in the same directory as `analysis.R`, then run:

```r
source("analysis.R")
