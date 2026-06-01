# Beyond Standard Summaries

Supplementary materials for the master's thesis *Beyond Standard Summaries: 
Evaluating and Visualizing the Impact of Information Loss in Clinical Trial 
Adverse Event Reporting*, submitted in partial fulfilment of the requirements 
for the degree of M.Sc. Survey Statistics and Data Analysis at 
Otto-Friedrich-Universität Bamberg, 2026.

## Abstract

Adverse event reporting in clinical trials routinely applies frequency-based 
thresholds and severity aggregation rules that reduce the dimensionality of 
safety data. This thesis develops a quantitative framework to measure the 
resulting information loss using entropy, divergence, coverage, Correspondence 
Analysis and Principal Component Analysis, and evaluates visualization-based 
alternatives through a structured expert study with eight biostatistics 
professionals.

## Repository Contents

| Folder / File | Description |
|---|---|
| `full_pipeline_thesis.R` | Complete analytical pipeline: data loading, cleaning, threshold application, information-loss metrics, multivariate analyses and figure generation |
| `app.R` | Shiny dashboard source code (SafetyLens) |
| `surveys/` | An exemplary survey version as administered via Microsoft Forms (PDF) |
| `data/` | All five snonymous response exports from Microsoft Forms (.xls); no personally identifiable information retained |
| `figures/` | High-resolution PNG exports of all nine figures from Chapters 4 and 5 |

## Reproducibility

All analyses were conducted in R version 4.3.1 (2023-06-16 ucrt) on Windows 10. 
The pipeline is self-contained and requires no external files beyond the 
`pharmaverseadam` package, which provides the synthetic CDISC ADaM dataset used 
throughout. Full package version information is documented in Appendix F of the 
thesis.

## Live Dashboard

The threshold-aware visualization framework is available as an interactive Shiny 
application at:

**https://violaschenk.shinyapps.io/SafetyLens/**

To run locally:
```r
shiny::runApp("app.R")
```

## Citation

Schenk, V. (2026). *Beyond Standard Summaries: Evaluating and Visualizing the 
Impact of Information Loss in Clinical Trial Adverse Event Reporting*. 
Master's thesis, Otto-Friedrich-Universität Bamberg.

## Contact

Viola Schenk — Otto-Friedrich-Universität Bamberg
