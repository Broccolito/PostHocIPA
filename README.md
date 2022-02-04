
# PostHocIPA

<!-- badges: start -->
<!-- badges: end -->

The goal of PostHocIPA is to further analyze output files from Ingenuity Pathway Analysis by creating heatmaps and boxplots

## Installation

You can install the development version of PostHocIPA like so:

``` r
library(devtools)
install_github("Broccolito/PostHocIPA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PostHocIPA)
IPA_PostHoc(
  ipa_input_file = "mint_flavor.csv",
  expression_matrix_file = "ecig_count_matrix",
  pathway_file = "mint_flavor_pathway.xls",
  
  sample_list = c(
    "Air_JUUL_1Month_noLPS_Female_S1.results",
    "Air_JUUL_1Month_noLPS_Female_S2.results",
    "Air_JUUL_1Month_noLPS_Female_S3.results",
    "Air_JUUL_1Month_noLPS_Female_S4.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S1.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S2.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S3.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S4.results"
  ),
  
  sample_names = c(
    "Air1",
    "Air2",
    "Air3",
    "Air4",
    "Mint1",
    "Mint2",
    "Mint3",
    "Mint4"
  ),
  
  group_names = c(
    "Air",
    "Air",
    "Air",
    "Air",
    "Mint",
    "Mint",
    "Mint",
    "Mint"
  ),
  
  pathway_header = TRUE)
```

