# Print Method for bbssr_alphaadj Objects

Prints the adjusted nominal significance levels of a design with blinded
sample size re-estimation and of the corresponding fixed-sample design.

## Usage

``` r
# S3 method for class 'bbssr_alphaadj'
print(x, digits = 6, ...)
```

## Arguments

- x:

  An object of class `bbssr_alphaadj`

- digits:

  Number of significant digits. Default is 6

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
# \donttest{
adj <- BinaryAlphaAdjBSSR(
  Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
  theta = seq(0.01, 0.99, by = 0.01)
)
print(adj)
#> Adjusted significance level for blinded sample size re-estimation
#> 
#>   Test            : Chisq
#>   Alternative     : greater
#>   Adjusted part   : final analysis only
#>   Maximum         : certified over theta in [0.01, 0.99]
#>   Target level    : 0.025
#> 
#>        Design   max.TIE alpha.adj max.TIE.adj
#>          BSSR 0.0272864 0.0204376   0.0237418
#>  Fixed sample 0.0293761 0.0207700   0.0247660
# }
```
