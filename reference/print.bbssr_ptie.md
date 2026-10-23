# Print Method for bbssr_ptie Objects

Prints the settings and the largest type I error rates of a
matched-pairs design with blinded sample size re-estimation and of the
fixed design.

## Usage

``` r
# S3 method for class 'bbssr_ptie'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_ptie` returned by
  [`PairedTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedTypeIErrorBSSR.md)

- digits:

  Number of significant digits of the rates. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
tie <- PairedTypeIErrorBSSR(Delta.A = 0.2, N = 40, n.interim = 20, alpha = 0.05,
                            tar.power = 0.8, Test = 'McNemar',
                            alternative = 'two.sided', ss.method = 'standard',
                            psi = seq(0.1, 0.9, by = 0.1))
print(tie)
#> Type I error rate of blinded sample size re-estimation for paired binary data
#> 
#>   Test            : McNemar
#>   Alternative     : two.sided
#>   Pairs           : N = 40 planned, 20 at the interim analysis
#>   Assumed effect  : 0.2
#>   Nominal level   : 0.05
#>   Grid            : 9 values of psi in [0.1, 0.9]
#>   Maximum         : certified over psi in [0.1, 0.9]
#> 
#> Largest type I error rate
#>        Design      psi     TIE
#>          BSSR 0.388446 0.05058
#>  Fixed sample 0.600769 0.05383
```
