# Print Method for bbssr_tie Objects

Prints the design settings and the largest type I error rate of a design
with blinded sample size re-estimation and of the corresponding
fixed-sample design.

## Usage

``` r
# S3 method for class 'bbssr_tie'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_tie`

- digits:

  Number of significant digits of the type I error rates. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
tie <- BinaryTypeIErrorBSSR(
  Delta.A = 0.3, N1 = 20, N2 = 20, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
  theta = seq(0.05, 0.95, by = 0.05)
)
print(tie)
#> Type I error rate of blinded sample size re-estimation for a binary endpoint
#> 
#>   Test            : Chisq
#>   Alternative     : greater
#>   Initial size    : N1 = 20, N2 = 20
#>   Interim size    : n1 = 10, n2 = 10
#>   Assumed effect  : 0.3 (RD)
#>   Nominal level   : 0.025
#>   Grid            : 19 values of theta in [0.05, 0.95], maxima refined
#> 
#> Largest type I error rate
#>        Design    theta     TIE
#>          BSSR 0.500000 0.02808
#>  Fixed sample 0.312886 0.02669
```
