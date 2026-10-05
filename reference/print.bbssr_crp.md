# Print Method for bbssr_crp Objects

Prints the design settings, the number of outcomes at which the
conditional rejection probabilities exceed the nominal level, their
largest values and, when the type I error rate has been decomposed, its
largest value over the supplied response probabilities.

## Usage

``` r
# S3 method for class 'bbssr_crp'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_crp`

- digits:

  Number of significant digits of the probabilities. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
crp <- BinaryCondRejectBSSR(
  Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Fisher', ss.method = 'standard',
  theta = c(0.3, 0.5)
)
print(crp)
#> Conditional rejection probabilities of blinded sample size re-estimation
#> for a binary endpoint, under equal response probabilities
#> 
#>   Test            : Fisher
#>   Alternative     : greater
#>   Initial size    : N1 = 12, N2 = 12
#>   Interim size    : n1 = 6, n2 = 6
#>   Assumed effect  : 0.3 (RD)
#>   Nominal level   : 0.025
#>   Outcomes        : 567 pairs (s, s2)
#>   Above the level : 8 for CRP, 0 for CRP.total
#> 
#> Largest conditional rejection probability
#>   Quantity   value s s2 N1 N2
#>        CRP 0.02597 2  3 24 24
#>  CRP.total 0.02500 3 15 32 32
#> 
#> Type I error rate over 2 values of theta: largest 0.01581 at theta = 0.5
```
