# Print Method for bbssr_scrp Objects

Prints the design settings, the number of outcomes at which the
conditional rejection probabilities of a stratified design exceed the
nominal level, their largest values and, when the type I error rate has
been decomposed, its largest value over the supplied response
probabilities, as returned by
[`StratCondRejectBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratCondRejectBSSR.md).

## Usage

``` r
# S3 method for class 'bbssr_scrp'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_scrp`

- digits:

  Number of significant digits of the probabilities. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Details

An object whose attributes were removed, for example by selecting some
of its columns, is printed as a data frame.

## Examples

``` r
crp <- StratCondRejectBSSR(Delta.A = 0.3, N1 = c(8, 8), N2 = c(8, 8), omega = 0.5,
                           r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                           ss.method = 'standard')
print(crp)
#> Conditional rejection probabilities of blinded sample size re-estimation
#> for a stratified design, under equal response probabilities in every stratum
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Planned sizes    : (8, 8) and (8, 8)
#>   Interim sizes    : (4, 4) and (4, 4)
#>   Assumed effect   : 0.3 (RD)
#>   Nominal level    : 0.025
#>   Interim outcomes : 81 pairs (s1, s2)
#>   Above the level  : 10605 outcomes for CRP, 17037 for CRP.total
#> 
#> Largest conditional rejection probability
#>   Quantity   value s1 s2 w1 w2
#>        CRP 0.06250  1  1  1  1
#>  CRP.total 0.05853  2  2  0  0
```
