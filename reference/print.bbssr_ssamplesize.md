# Print Method for bbssr_ssamplesize Objects

Prints the sample size of a stratified test returned by
[`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md).

## Usage

``` r
# S3 method for class 'bbssr_ssamplesize'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_ssamplesize`

- digits:

  Number of decimal places of the power. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
ss <- StratSampleSize(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.5, 0.5),
                      alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                      method = 'standard')
print(ss)
#> Sample size of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Method           : normal approximation
#>   Alpha            : 0.025, target power 0.8
#>   Normal approx.   : 186 (unrounded total)
#>   Total size       : N = 188, power 0.7985
#> 
#>  Stratum frac r  p1  p2 N1 N2
#>        1  0.5 1 0.5 0.3 47 47
#>        2  0.5 1 0.7 0.5 47 47
```
