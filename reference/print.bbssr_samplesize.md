# Print a Sample Size Calculation

Print a Sample Size Calculation

## Usage

``` r
# S3 method for class 'bbssr_samplesize'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_samplesize` returned by
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)

- digits:

  Number of significant digits used for the power

- ...:

  Further arguments, currently ignored

## Value

The object `x`, returned invisibly

## Examples

``` r
ss <- BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
                       tar.power = 0.8, Test = 'Chisq')
print(ss)
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Chisq
#>   Alternative      : greater
#>   Response rates   : p1 = 0.4, p2 = 0.2
#>   Allocation ratio : 1 to 1
#>   Alpha            : 0.025
#>   Target power     : 0.8
#> 
#>   Required sample size: N1 = 80, N2 = 80, total N = 160
#>   Attained power      : 0.8009
```
