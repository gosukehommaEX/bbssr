# Print Method for bbssr_pbssr Objects

Prints the result of a blinded sample size re-estimation from paired
interim data.

## Usage

``` r
# S3 method for class 'bbssr_pbssr'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_pbssr` returned by
  [`PairedBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedBSSR.md)

- digits:

  Number of decimal places. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Details

An object whose attributes were removed, for example by selecting some
of its rows or columns, is printed as a data frame.

## Examples

``` r
res <- PairedBSSR(n1 = 25, D1 = 5, Delta.A = 0.15, alpha = 0.05, tar.power = 0.8,
                  Test = 'McNemar', alternative = 'two.sided', ss.method = 'standard')
print(res)
#> Blinded sample size re-estimation from paired interim data
#> 
#>   Test             : McNemar
#>   Alternative      : two.sided
#>   Re-estimation    : normal approximation (unrestricted)
#>   Assumed effect   : 0.15
#>   Alpha            : 0.05, target power 0.8
#> 
#> Interim data
#>   Pairs            : n1 = 25, discordant D1 = 5
#>   Discordance      : estimate 0.2, used 0.2
#>   Planning values  : p10 = 0.175, p01 = 0.025
#> 
#> Re-estimation
#>   Re-estimated N   : 68
#>   Still to enrol   : 43 pairs
#>   Final N          : 68
#>   Power at final N : 0.8427
#>   Conditional power: 0.8519
```
