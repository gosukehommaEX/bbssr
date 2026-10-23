# Print Method for bbssr_psamplesize Objects

Prints the number of pairs required for a test for paired binary data.

## Usage

``` r
# S3 method for class 'bbssr_psamplesize'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_psamplesize` returned by
  [`PairedSampleSize`](https://gosukehommaex.github.io/bbssr/reference/PairedSampleSize.md)

- digits:

  Number of decimal places of the power. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
ss <- PairedSampleSize(p10 = 0.2, p01 = 0.1, alpha = 0.05, tar.power = 0.9,
                       Test = 'McNemar', method = 'standard')
print(ss)
#> Number of pairs for paired binary data
#> 
#>   Test           : McNemar
#>   Alternative    : greater
#>   Discordance    : p10 = 0.2, p01 = 0.1
#>   Alpha          : 0.05
#>   Target power   : 0.9
#>   Method         : normal approximation
#> 
#>   Required number of pairs: N = 254
#>   Attained power          : 0.9062
```
