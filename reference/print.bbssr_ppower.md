# Print Method for bbssr_ppower Objects

Prints the exact power of a test for paired binary data.

## Usage

``` r
# S3 method for class 'bbssr_ppower'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_ppower` returned by
  [`PairedPower`](https://gosukehommaex.github.io/bbssr/reference/PairedPower.md)

- digits:

  Number of decimal places of the power. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
pw <- PairedPower(p10 = c(0.3, 0.35), p01 = 0.1, N = 40, alpha = 0.025,
                  Test = 'McNemar-exact')
print(pw)
#> Exact power for paired binary data
#> 
#>   Test        : McNemar-exact
#>   Alternative : greater
#>   Pairs       : N = 40
#>   Alpha       : 0.025
#> 
#>   p10 p01  Power
#>  0.30 0.1 0.4397
#>  0.35 0.1 0.6052
```
