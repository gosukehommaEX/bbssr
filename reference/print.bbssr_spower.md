# Print Method for bbssr_spower Objects

Prints the exact power of a stratified test returned by
[`StratPower`](https://gosukehommaex.github.io/bbssr/reference/StratPower.md).

## Usage

``` r
# S3 method for class 'bbssr_spower'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_spower`

- digits:

  Number of decimal places. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Details

An object whose attributes were removed, for example by selecting some
of its columns, is printed as a data frame.

## Examples

``` r
pw <- StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), N1 = c(40, 40), N2 = c(40, 40),
                 alpha = 0.025, Test = 'CMH')
print(pw)
#> Exact power of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Group sizes      : N1 = (40, 40), N2 = (40, 40)
#>   Alpha            : 0.025
#> 
#>  p1.S1 p1.S2 p2.S1 p2.S2  Power
#>    0.5   0.7   0.3   0.5 0.7322
```
