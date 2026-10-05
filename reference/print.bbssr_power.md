# Print Exact Power Results

Print Exact Power Results

## Usage

``` r
# S3 method for class 'bbssr_power'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_power` returned by
  [`BinaryPower`](https://gosukehommaex.github.io/bbssr/reference/BinaryPower.md)

- digits:

  Number of significant digits used for the power

- ...:

  Further arguments, currently ignored

## Value

The object `x`, returned invisibly

## Examples

``` r
pw <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
print(pw)
#> Exact power for a two-arm trial with a binary endpoint
#> 
#>   Test         : Chisq
#>   Alternative  : greater
#>   Sample sizes : N1 = 5, N2 = 5
#>   Alpha        : 0.025
#> 
#>   p1  p2 Power
#>  0.5 0.2 0.183
```
