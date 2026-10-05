# Print a Rejection Region

Print a Rejection Region

## Usage

``` r
# S3 method for class 'bbssr_rr'
print(x, show.map = NULL, ...)
```

## Arguments

- x:

  An object of class `bbssr_rr` returned by
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)

- show.map:

  Logical. If `TRUE`, the outcome grid is printed as a map in which `X`
  marks rejection of the null hypothesis. The default prints the map
  only when the grid has at most 400 cells

- ...:

  Further arguments, currently ignored

## Value

The object `x`, returned invisibly

## Examples

``` r
RR <- BinaryRR(N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
print(RR)
#> Rejection region for a two-arm trial with a binary endpoint
#> 
#>   Test          : Chisq
#>   Alternative   : greater
#>   Sample sizes  : N1 = 5, N2 = 5
#>   Alpha         : 0.025
#> 
#>   Rejected outcomes: 5 of 36
#> 
#>      x2=0 x2=1 x2=2 x2=3 x2=4 x2=5
#> x1=0 .    .    .    .    .    .   
#> x1=1 .    .    .    .    .    .   
#> x1=2 .    .    .    .    .    .   
#> x1=3 X    .    .    .    .    .   
#> x1=4 X    .    .    .    .    .   
#> x1=5 X    X    X    .    .    .   
#> 
#>   X marks rejection of the null hypothesis
```
