# Print Method for bbssr_srr Objects

Prints the critical values of a stratified test returned by
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md).

## Usage

``` r
# S3 method for class 'bbssr_srr'
print(x, n = 20, ...)
```

## Arguments

- x:

  An object of class `bbssr_srr`

- n:

  Number of rows of the table to print. Default is 20

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Details

An object whose attributes were removed, for example by selecting some
of its columns, is printed as a data frame.

## Examples

``` r
rr <- StratRR(N1 = c(10, 12), N2 = c(10, 12), alpha = 0.025, Test = 'CMH')
print(rr, n = 10)
#> Rejection region of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Group sizes      : N1 = (10, 12), N2 = (10, 12)
#>   Alpha            : 0.025
#>   Combinations     : 525, of which 485 have a rejected value of S
#> 
#>  z.S1 z.S2 S.min S.max lower upper
#>     4    0     0     4    NA     4
#>     5    0     0     5    NA     5
#>     6    0     0     6    NA     6
#>     7    0     0     7    NA     6
#>     8    0     0     8    NA     7
#>     9    0     0     9    NA     7
#>    10    0     0    10    NA     8
#>    11    0     1    10    NA     8
#>    12    0     2    10    NA     9
#>    13    0     3    10    NA     9
#>   ... 475 more rows
```
