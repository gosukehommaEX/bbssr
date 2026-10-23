# Print Method for bbssr_sbssr Objects

Prints the re-estimated group sizes of a stratified trial returned by
[`StratBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratBSSR.md).

## Usage

``` r
# S3 method for class 'bbssr_sbssr'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_sbssr`

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
re <- StratBSSR(n1 = c(15, 15), n2 = c(15, 15), S = c(12, 18), Delta.A = 0.2, r = 1,
                alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(30, 30),
                N2 = c(30, 30), ss.method = 'standard')
print(re)
#> Blinded sample size re-estimation for a stratified trial
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Re-estimation    : normal approximation (unrestricted)
#>   Effect           : RD, assumed 0.2
#>   Final total      : 188, power 0.7985 at the planning probabilities
#> 
#>  Stratum n1 n2  S hat.p hat.p1 hat.p2 N1.final N2.final n1.stage2 n2.stage2
#>        1 15 15 12   0.4    0.5    0.3       47       47        32        32
#>        2 15 15 18   0.6    0.7    0.5       47       47        32        32
```
