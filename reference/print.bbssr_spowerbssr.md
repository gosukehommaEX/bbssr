# Print Method for bbssr_spowerbssr Objects

Prints the power and the expected final total sample size of a
stratified design with blinded sample size re-estimation returned by
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md).

## Usage

``` r
# S3 method for class 'bbssr_spowerbssr'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_spowerbssr`

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
res <- StratPowerBSSR(p = c(0.3, 0.5), Delta.A = 0.2, Delta.T = 0.2,
                      N1 = c(30, 30), N2 = c(30, 30), omega = 0.5, r = 1,
                      alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                      ss.method = 'standard')
print(res)
#> Blinded sample size re-estimation for a stratified design
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Re-estimation    : normal approximation (unrestricted)
#>   Planned sizes    : (30, 30) and (30, 30)
#>   Interim sizes    : (15, 15) and (15, 15)
#>   Effect           : RD, assumed 0.2, true 0.2
#>   Alpha            : 0.025, target power 0.8
#> 
#>  p.S1 p.S2 p1.S1 p1.S2 p2.S1 p2.S2 power.BSSR power.TRAD    E.N
#>   0.3  0.5   0.4   0.6   0.2   0.4     0.7892     0.6274 174.55
```
