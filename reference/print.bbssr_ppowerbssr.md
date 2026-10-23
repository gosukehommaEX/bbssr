# Print Method for bbssr_ppowerbssr Objects

Prints the power and the expected final number of pairs of a
matched-pairs design with blinded sample size re-estimation.

## Usage

``` r
# S3 method for class 'bbssr_ppowerbssr'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_ppowerbssr` returned by
  [`PairedPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedPowerBSSR.md)

- digits:

  Number of decimal places. Default is 4

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Examples

``` r
res <- PairedPowerBSSR(psi = c(0.2, 0.3), Delta.T = 0.15, Delta.A = 0.15, N = 72,
                       n.interim = 30, alpha = 0.05, tar.power = 0.8,
                       Test = 'McNemar', alternative = 'two.sided',
                       ss.method = 'standard')
print(res)
#> Blinded sample size re-estimation for paired binary data
#> 
#>   Test             : McNemar
#>   Alternative      : two.sided
#>   Re-estimation    : normal approximation (unrestricted)
#>   Pairs            : N = 72 planned, 30 at the interim analysis
#>   Effect           : assumed 0.15, true 0.15
#>   Alpha            : 0.05, target power 0.8
#> 
#>  psi   p10   p01 power.BSSR power.TRAD    E.N
#>  0.2 0.175 0.025     0.7986     0.8654  71.22
#>  0.3 0.225 0.075     0.7799     0.6773 102.98
```
