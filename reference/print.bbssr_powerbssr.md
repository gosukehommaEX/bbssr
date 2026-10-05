# Print the Operating Characteristics of a BSSR Design

Print the Operating Characteristics of a BSSR Design

## Usage

``` r
# S3 method for class 'bbssr_powerbssr'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_powerbssr` returned by
  [`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)

- digits:

  Number of significant digits used for the power

- ...:

  Further arguments, currently ignored

## Value

The object `x`, returned invisibly

## Examples

``` r
res <- BinaryPowerBSSR(
  p = 0.45,
  Delta.A = 0.3, Delta.T = 0.3,
  N1 = 5, N2 = 5, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
)
print(res)
#> Blinded sample size re-estimation for a binary endpoint
#> 
#>   Test            : Chisq
#>   Alternative     : greater
#>   Design rule     : unrestricted
#>   Initial size    : N1 = 5, N2 = 5
#>   Interim fraction: 0.5, giving n1 = 3 and n2 = 3
#>   Treatment effect: assumed 0.3, true 0.3 (risk difference)
#>   Re-estimation   : exact power of Chisq at level 0.025
#>   Alpha           : 0.025, target power 0.8
#> 
#>     p  p1  p2 power.BSSR power.TRAD  E.N
#>  0.45 0.6 0.3      0.719     0.1667 73.3
```
