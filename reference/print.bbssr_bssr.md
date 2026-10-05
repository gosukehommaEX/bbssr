# Print a Sample Size Re-estimation from Interim Data

Print a Sample Size Re-estimation from Interim Data

## Usage

``` r
# S3 method for class 'bbssr_bssr'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bbssr_bssr` returned by
  [`BinaryBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md)

- digits:

  Number of significant digits used for the proportions and the power

- ...:

  Further arguments, currently ignored

## Value

The object `x`, returned invisibly

## Examples

``` r
res <- BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
                  alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
print(res)
#> Blinded sample size re-estimation from interim data
#> 
#>   Test             : Chisq
#>   Alternative      : greater
#>   Design rule      : unrestricted
#>   Assumed effect   : 0.3
#>   Alpha            : 0.025, target power 0.8
#> 
#> Interim data
#>   Patients         : n1 = 20, n2 = 20, total n = 40
#>   Responders       : S = 11, blinded pooled rate = 0.275
#>   Recovered rates  : hat.p1 = 0.425, hat.p2 = 0.125
#> 
#> Re-estimation
#>   Required total   : N1 = 33, N2 = 33, total N = 66
#>   Still to enrol   : group 1 = 13, group 2 = 13, total = 26
#>   Final size       : N1 = 33, N2 = 33, total N = 66
#>   Power at final N : 0.8136
```
