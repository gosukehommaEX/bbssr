# Print Method for bbssr_stie Objects

Prints the largest type I error rates of a stratified design with
blinded sample size re-estimation and of the fixed design, returned by
[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md).

## Usage

``` r
# S3 method for class 'bbssr_stie'
print(x, digits = 6, ...)
```

## Arguments

- x:

  An object of class `bbssr_stie`

- digits:

  Number of significant digits. Default is 6

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Details

An object whose attributes were removed, for example by selecting some
of its columns, is printed as a data frame.

## Examples

``` r
tie <- StratTypeIErrorBSSR(Delta.A = 0.2, N1 = c(20, 20), N2 = c(20, 20), omega = 0.5,
                           r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                           ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.2))
print(tie)
#> Type I error rate of blinded sample size re-estimation for a stratified design
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Planned sizes    : (20, 20) and (20, 20)
#>   Points           : 25
#>   Maximum          : certified over (theta.S1, theta.S2) in [0.1, 0.9] x [0.1, 0.9]
#>   Alpha            : 0.025
#> 
#>        Design theta.S1 theta.S2       TIE     bound
#>          BSSR      0.5      0.5 0.0250852 0.0250852
#>  Fixed sample      0.5      0.5 0.0283368 0.0283368
```
