# Summary Method for bbssr_grid Objects

Summarizes each design of
[`BinaryGridBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryGridBSSR.md)
over the true pooled response probabilities at which it was evaluated.

## Usage

``` r
# S3 method for class 'bbssr_grid'
summary(object, ...)
```

## Arguments

- object:

  An object of class `bbssr_grid`

- ...:

  Further arguments, ignored

## Value

A data frame with one row per design, containing the attribute `designs`
of `object` followed by the smallest and the largest power of the BSSR
design (`power.BSSR.min`, `power.BSSR.max`) and of the fixed-sample
design (`power.TRAD.min`, `power.TRAD.max`), and the largest expected
total sample size `E.N.max`

## Examples

``` r
design <- data.frame(Test = c('Chisq', 'Fisher'))
grid <- BinaryGridBSSR(design, p = c(0.3, 0.4, 0.5), Delta.A = 0.3, N1 = 30, N2 = 30,
                       omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                       ss.method = 'standard')
summary(grid)
#>   design   Test N1 N2 n1.interim n2.interim Delta.T power.BSSR.min
#> 1      1  Chisq 30 30         15         15     0.3      0.7885735
#> 2      2 Fisher 30 30         15         15     0.3      0.7026320
#>   power.BSSR.max power.TRAD.min power.TRAD.max  E.N.max
#> 1      0.7982702      0.6617280      0.7369757 83.27159
#> 2      0.7265063      0.5593173      0.6447901 83.27159
```
