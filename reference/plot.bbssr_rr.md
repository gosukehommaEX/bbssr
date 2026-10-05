# Plot a Rejection Region

Displays the outcome grid of a two-arm trial with a binary endpoint and
shades the outcomes for which the null hypothesis is rejected.

## Usage

``` r
# S3 method for class 'bbssr_rr'
plot(
  x,
  main = NULL,
  sub = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.title = NULL,
  legend.labels = NULL,
  colours = NULL,
  base_size = 11,
  ...
)
```

## Arguments

- x:

  An object of class `bbssr_rr` returned by
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)

- main:

  Title of the plot. The default of `NULL` builds the title from the
  object and a single `NA` drops it

- sub:

  Subtitle of the plot, following the convention of `main`

- xlab:

  Label of the horizontal axis, following the convention of `main`

- ylab:

  Label of the vertical axis, following the convention of `main`

- legend.title:

  Title of the legend. The default of `NULL` leaves the legend untitled

- legend.labels:

  Character vector of length two replacing the entries of the legend, in
  the order retained then rejected

- colours:

  Character vector of length two giving the fill of the tiles, in the
  order retained then rejected

- base_size:

  Base font size of the theme, in points

- ...:

  Further arguments, currently ignored

## Value

A `ggplot` object

## Details

Both axes count responders, so the breaks are restricted to integers.
The vertical axis runs downwards, which places the outcome with no
responders in the top left corner and matches the layout of the map
printed by
[`print.bbssr_rr`](https://gosukehommaex.github.io/bbssr/reference/print.bbssr_rr.md).

## Examples

``` r
RR <- BinaryRR(N1 = 10, N2 = 10, alpha = 0.025, Test = 'Chisq')
plot(RR)


# Enlarge the type and relabel the legend
plot(RR, base_size = 14, legend.title = 'Decision',
     legend.labels = c('do not reject', 'reject'),
     colours = c('white', 'grey30'))

```
