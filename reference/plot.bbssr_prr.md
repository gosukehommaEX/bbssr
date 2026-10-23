# Plot a Rejection Region for Paired Binary Data

Displays the possible outcomes of a matched-pairs trial and shades the
outcomes for which the null hypothesis is rejected.

## Usage

``` r
# S3 method for class 'bbssr_prr'
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

  An object of class `bbssr_prr` returned by
  [`PairedRR`](https://gosukehommaex.github.io/bbssr/reference/PairedRR.md)

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

The vertical axis counts the pairs in which only treatment 1 responds
and runs downwards, and the horizontal axis counts the pairs in which
only treatment 2 responds, as in the map printed by
[`print.bbssr_prr`](https://gosukehommaex.github.io/bbssr/reference/print.bbssr_prr.md).
Only the possible outcomes, with at most `N` discordant pairs, are
drawn.

## Examples

``` r
RR <- PairedRR(N = 20, alpha = 0.025, Test = 'McNemar-midP')
plot(RR)

```
