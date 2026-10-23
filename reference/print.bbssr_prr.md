# Print Method for bbssr_prr Objects

Prints the settings of a rejection region for paired binary data and,
for a small number of pairs, a map of the outcomes.

## Usage

``` r
# S3 method for class 'bbssr_prr'
print(x, show.map = NULL, ...)
```

## Arguments

- x:

  An object of class `bbssr_prr` returned by
  [`PairedRR`](https://gosukehommaex.github.io/bbssr/reference/PairedRR.md)

- show.map:

  Logical. Whether the map of the outcomes is printed. The default of
  `NULL` prints it when there are at most 20 pairs

- ...:

  Further arguments, ignored

## Value

`x`, invisibly

## Details

The map has one row for each number `a` of pairs in which only treatment
1 responds and one column for each number `b` of pairs in which only
treatment 2 responds. `X` marks a rejected outcome, `.` a retained
outcome, and the impossible outcomes `a + b > N` are left blank.

## Examples

``` r
RR <- PairedRR(N = 15, alpha = 0.05, Test = 'McNemar')
print(RR)
#> Rejection region for paired binary data
#> 
#>   Test          : McNemar
#>   Alternative   : greater
#>   Pairs         : N = 15
#>   Alpha         : 0.05
#> 
#>   Rejected outcomes: 34 of 136
#> 
#>      b=0 b=1 b=2 b=3 b=4 b=5 b=6 b=7 b=8 b=9 b=10 b=11 b=12 b=13 b=14 b=15
#> a=0  .   .   .   .   .   .   .   .   .   .   .    .    .    .    .    .   
#> a=1  .   .   .   .   .   .   .   .   .   .   .    .    .    .    .        
#> a=2  .   .   .   .   .   .   .   .   .   .   .    .    .    .             
#> a=3  X   .   .   .   .   .   .   .   .   .   .    .    .                  
#> a=4  X   .   .   .   .   .   .   .   .   .   .    .                       
#> a=5  X   .   .   .   .   .   .   .   .   .   .                            
#> a=6  X   X   .   .   .   .   .   .   .   .                                
#> a=7  X   X   X   .   .   .   .   .   .                                    
#> a=8  X   X   X   .   .   .   .   .                                        
#> a=9  X   X   X   X   .   .   .                                            
#> a=10 X   X   X   X   .   .                                                
#> a=11 X   X   X   X   X                                                    
#> a=12 X   X   X   X                                                        
#> a=13 X   X   X                                                            
#> a=14 X   X                                                                
#> a=15 X                                                                    
#> 
#>   a: pairs in which only treatment 1 responds
#>   b: pairs in which only treatment 2 responds
#>   X marks rejection of the null hypothesis
```
