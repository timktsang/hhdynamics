# Example of input data

This is an example of the input data used in the `hhdynamics` function.
This data frame illustrates the format of the input data. This is a
simulated data.

## Usage

``` r
data(inputdata)
```

## Format

A example data with 8 variables, where each row represents an
individual. For user's data, more predictors could be added by adding
columns. For categorical variables, it should be specificed as factor,
by using as.factor() function. The date is the dataset is an integer, by
selecting a reference date as day 1. In the example dataset, 2008-01-01
is day 1.

- hhID:

  The household id

- member:

  The member id in a households, use 0 to index, 1 and so on as
  household contact

- size:

  The number of individuals in the household

- end:

  The end date of follow-up of that individual

- inf:

  The infection status of the member. By defintion, index must be
  infected

- onset:

  The onset time of infected individual

- age_group:

  The age group of individual. 0: 0-19, 1: 20-64, 2: 65 or above

- sex:

  The sex of the individual. 0: Female, 1: Male

## Examples

``` r
data(inputdata)
str(inputdata)
#> 'data.frame':    1533 obs. of  8 variables:
#>  $ hhID  : int  1 1 1 100 100 100 100 10005 10005 10005 ...
#>  $ member: int  0 1 2 0 1 2 3 0 1 2 ...
#>  $ size  : int  3 3 3 4 4 4 4 3 3 3 ...
#>  $ end   : int  15 15 15 66 66 66 66 922 922 922 ...
#>  $ inf   : int  1 0 0 1 0 0 0 1 0 0 ...
#>  $ onset : int  8 NA NA 57 NA NA NA 913 NA NA ...
#>  $ age   : Factor w/ 3 levels "0","1","2": 1 1 3 2 1 1 3 1 3 2 ...
#>  $ sex   : Factor w/ 2 levels "0","1": 2 2 2 1 2 1 1 1 2 2 ...
head(inputdata)
#>   hhID member size end inf onset age sex
#> 1    1      0    3  15   1     8   0   1
#> 2    1      1    3  15   0    NA   0   1
#> 3    1      2    3  15   0    NA   2   1
#> 4  100      0    4  66   1    57   1   0
#> 5  100      1    4  66   0    NA   0   1
#> 6  100      2    4  66   0    NA   0   0
table(inputdata$inf)
#> 
#>    0    1 
#> 1055  478 
```
