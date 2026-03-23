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
