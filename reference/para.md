# Example of parameter vector for the main model

This is an example of the parameter vector for the main model used in
the `hhdynamics` function. This vector specifies the format of the
parameter vector for the main model.

## Usage

``` r
data(para)
```

## Format

A vector with 7 elements, where each of them is a model parameter:

- element 1:

  the random effect of individual infectivity (not available in this
  version. Models assumed no indiviudal heterogeneity after accounting
  for factors affecting infectivity)

- element 2:

  the probability of infection from the community

- element 3:

  the probability of person-to-person transmission in households

- element 4:

  the parameter of the relationship between the number of household
  contacts and transmission (default 0 in this version)

- element 5:

  the parameter for infectivity of male (Reference group: male)

- element 6:

  the relative susceptibility of age group 1 (Reference group: age group
  0)

- element 7:

  the relative susceptibility of age group 2 (Reference group: age group
  0)

## See also

Other example_data:
[`SI`](https://timktsang.github.io/hhdynamics/reference/SI.md)
