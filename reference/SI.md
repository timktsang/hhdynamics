# Example of serial interval distribution of influenza

This is an example of the serial interval distribution used in the
`hhdynamics` function. This vector specifies the format of the serial
interval distribution. This is estimated from Tsang et al. Association
between antibody titers and protection against influenza virus infection
within households. J Infect Dis. 2014 Sep 1;210(5):684-92.

## Usage

``` r
data(SI)
```

## Format

A numeric vector of length 14. Element `i` gives the probability that
symptom onset occurs `i` days after the infector's onset. The vector
sums to 1.

## See also

Other example_data:
[`para`](https://timktsang.github.io/hhdynamics/reference/para.md)
