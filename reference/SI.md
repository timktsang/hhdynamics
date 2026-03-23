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

## Examples

``` r
data(SI)
print(SI)
#>  [1] 5.872190e-02 3.155084e-01 4.863670e-01 1.369000e-01 2.502287e-03
#>  [6] 3.537313e-07 1.257039e-14 0.000000e+00 0.000000e+00 0.000000e+00
#> [11] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
barplot(SI, names.arg = seq_along(SI),
        xlab = "Days since infector onset",
        ylab = "Probability",
        main = "Serial interval distribution (influenza)")
```
