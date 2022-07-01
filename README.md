# clusterpval: Inference for Estimated Clusters in R

See [http://lucylgao.com/clusterpval](http://lucylgao.com/clusterpval) for tutorials and examples. See [https://arxiv.org/abs/2012.02936](https://arxiv.org/abs/2012.02936) for the preprint.

Install 
-----

Make sure that ``devtools`` is installed by running ``install.packages("devtools")``, then type

```R
devtools::install_github("lucylgao/clusterpval")
```

## Forking to try to prevent float overflow
```R
> i = 1
> while (gamma(i/2) != Inf) {i = i + 1}
> i
[1] 344
```

Any dataset with more than 343 cols will have a log_survives of -Infinity for every perturbed dataset, and NaN values for pval and stderr. I propose writing a function that does the log and the gamma at the same time to avoid overflowing a float.
