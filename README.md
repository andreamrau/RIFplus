# RIFplus

To install package, use the following command

```
library(devtools)
install_github("andreamrau/RIFplus")
```

Here's an example of some code to run:
```
library(RIFplus)
data(RIF_data)
exprs <- RIF_data
conds <- factor(c(rep("nain", 24), rep("wt", 24)))
data(RIF_genenames)
rownames(exprs) <- RIF_genenames$V2
TFnames <- rownames(exprs)[1:754]
rif_calc <- RIFplus(exprs=exprs, TFnames=TFnames, conds=conds)
rif <- cbind(rif_calc$rif1, rif_calc$rif2)
```