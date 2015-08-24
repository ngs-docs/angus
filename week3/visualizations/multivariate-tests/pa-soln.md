```R
adonis(decostand(dataset[,-c(1:3)],method="pa")~dataset$fly*dataset$type,method="jaccard")

Call:
adonis(formula = decostand(dataset[, -c(1:3)], method = "pa") ~      dataset$fly * dataset$type, method = "jaccard") 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                         Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
dataset$fly               2  0.016070 0.0080348  1.5635 0.25265  0.001 ***
dataset$type              1  0.006278 0.0062776  1.2216 0.09870  0.115    
dataset$fly:dataset$type  2  0.010423 0.0052115  1.0141 0.16387  0.433    
Residuals                 6  0.030833 0.0051389         0.48478           
Total                    11  0.063603                   1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```