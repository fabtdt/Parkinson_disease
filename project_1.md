Network-based Data Analysis
================

# Introduction

The following work aims at understanding the characterization of
candidate biomarkers with respect to Parkinson’s disease by using gene
expression scans in blood.

This task will be pursued by analyzing two similar Dataset: *GSE6613*
and *GSE72267*.

In particular:

-   The Dataset *GSE6613* was produced in 2006 by performing a
    trascriptome scan in 105 individuals (50 with Parkinson’s disease,
    33 with neurogenerative diseases and 23 healthy controls). The gene
    expression profile was computed using the platform: \[HG-U133A\]
    Affymetrix Human Genome U133A Array.
-   The Dataset *GSE72267* was produces in 2015 and as well as the other
    it tries to identify changes in gene expression in 59 individuals
    (40 with Parkinson’s Disease and 19 controls). The gene expression
    profile was computed using the platform: \[HG-U133A\_2\] Affymetrix
    Human Genome U133A 2.0 Array.

### Libraries

All the libraries that will be used are here reported in the table:

| Library    | Version |
|------------|---------|
| GEOquery   | 2.58.0  |
| ggplot2    | 3.3.3   |
| ggfortify  | 0.4.11  |
| useful     | 1.2.6   |
| ggpubr     | 0.4.0   |
| factoextra | 1.0.7   |

# Data Exploration

The Data are retrieved from *Gene Expression Omnibus (GEO)*. Dataset has
been queried only the first time and then saved locally.

**Loading the data**

``` r
library ("GEOquery")
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
# The following lines are retrieve the Datset and store them locally, hence there is the need of running them only the first time

#gse_06<- getGEO("GSE6613", destdir = ".", getGPL = FALSE)
#gse_15<- getGEO("GSE72267", destdir = ".", getGPL = FALSE)

# Loading the Dataset

gse_06<- getGEO(file = "GSE6613_series_matrix.txt.gz", getGPL = FALSE)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   .default = col_double(),
    ##   ID_REF = col_character()
    ## )
    ## i Use `spec()` for the full column specifications.

``` r
gse_15<- getGEO(file = "GSE72267_series_matrix.txt.gz", getGPL = FALSE)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   .default = col_double(),
    ##   ID_REF = col_character()
    ## )
    ## i Use `spec()` for the full column specifications.

**Inspecting the Dataset**

``` r
show(gse_06)
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 22283 features, 105 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM153404 GSM153405 ... GSM153508 (105 total)
    ##   varLabels: title geo_accession ... data_row_count (29 total)
    ##   varMetadata: labelDescription
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 17215369
    ## 18669654 
    ## Annotation: GPL96

``` r
show(gse_15) 
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 22277 features, 59 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM1859079 GSM1859080 ... GSM1859137 (59 total)
    ##   varLabels: title geo_accession ... tissue:ch1 (33 total)
    ##   varMetadata: labelDescription
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 26510930 
    ## Annotation: GPL571

*exprs()* allows to get the data about the level of expression from the
Large ExpressionSet data structure.

``` r
# ex_06 and ex_15 will be the Dataset that will contain the data of the level of expression
ex_06 <- as.data.frame(exprs(gse_06))
ex_15 <- as.data.frame(exprs(gse_15))

dim(ex_06) # rows: 22 283, column: 105
```

    ## [1] 22283   105

``` r
dim(ex_15) # rows: 22 277, column: 59
```

    ## [1] 22277    59

``` r
#head(ex_06)
#head(ex_15)
```

**Harmonize the two Dataset**

The two Dataset were produced by the means of two slightly different
Affimetrix array, hence there is the need of making the two two Dataset
compatible.

``` r
# Discard the rows that has not a corresponding gene in the other Dataset
ex_06<-ex_06[which(rownames(ex_06) %in% rownames(ex_15)) ,]
ex_15<-ex_15[which(rownames(ex_15) %in% rownames(ex_06)),]
```

Now the two Datasets will be made of:

-   *GSE6613* –&gt; rows(genes): 22 227 and columns(individuals): 105
-   *GSE72267* –&gt; rows(genes): 22 227 and columns(individuals): 59

That means 6 rows has been discarded from *GSE6613*.

The following step is to check for the presence of *NA*.

``` r
sum(is.na.data.frame(ex_06))
```

    ## [1] 0

``` r
sum(is.na.data.frame(ex_06))
```

    ## [1] 0

Luckily there are no missing values in either Datasets.

**Boxplot**

As part of the data exploration it is important to boxplot the
distribution of the level of expression for each individual. This
procedure is pivotal for:

1.  check for previous transformation/normalization of the data
2.  check the data distribution in order to verify if normalization or
    transformation is needed (important procedure for reducing
    non-biological noise in the data)

``` r
par(mfrow=c(1,2))
boxplot(ex_06, main = "GSE6613", xlab = "Individuals", ylab = "Level of expression")
boxplot(ex_15, main = "GSE72267", xlab = "Individuals", ylab = "Level of expression")
```

![](project_1_files/figure-gfm/boxplot-1.png)<!-- -->

The two boxplot seems to be quite different one another. In particular
the plot regarding the Dataset *GSE72267* seems to be log-transformed
and cleaned from the outliers.

On the other hand the plot regrading the Dataset *GSE6613* seems to be
noisy and it will require some cleaning.

In order to reduce

``` r
for(i in 1:ncol(ex_06)) {       # for-loop over columns
  ex_06[ , i] <- log(ex_06[ , i])
}

boxplot(ex_06)
```

![](project_1_files/figure-gfm/log-transofmr%20GSE6613-1.png)<!-- -->

After the log-transformation now the boxplot are aligned. Nevertheless I
decided to add a constant 3 to all the data in order to have only
positive data and in almost the same rang as the one of the other
Dataset.

``` r
ex_06 <- ex_06 + 3

boxplot(ex_06)
```

![](project_1_files/figure-gfm/translate%20+3%20GSE6613-1.png)<!-- -->

Now the two Dataset are homogeneous.

# Data Analysis

**Data Preparation**

Preparing the Dataset for the PCA. In particular:

1.  transposing the Dataset
2.  adding a column that will tel whether the sample comes from an
    healthy person or not

*GSE72267*

``` r
t_ex_15 <- as.data.frame(t(na.omit(ex_15))) # transposing the Dataset and excluding NA

# creating a column containing the type of the individual
t_ex_15$Diagnosis <- c(gse_15$characteristics_ch1)

# renaming the type 
t_ex_15$Diagnosis[t_ex_15$Diagnosis == "diagnosis: Healthy"] <- "Control"
t_ex_15$Diagnosis[t_ex_15$Diagnosis == "diagnosis: Parkinson's disease"] <- "Parkinson's disease"

t_ex_15$Diagnosis <- as.factor(t_ex_15$Diagnosis)

#head(t_ex_15)
```

*GSE6613* This Dataset contains also 33 individuals with other
neurological disease other than Parkinson, this individuals will be
discarded.

``` r
t_ex_06 <- as.data.frame(t(na.omit(ex_06))) # transposing the Dataset and excluding NA

# creating a column containing the type of the individual
t_ex_06$Diagnosis <- c(gse_06$characteristics_ch1)
t_ex_06 <-t_ex_06[!(t_ex_06$Diagnosis == "neurological disease control"),]
t_ex_06$Diagnosis[t_ex_06$Diagnosis == "healthy control"] <- "Control"

t_ex_06$Diagnosis <- as.factor(t_ex_06$Diagnosis)

#head(t_ex_06)
```

## Dimensionality reduction -&gt; PCA

PCA is a method of obtaining important variables (in form of components)
from a large set of variables available in a data set. It extracts low
dimensional set of features by taking a projection of irrelevant
dimensions from a high dimensional data set with a motive to capture as
much information as possible.

*GSE72267*

``` r
library("ggfortify")
```

    ## Loading required package: ggplot2

``` r
library("ggplot2") # library needed for autoplot
pca_15 <- prcomp(t_ex_15[,-which(colnames(t_ex_15)=="Diagnosis")], scale = T) 
# not considering diagnosis because it is not numeric

par(mfrow = c(2, 1))
autoplot(pca_15,choices = c(1,3) , data = t_ex_15, colour = "Diagnosis")
```

    ## Warning: `select_()` was deprecated in dplyr 0.7.0.
    ## Please use `select()` instead.

![](project_1_files/figure-gfm/PCA%20GSE72267-1.png)<!-- -->

``` r
autoplot(pca_15,choices = c(1,3) , data = t_ex_15, x= 1, y = 3, colour = "Diagnosis")
```

![](project_1_files/figure-gfm/PCA%20GSE72267-2.png)<!-- -->

*GSE6613*

``` r
pca_06 <- prcomp(t_ex_06[,-which(colnames(t_ex_06)=="Diagnosis")], scale = T) 
# not considering diagnosis because it is not numeric

par(mfrow = c(2, 1))
autoplot(pca_06,choices = c(1,3) , data = t_ex_06, colour = "Diagnosis")
```

![](project_1_files/figure-gfm/PCA%20GSE6613-1.png)<!-- -->

``` r
autoplot(pca_06,choices = c(1,3) , data = t_ex_06, x= 1, y = 3, colour = "Diagnosis")
```

![](project_1_files/figure-gfm/PCA%20GSE6613-2.png)<!-- -->

PCA has not highlighted any underlying structure of the data. This could
be due to the fact that PCA is a linear algorithm, hence it cannot
represent complex relationship between features.

## Clustering

Clustering is an unsupervised learning algorithm and it aims at
identifying cluster within the data point by computing (dis)similarity
measure between clusters and datapoints.

The types of clustering that will be seen are:

-   K-means: here the number of cluster is represented by the
    hyperparameter k indeed k is choose via hyperparameter tuning or
    give any prior knowledge (as here).
-   hierarchical clustering: here it will be used top-down approach

#### K-means Clustering

*GSE72267*

``` r
library("useful")
```

    ## Registered S3 methods overwritten by 'useful':
    ##   method         from     
    ##   autoplot.acf   ggfortify
    ##   fortify.acf    ggfortify
    ##   fortify.kmeans ggfortify
    ##   fortify.ts     ggfortify

``` r
library("factoextra")
```

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
kmeans_15 <- kmeans(t_ex_15[,-which(colnames(t_ex_15)=="Diagnosis")], 2)
#table(kmeans_15$cluster)


fviz_cluster(kmeans_15, data=t_ex_15[,-which(colnames(t_ex_15)=="Diagnosis")], 
             geom = "point", shape = t_ex_15$Diagnosis, pointsize = 3,
             show.clust.cent = FALSE, ellipse.type = "convex", ellipse = TRUE, ellipse.level = 0.95,  ellipse.alpha = 0.06, )  +
   theme_classic()+ theme(panel.border = element_rect(colour = "black", fill=NA))
```

    ## Warning in if (shape %in% colnames(data)) {: la condizione la lunghezza > 1 e
    ## solo il promo elemento verrà utilizzato

![](project_1_files/figure-gfm/k-means%20GSE72267-1.png)<!-- -->

*GSE6613*

``` r
kmeans_06 <- kmeans(t_ex_06[,-which(colnames(t_ex_06)=="Diagnosis")], 2)
table(kmeans_06$cluster)
```

    ## 
    ##  1  2 
    ## 55 17

``` r
fviz_cluster(kmeans_06, data=t_ex_06[,-which(colnames(t_ex_06)=="Diagnosis")], 
             geom = "point", shape = t_ex_06$Diagnosis, pointsize = 3,
             show.clust.cent = FALSE, ellipse.type = "convex", ellipse = TRUE, ellipse.level = 0.95,  ellipse.alpha = 0.06, )  +   
      theme_classic()+ theme(panel.border = element_rect(colour = "black", fill=NA))
```

    ## Warning in if (shape %in% colnames(data)) {: la condizione la lunghezza > 1 e
    ## solo il promo elemento verrà utilizzato

![](project_1_files/figure-gfm/k_means%20GSE6613-1.png)<!-- -->

Points and Triangles represent the diagnosis (Control or Parkinson’s
Disease). It is possible to see that for both the dataset K-means
clustering is not able to effectivelly separate healthy patient from the
other ones.

### Hierarchical Clustering

*GSE72267*

``` r
hc_15 <- hclust(dist(t_ex_15[,-which(colnames(t_ex_15)=="Diagnosis")]), method = "ave")
hc_15$labels <- t_ex_15$Diagnosis # assign as label the diagnosis

fviz_dend(hc_15, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
```

![](project_1_files/figure-gfm/hierarchical%20GSE72267-1.png)<!-- -->

*GSE6613*

``` r
hc_06 <- hclust(dist(t_ex_06[,-which(colnames(t_ex_06)=="Diagnosis")]), method = "ave")
hc_06$labels <- t_ex_06$Diagnosis # assign as label the diagnosis

fviz_dend(hc_06, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
```

![](project_1_files/figure-gfm/hierarchical%20GSE6613-1.png)<!-- --> I
decided to use 4 groups because the algorithm identify two clusters that
are compose only by a single individual. It is possible to state that,
as well as k-means clustering, also hierarchical clustering is not able
to correctly identify the diagnosis of patient. This is probably due to
the fact that it uses the output of PCA as a starting point for
clustering.

## Random forest

Random forest is basically an ensemble of trees consisting of bagging of
un-pruned decision tree learners with a randomized selection at each
step. Here one important hyperparameter that needs to be tuned is the
number of trees to grow.

``` r
t_ex_15$Diagnosis <- as.factor(t_ex_15$Diagnosis)
summary(t_ex_15$Diagnosis)
```

    ##             Control Parkinson's disease 
    ##                  19                  40

In order to evaluate how well a classification tree performs, we first
randomly split the Dataset into two different sets, one used for fitting
the tree on our train data and the other used for validation and test
error estimation.

Creating a random sample from the dataset. This random sample will
contain 40 individuals and it will be used as training set. The other 19
individual will be used as a validation set.

``` r
set.seed(10)

train <- sample(nrow(t_ex_15), 40)
```

``` r
library("randomForest")
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
rf <- randomForest(x=t_ex_15[train,-which(colnames(t_ex_06)=="Diagnosis")], y=t_ex_15[train,]$Diagnosis, ntree=1000)
print("rf done")# it takes about 20 seconds for ntree=10000
```

    ## [1] "rf done"

Because random forest use a sample of the feature each time it trains a
tree (in order to improve the randomness), this metho performs also a
sort of feautre selection. Indeed it is possible to rank feature
according, for example, to the mean decrease in Gini index.

``` r
varImpPlot(rf, main = "Random Forests - Variable importance")
```

![](project_1_files/figure-gfm/unnamed-chunk-1-1.png)<!-- --> Now it is
possible to test the goodness of the random forest by testing it on
unseen data.

``` r
rf_error <- rf$err.rate[nrow(rf$err.rate), 1]
cat("Random forest test error is: ", rf_error, "\n")
```

    ## Random forest test error is:  0.375

Random forest does not perform so well since it has a quite high error
rate.

For RF there is no need for computing the test error with CV, instead it
is possible to look at the OOB classification error to estimate it.
Indeed, for any observation, the OOB predicts the response using all
models that do not include that observation.

``` r
rf
```

    ## 
    ## Call:
    ##  randomForest(x = t_ex_15[train, -which(colnames(t_ex_06) == "Diagnosis")],      y = t_ex_15[train, ]$Diagnosis, ntree = 1000) 
    ##                Type of random forest: classification
    ##                      Number of trees: 1000
    ## No. of variables tried at each split: 149
    ## 
    ##         OOB estimate of  error rate: 37.5%
    ## Confusion matrix:
    ##                     Control Parkinson's disease class.error
    ## Control                   6                  11   0.6470588
    ## Parkinson's disease       4                  19   0.1739130

The OOB error rate is high, that means RF does not perform so well.

``` r
plot(rf)
```

![](project_1_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
imp <- importance(rf)
impsor <- sort(imp[, 1], decreasing=TRUE)
plot(impsor)
```

![](project_1_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Linear Discriminant Functions

LDA finds the most discriminant projection by maximizing between-class
distance and minimizing within-class distance. It allows to reduce the
dimensionality and still preserve the ability to discriminate.

At first there is the need to implement some sort of gene filtering.

``` r
library("genefilter")

ex_15.2 <- exprs(gse_15)

tt_40 <- rowttests(ex_15.2, factor(t_ex_15$Diagnosis))

summary(tt_40$p.value)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 0.0000091 0.2198640 0.4743269 0.4799073 0.7333589 0.9999995

According to the values reported on the summary, the p-value threshold
is choosen.

``` r
keepers <- which(tt_40$p.value<0.1) # Change the P-value for choosing different subset of genes

t_ex_15.2 <- as.data.frame(t(ex_15.2[keepers,]))
t_ex_15.2$Diagnosis <- factor(t_ex_15$Diagnosis)

dim(t_ex_15.2)
```

    ## [1]   59 2790

The number of genes has been reduced

``` r
set.seed(2)

train <- sample(nrow(t_ex_15.2), 40)
```

``` r
library("MASS")
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:genefilter':
    ## 
    ##     area

``` r
mod <- lda(Diagnosis ~ ., data=t_ex_15.2, prior = c(0.5,0.5), subset = train)
```

    ## Warning in lda.default(x, grouping, ...): variables are collinear

``` r
plot(mod)
```

![](project_1_files/figure-gfm/perform%20LDA-1.png)<!-- -->

From this plot is possible to see that LDA seems to be able of
projecting LDA and Parkinson’s disease in a fine manner.

``` r
mod.values <- predict(mod, t_ex_15.2[train,])
mod.values$class
```

    ##  [1] Control             Control             Parkinson's disease
    ##  [4] Parkinson's disease Parkinson's disease Parkinson's disease
    ##  [7] Parkinson's disease Control             Parkinson's disease
    ## [10] Control             Parkinson's disease Parkinson's disease
    ## [13] Control             Parkinson's disease Parkinson's disease
    ## [16] Parkinson's disease Parkinson's disease Control            
    ## [19] Control             Parkinson's disease Control            
    ## [22] Parkinson's disease Parkinson's disease Parkinson's disease
    ## [25] Control             Control             Parkinson's disease
    ## [28] Control             Parkinson's disease Parkinson's disease
    ## [31] Parkinson's disease Parkinson's disease Parkinson's disease
    ## [34] Control             Control             Control            
    ## [37] Control             Parkinson's disease Parkinson's disease
    ## [40] Parkinson's disease
    ## Levels: Control Parkinson's disease

``` r
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1],
col=c(as.numeric(t_ex_15.2[train,"Diagnosis"])+10))
```

![](project_1_files/figure-gfm/unnamed-chunk-6-1.png)<!-- --> CV in
roder to test the performance of LDA.

``` r
preds<-predict(mod, t_ex_15.2[-train,])
#preds$class

table(as.numeric(preds$class),
as.numeric(t_ex_15.2[-train, "Diagnosis"]) )
```

    ##    
    ##      1  2
    ##   1  5  1
    ##   2  1 12

By looking at the confusion matrix, LDA seems to perform well with only
two misclassification example out of the 19 test individuals.

Plot the ROC curve.

``` r
library("pROC")
```

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     var

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
roc_lda <- plot.roc(as.numeric(preds$class),
as.numeric(t_ex_15.2[-train, "Diagnosis"]) )
```

    ## Setting levels: control = 1, case = 2

    ## Setting direction: controls < cases

![](project_1_files/figure-gfm/ROC%20curve-1.png)<!-- -->

Try to perform 10 fold CV by using caret

``` r
library("caret")
```

    ## Loading required package: lattice

``` r
library("e1071")

control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

fit.lda <- train(Diagnosis~., data=t_ex_15.2, method="lda", metric=metric, trControl=control)
```

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

``` r
fit.rf <- train(Diagnosis ~., data=t_ex_15.2, method="rf", metric=metric, trControl=control)

results <- resamples(list(LDA=fit.lda, RF=fit.rf))

summary(results)
```

    ## 
    ## Call:
    ## summary.resamples(object = results)
    ## 
    ## Models: LDA, RF 
    ## Number of resamples: 10 
    ## 
    ## Accuracy 
    ##          Min.   1st Qu.    Median      Mean   3rd Qu. Max. NA's
    ## LDA 0.8333333 0.8333333 1.0000000 0.9333333 1.0000000    1    0
    ## RF  0.8000000 0.8333333 0.8333333 0.8466667 0.8333333    1    0
    ## 
    ## Kappa 
    ##          Min.   1st Qu.    Median      Mean   3rd Qu. Max. NA's
    ## LDA 0.5714286 0.5952381 1.0000000 0.8380952 1.0000000    1    0
    ## RF  0.0000000 0.5714286 0.5714286 0.5571429 0.5714286    1    0

``` r
ggplot(results) + labs(y = "Accuracy") 
```

![](project_1_files/figure-gfm/unnamed-chunk-8-1.png)<!-- --> run
algorithm using 10-fold CV, 10 times

``` r
control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

fit.lda.2 <- train(Diagnosis~., data=t_ex_15.2, method="lda", metric=metric, trControl=control)
```

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

    ## Warning in lda.default(x, grouping, ...): variables are collinear

``` r
fit.rf.2 <- train(Diagnosis~., data=t_ex_15.2, method="rf", metric=metric, trControl=control)

results <- resamples(list(LDA=fit.lda.2, RF=fit.rf.2))
```

``` r
ggplot(results) + labs(y = "Accuracy")
```

![](project_1_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
