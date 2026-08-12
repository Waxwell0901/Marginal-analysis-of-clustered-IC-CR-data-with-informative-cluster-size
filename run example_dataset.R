##=============================================================
# Example analysis using the proposed method
#
# Files used:
#   - analysis.R
#   - example_dataset.csv
#
# This script provides an example implementation of the proposed method
#   using the simulated dataset.
#
# Tip: Before running this script, set the working directory to the folder
#   containing analysis.R and example_dataset.csv.
#==============================================================

### 1. Load the estimation function

source("analysis.R")

### 2. Load the example dataset

example_data<-read.csv("example_dataset.csv")
 
### 3. Inspect the example dataset

head(example_data)

## Example output:
#
#     identi  id  sizenum     v          u        T    J  X1          X2
# 1      1    1    12      4.478684  5.834566 4.972040 1  0    -0.33470558
# 2      1    2    12     13.979183      Inf      Inf  0  0    -0.23870859
# 3      1    3    12     13.917706      Inf      Inf  0  0    -0.13696193
# 4      1    4    12      6.707961  7.267917 6.838714 1  1    -0.04307298
# ...

### 4. Specify the transformation parameters
## Here, alpha1 = alpha2 = 0 is used as an illustrative example,
##   corresponding to the Fine-Gray model.

alpha1<-0
alpha2<-0

### 5. Fit the proposed model
  
fit<-analysis(simulation=example_data,alpha1=alpha1,alpha2=alpha2)

## Example output:
# "Convergence successful"

### 6. Display the regression coefficient estimates 
p<-length(grep("^X[0-9]+$",names(example_data)))

beta1_hat<-fit[1:p]
beta2_hat<-fit[(p+1):(2*p)]

cat("Regression coefficient estimates for cause 1:\n",beta1_hat)

cat("Regression coefficient estimates for cause 2:\n",beta2_hat)

## Example output:
# Regression coefficient estimates for cause 1:
# 0.5011788 -0.3372753
# Regression coefficient estimates for cause 2:
# -0.4701433 0.3490244 