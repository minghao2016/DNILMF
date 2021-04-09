# DNILMF

2021-04-09

I modified the demo_DNILMF.R file which originally used the information at time (t - 1) when computing similarity matrix at time t in cv. The better way is to compute similarity matrix in each fold independently (we appreciate Dr. Soheil Jahangiri who noted this issue). Note the new modified demo file will give the results with lower accuracy compared with the original one (calculated by demo_DNILMF_old.R). Other demo files were also modified accordingly.

One should try to tune the parameters to see if better performance can be achived.


2017-07-18

Add two cpp functions (sigmoid.cpp and log1pexp.cpp) for numerical stability.



The R codes were used to predict drug-target interactions (DTI) by the proposed dual-network integrated logistic matrix factorization [1].

How To Run:

(1) Install required packages

(2) setwd("dir including source codes")

(3) source("demo_DNILMF.R") ## GPCR data was used as an example

These codes were tested on Windows 7, but should also work on Linux.

For any problem about these codes, please feel free to contact Dr. Ming Hao at: kevin.m.hao@gmail.com.

Reference:

[1] M Hao, et al., Predicting drug-target interactons by dual-network integrated logistic matrix factorization (accepted).
