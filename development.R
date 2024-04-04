
library("devtools")
load_all()


## example 1
model1 <- system.file("extdata", "rb.fig1.model", package = "MPTinR")
m1_1 <- read_mpt(model1)

model1.eqn <- system.file("extdata", "rb.fig1.model.eqn", package = "MPTinR")
m1_2 <- read_mpt(model1.eqn)
all.equal(m1_1$path, m1_2$path)

model1.txt <- "p * q * r
p * q * (1-r)
p * (1-q) * r
p * (1-q) * (1-r) + (1-p)"

m1_3 <- read_mpt(text = model1.txt)
all.equal(m1_1, m1_3)

mp <- parse_model_df(m1_3)
mp

find.MPT.params(mp)
check.MPT.probabilities(mp)

### Example 2
?TreeBUGS::traitMPT
EQNfile <- system.file("MPTmodels/2htsm.eqn", package = "TreeBUGS")
tmp <- read_mpt(EQNfile)

mp <- parse_model_df(tmp)
find.MPT.params(mp)
check.MPT.probabilities(mp)

##  install.packages("randtoolbox")
## install.packages("fastGHQuad")

install.packages("rtmb_try/mptmem_0.04.tar.gz", repos = NULL, type = "source")
