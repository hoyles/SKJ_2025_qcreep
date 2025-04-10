# run the diagnostic assessment
# S2M1D2R1G1
# install.packages("iterators")
# install.packages("FLCore", repos="http://flr-project.org/R")
# install.packages("remotes")
# install.github("remotes")

library(tidyverse)
library(remotes)
#install_github("PacificCommunity/ofp-sam-flr4mfcl")
#install_github("PacificCommunity/ofp-sam-r4mfcl")
library(FLR4MFCL)
library(R4MFCL)
library(mgcv)
library(data.table)

# Run Linux mfclo64 in wsl environment
skj22_dir <- getwd()

basedir <- paste0(skj22_dir, "/Runs/")
testdir <- paste0(basedir, "/000_test/")
setwd(testdir)

shell("wsl ./doitall.skj >stdout.log", wait = FALSE)
shell("wsl ./doitall.skj2 >stdout2.log", wait = FALSE)
setwd(skj22_dir)


# Check which quarter has highest CPUE in each region, based on 2012 CPUE. Turns out to be the 3rd quarter in every region.  
#### new runs with WCPFC16
source_basedir <- paste0(basedir, "/000_files/")
sourcenames <- c("skj.frq", "skj.ini", "00.par", "doitall.skj", "mfcl.cfg", "mfclo64")
runfiles <- paste0(source_basedir, "/", c("skj.frq", "skj.tag", "skj.ini", "doitall.skj", "mfcl.cfg", "mfclo64", "labels.tmp"))

diagn_dir <- paste0(basedir, '00_diagnostic/')

############
run01_drop_SSAP <- list()
run01_drop_SSAP$dir <- paste0(basedir, "/run01_drop_SSAP")
dir.create(run01_drop_SSAP$dir)
setwd(run01_drop_SSAP$dir)

file.copy(from=runfiles, to=run01_drop_SSAP$dir, overwrite=TRUE)
run01_drop_SSAP$tag <- read.tag(paste(run01_drop_SSAP$dir,"skj.tag",sep="/"))     
dropgrps <- (1:328)[run01_drop_SSAP$tag$rel$y < 1985]
run01_drop_SSAP$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run01_drop_SSAP$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run01_drop_SSAP$tag$rel.recov <- run01_drop_SSAP$tag$rel.recov[!run01_drop_SSAP$tag$rel.recov$grp %in% dropgrps,]

#a <- run01_drop_SSAP$doitall <- readLines("doitall.skj")
#write(a, "doitall.skj")
write.tag("skj.tag", run01_drop_SSAP$tag)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

############
run02_drop_earlysizes <- list()
run02_drop_earlysizes$dir <- paste0(basedir, "/run02_drop_earlysizes1/")
dir.create(run02_drop_earlysizes$dir)
setwd(run02_drop_earlysizes$dir)

file.copy(from=runfiles, to=run02_drop_earlysizes$dir, overwrite=TRUE)
run02_drop_earlysizes$frq <- read.frq(paste0(run02_drop_earlysizes$dir, "skj.frq"))
a <- run02_drop_earlysizes$frq
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
a$mat[pos,9:61] <- 0
# pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(17)
# a$mat[pos,8] <- -1
# a$mat[pos,9:61] <- 0
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
a$mat[pos,9:61] <- 0

write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

############
run03_SSAP_sizes <- list()
run03_SSAP_sizes$dir <- paste0(basedir, "/run03_SSAP_sizes/")
dir.create(run03_SSAP_sizes$dir)
setwd(run03_SSAP_sizes$dir)

file.copy(from=runfiles, to=run03_SSAP_sizes$dir, overwrite=TRUE)
run03_SSAP_sizes$tag <- read.tag(paste(run03_SSAP_sizes$dir,"skj.tag",sep="/"))     
dropgrps <- (1:328)[run03_SSAP_sizes$tag$rel$y < 1985]
run03_SSAP_sizes$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run03_SSAP_sizes$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run03_SSAP_sizes$tag$rel.recov <- run03_SSAP_sizes$tag$rel.recov[!run03_SSAP_sizes$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run03_SSAP_sizes$tag)
a <- run03_SSAP_sizes$frq <- read.frq(paste0(run03_SSAP_sizes$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
a$mat[pos,9:61] <- 0
# pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(17)
# a$mat[pos,8] <- -1
# a$mat[pos,9:61] <- 0
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
a$mat[pos,9:61] <- 0
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

############
run04_add_latesizes <- list()
run04_add_latesizes$dir <- paste0(basedir, "/run04_add_latesizes/")
dir.create(run04_add_latesizes$dir)
setwd(run04_add_latesizes$dir)

file.copy(from=runfiles, to=run04_add_latesizes$dir, overwrite=TRUE)
run04_add_latesizes$tag <- read.tag(paste(run04_add_latesizes$dir,"skj.tag",sep="/"))     
dropgrps <- (1:328)[run04_add_latesizes$tag$rel$y < 1985]
run04_add_latesizes$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run04_add_latesizes$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run04_add_latesizes$tag$rel.recov <- run04_add_latesizes$tag$rel.recov[!run04_add_latesizes$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run04_add_latesizes$tag)
a <- run04_add_latesizes$frq <- read.frq(paste0(run04_add_latesizes$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
#a$mat[pos,9:61] <- 0
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
#a$mat[pos,9:61] <- 0
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
#a$mat[pos,9:61] <- 0
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############
run05_add_creep <- list()
run05_add_creep$dir <- paste0(basedir, "/run05_add_creep1/")
dir.create(run05_add_creep$dir)
setwd(run05_add_creep$dir)

file.copy(from=runfiles, to=run05_add_creep$dir, overwrite=TRUE)
run05_add_creep$tag <- read.tag(paste(run05_add_creep$dir,"skj.tag",sep="/"))     
dropgrps <- (1:328)[run05_add_creep$tag$rel$y < 1985]
run05_add_creep$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run05_add_creep$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run05_add_creep$tag$rel.recov <- run05_add_creep$tag$rel.recov[!run05_add_creep$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run05_add_creep$tag)
a <- run05_add_creep$frq <- read.frq(paste0(run05_add_creep$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
fcp <- c(32:35,37:38)
for (f in fcp) {
  pos <- a$mat[,4] == f
  a$mat[pos,6] <- a$mat[pos,6] * 1.02^((a$mat[pos,1]+a$mat[pos,2]/12)-1970)
}
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############
run06_add_cpuecv <- list()
run06_add_cpuecv$dir <- paste0(basedir, "/run06_add_cpuecv/")
dir.create(run06_add_cpuecv$dir)
setwd(run06_add_cpuecv$dir)

file.copy(from=runfiles, to=run06_add_cpuecv$dir, overwrite=TRUE)
run06_add_cpuecv$tag <- read.tag(paste(run06_add_cpuecv$dir,"skj.tag",sep="/"))     
dropgrps <- (1:328)[run06_add_cpuecv$tag$rel$y < 1985]
run06_add_cpuecv$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run06_add_cpuecv$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run06_add_cpuecv$tag$rel.recov <- run06_add_cpuecv$tag$rel.recov[!run06_add_cpuecv$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run06_add_cpuecv$tag)
a <- run06_add_cpuecv$frq <- read.frq(paste0(run06_add_cpuecv$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
fcp <- c(32:35,37:38)
for (f in fcp) {
  pos <- a$mat[,4] == f
  a$mat[pos,6] <- a$mat[pos,6] * 1.02^((a$mat[pos,1]+a$mat[pos,2]/12)-1970)
}
write.frq("skj.frq", a)
doit <- readLines("doitall.skj")
#doit <- change.fishflag(doit, fisheries=c(1:16),flagnum=4, newvals=2) # phase 2
doit <- change.fishflag(doit, fisheries=c(32),flagnum=92, newvals=30) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(33),flagnum=92, newvals=30) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(34),flagnum=92, newvals=30) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(35),flagnum=92, newvals=30) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(37),flagnum=92, newvals=30) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(38),flagnum=92, newvals=30) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(39),flagnum=92, newvals=40) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(40),flagnum=92, newvals=40) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(41),flagnum=92, newvals=40) # phase 8  #
doit <- change.fishflag(doit, fisheries=c(36),flagnum=92, newvals=50) # phase 8  #
doit <- add.flag(doit, flagtype=c(-(32:41)),flagnum=66, newval=0, phase=1) # phase 8  #

outfile <- file("./doitallx.skj", "wb")
writeLines(doit, outfile)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############
run07_only_RTTP <- list()
run07_only_RTTP$dir <- paste0(basedir, "/run07_only_RTTP/")
dir.create(run07_only_RTTP$dir)
setwd(run07_only_RTTP$dir)

file.copy(from=runfiles, to=run07_only_RTTP$dir, overwrite=TRUE)
run07_only_RTTP$tag <- read.tag(paste(run07_only_RTTP$dir,"skj.tag",sep="/"))     
dropgrps <- c(1:23, 53:328)
run07_only_RTTP$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run07_only_RTTP$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run07_only_RTTP$tag$rel.recov <- run07_only_RTTP$tag$rel.recov[!run07_only_RTTP$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run07_only_RTTP$tag)
a <- run07_only_RTTP$frq <- read.frq(paste0(run07_only_RTTP$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
fcp <- c(32:35,37:38)
for (f in fcp) {
  pos <- a$mat[,4] == f
  a$mat[pos,6] <- a$mat[pos,6] * 1.02^((a$mat[pos,1]+a$mat[pos,2]/12)-1970)
}
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############
run08_only_PTTP <- list()
run08_only_PTTP$dir <- paste0(basedir, "/run08_only_PTTP/")
dir.create(run08_only_PTTP$dir)
setwd(run08_only_PTTP$dir)

file.copy(from=runfiles, to=run08_only_PTTP$dir, overwrite=TRUE)
run08_only_PTTP$tag <- read.tag(paste(run08_only_PTTP$dir,"skj.tag",sep="/"))     
dropgrps <- c(1:52, 96:328)
run08_only_PTTP$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run08_only_PTTP$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run08_only_PTTP$tag$rel.recov <- run08_only_PTTP$tag$rel.recov[!run08_only_PTTP$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run08_only_PTTP$tag)
a <- run08_only_PTTP$frq <- read.frq(paste0(run08_only_PTTP$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
fcp <- c(32:35,37:38)
for (f in fcp) {
  pos <- a$mat[,4] == f
  a$mat[pos,6] <- a$mat[pos,6] * 1.02^((a$mat[pos,1]+a$mat[pos,2]/12)-1970)
}
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############
run09_only_JPTP <- list()
run09_only_JPTP$dir <- paste0(basedir, "/run09_only_JPTP/")
dir.create(run09_only_JPTP$dir)
setwd(run09_only_JPTP$dir)

file.copy(from=runfiles, to=run09_only_JPTP$dir, overwrite=TRUE)
run09_only_JPTP$tag <- read.tag(paste(run09_only_JPTP$dir,"skj.tag",sep="/"))     
dropgrps <- c(1:95)
run09_only_JPTP$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run09_only_JPTP$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run09_only_JPTP$tag$rel.recov <- run09_only_JPTP$tag$rel.recov[!run09_only_JPTP$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run09_only_JPTP$tag)
a <- run09_only_JPTP$frq <- read.frq(paste0(run09_only_JPTP$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
fcp <- c(32:35,37:38)
for (f in fcp) {
  pos <- a$mat[,4] == f
  a$mat[pos,6] <- a$mat[pos,6] * 1.02^((a$mat[pos,1]+a$mat[pos,2]/12)-1970)
}
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############
run10_only_SSAP <- list()
run10_only_SSAP$dir <- paste0(basedir, "/run10_only_SSAP/")
dir.create(run10_only_SSAP$dir)
setwd(run10_only_SSAP$dir)

file.copy(from=runfiles, to=run10_only_SSAP$dir, overwrite=TRUE)
run10_only_SSAP$tag <- read.tag(paste(run10_only_SSAP$dir,"skj.tag",sep="/"))     
dropgrps <- c(24:328)
run10_only_SSAP$tag$rel.lens[dropgrps,] <- 0 # Set releases to 0
run10_only_SSAP$tag$rel.lens[dropgrps,10] <- 1  # Except 1 small fish which is unlikely to survive to be recaptures  
run10_only_SSAP$tag$rel.recov <- run10_only_SSAP$tag$rel.recov[!run10_only_SSAP$tag$rel.recov$grp %in% dropgrps,]
write.tag("skj.tag", run10_only_SSAP$tag)
a <- run10_only_SSAP$frq <- read.frq(paste0(run10_only_SSAP$dir, "skj.frq"))
llfs <- c(3,6,9,17,21,23,27,31)
pos <- a$mat[,1] < 1980 & a$mat[,4] %in% llfs
a$mat[pos,8] <- -1
pos <- a$mat[,1] > 2019 & a$mat[,4] %in% c(3,9,17,21,27)
a$mat[pos,8] <- -1
r5ph <- c(12)
pos <- a$mat[,1] < 1982 & a$mat[,4] %in% r5ph
a$mat[pos,8] <- -1
fcp <- c(32:35,37:38)
for (f in fcp) {
  pos <- a$mat[,4] == f
  a$mat[pos,6] <- a$mat[pos,6] * 1.02^((a$mat[pos,1]+a$mat[pos,2]/12)-1970)
}
write.frq("skj.frq", a)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

#############

# filepath <- '/media/sf_SLOTH/skj_MSE/skj_21/MSE_Inputs/'
# model    <- 'S2M1D2R1G2'
model    <- '00_diagnostic/'

ini <- read.MFCLIni(paste0(basedir, model, 'skj.ini'))
frq <- read.MFCLFrq(paste0(basedir, model, 'skj.frq'))
par <- read.MFCLPar(paste0(basedir, model, '08.par'), first.yr = 1972)   #need to specify first year when reading the par file - annoyingly
rep <- read.MFCLRep(paste0(basedir, model, 'plot-08.par.rep'))

# Compare par files
stepmod5 <- stepmod6 <- list()
stepmod5$par <- read.par(paste0(skj22_dir, "/stepwiseskj21/5_Notag4/","13.par"))
stepmod6$par <- read.par(paste0(skj22_dir, "/stepwiseskj21/6_WCPO21_k3_2_2016/","13.par"))
compare_par_flags(stepmod5$par, stepmod6$par)

stepmod5$par <- read.par(paste0(skj22_dir, "/stepwiseskj21/5_Notag4/","01.par"))
stepmod6$par <- read.par(paste0(skj22_dir, "/stepwiseskj21/6_WCPO21_k3_2_2016/","01.par"))
compare_par_flags(stepmod5$par, stepmod6$par)


## depletion
xyplot(data~year, data= SBSBF0latest(rep),  type="l", ylab="SBSBF0_latest", xlab="Year", ylim=c(0,1.2))

## recruitment
xyplot(data~year, data=areaSums(seasonSums(rec(rep))),  type="l", ylab="Recruitment", xlab="Year", ylim=c(0,6e+09))

# Check SS n size data
frq <- read.frq(paste0(rundir, "skj.frq"))
a <- data.frame(frq$mat)
a$yq <- a[,1] + a[,2]/12
head(a)
a$nfreq <- apply(a[,8:107],1,sum)
names(a[,4]) <- "fs"
a[a$fishery==1,c("fishery", "yq", "nfreq")]
a[a$fishery==2,c("fishery", "yq", "nfreq")]
a[a$fishery==3,c("fishery", "yq", "nfreq")]
a[a$fishery==4,c("fishery", "yq", "nfreq")]
a[a$fishery==5,c("fishery", "yq", "nfreq")]
a[a$fishery==6,c("fishery", "yq", "nfreq")]
a[a$fishery==7,c("fishery", "yq", "nfreq")]
a[a$fishery==8,c("fishery", "yq", "nfreq")]
a[a$fishery==9,c("fishery", "yq", "nfreq")]
a[a$fishery==10,c("fishery", "yq", "nfreq")]
a[a$fishery==11,c("fishery", "yq", "nfreq")]
a[a$fishery==12,c("fishery", "yq", "nfreq")]
a[a$fishery==13,c("fishery", "yq", "nfreq")]
a[a$fishery==14,c("fishery", "yq", "nfreq")]
a[a$fishery==15,c("fishery", "yq", "nfreq")]
a[a$fishery==16,c("fishery", "yq", "nfreq")]
a[a$fishery==17,c("fishery", "yq", "nfreq")]
a[a$fishery==18,c("fishery", "yq", "nfreq")]
a[a$fishery==19,c("fishery", "yq", "nfreq")]
a[a$fishery==20,c("fishery", "yq", "nfreq")]
a[a$fishery==21,c("fishery", "yq", "nfreq")]
a[a$fishery==22,c("fishery", "yq", "nfreq")]
a[a$fishery==23,c("fishery", "yq", "nfreq")]
a[a$fishery==24,c("fishery", "yq", "nfreq")]
a[a$fishery==15,c("fishery", "yq", "nfreq")]
a[a$fishery==16,c("fishery", "yq", "nfreq")]
a[a$fishery==17,c("fishery", "yq", "nfreq")]
a[a$fishery==21,c("fishery", "yq", "nfreq")]

modnames <- c("1_DCskj18", "2_MFCL208", "3_fixgrowth_NH", "4_CPUE4", "5_Notag4", "6_WCPO21_k3_2_2016", 
              "7_WCPO21_k3_2_25", "8_S2M1D2R1G1_diagcase")
tstnames <- c("run01_est_tzero","run02_S2M2D2R1G1","run03_merge_b42000","run04__rm_earlier_R3LLsize",
              "run05__add_later_R3LLsize","run06_WCPFC16_last10troll","run07_diagn_last10troll","run07b_WCPFC16")


# Estimate all but last 4 recruitments instead of 9
rundir <- paste0(basedir, "/run21_diag_last4recs/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
setwd(rundir)
doit <- readLines(paste0(rundir, "doitall.skj"))
doit <- change.flag(doit, 1, 400, 4)
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# More flexible selectivity
rundir <- paste0(basedir, "/run24_select_flex4/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
setwd(rundir)
doit <- readLines(paste0(rundir, "doitall.skj"))
doit <- change.fishflag(doit, 999, 61, 4)
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# fix recruitments constant since 2014 (24)
rundir <- paste0(basedir, "/run27b_diag_last200recs/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
setwd(rundir)
doit <- readLines(paste0(rundir, "doitall.skj"))
doit <- change.flag(doit, 1, 400, 200)
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# Downweight size data and fix selectivity at current values
rundir <- paste0(basedir, "/run39_downwt_size_data/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
file.copy(from=paste0(diagn_dir, "06.par"), to=rundir, overwrite=TRUE)
setwd(rundir)
doit <- readLines(paste0(rundir, "doitall.skj"))
grep("PHASE6", doit)
doit <- c(
  " ./mfclo64 skj.frq 06.par 07.par -file - <<PHASE7",
  " -999 48 0    # turn off selectivity estimation",
  " -999 49 5000 # reduce weight on size data",
  "PHASE7")
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# Average biomass set to 2e6
rundir <- paste0(basedir, "/run48_fixed_biomass_scale/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
file.copy(from=paste0(diagn_dir, "06.par"), to=rundir, overwrite=TRUE)
setwd(rundir)
doit <- c(
  " ./mfclo64 skj.frq 06.par 07.par -file - <<PHASE7",
  " 1 346 2  ",
  " 1 347 4000000  ",
  " 1 173 230  ",
  " 1 174 5  ",
  " 1 187 0",
  " 1 188 0 ",
  " 1 348 10000000",
  "PHASE7")
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# Triple totpop
rundir <- paste0(basedir, "/run49_triple_totpop/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
file.copy(from=paste0(diagn_dir, "06.par"), to=rundir, overwrite=TRUE)
setwd(rundir)
parx <- readLines("06.par")
pos <- grep("# total populations scaling parameter", parx)
tp <- as.numeric(parx[pos+1])
newtp <- log(3*exp(tp))
parx[pos+1] <- newtp
writeLines(parx, paste0(rundir, "06b.par"))

doit <- c(
  " ./mfclo64 skj.frq 06b.par 07.par -file - <<PHASE7",
  " 2 31 0  ",
  " 2 32 0  ",
  "PHASE7")
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# Remove last 10 years of size data in F1, F3, F22, F23 and change last 4 quarters of cpue
# trying to identify data that make last year depleted in R1. 
rundir <- paste0(basedir, "/run50_diag_last10size_R4fisheries/")
dir.create(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
setwd(rundir)
frq <- read.frq(paste0(rundir, "skj.frq"))
a <- frq$mat
pos <- a[,1] >= 2010 & a[,4] %in% c(1,3,22,23) # Sizes
a[pos,8] <- -1
a[pos,9:107] <- 0
pos <- (1:dim(a)[1])[a[,1] > 2017 & a[,4] %in% c(25)] # CPUE
a[pos,6] <- a[pos-8,6]
a[pos,7] <- 0.1
frq$mat <- a
write.frq(paste0(rundir, '/skj.frq'), frq)
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)

# Downweight badly fitting R4 in final phase to reduce downward scaling pressure
rundir <- paste0(basedir, "/run51_downwt_r4/")
dir.create(rundir)
setwd(rundir)
file.copy(from=diag_sourcefiles, to=rundir, overwrite=TRUE)
file.copy(from=paste0(diagn_dir, "06.par"), to=rundir, overwrite=TRUE)
doit <- c(
  " ./mfclo64 skj.frq 06.par 07.par -file - <<PHASE7",
  " -18 49 500  ",
  " -19 49 500  ",
  " -20 49 500  ",
  " -25 49 500  ",
  " -18 48 0    # turn off selectivity estimation",
  " -19 48 0    # turn off selectivity estimation",
  " -20 48 0    # turn off selectivity estimation",
  "PHASE7")
writeLines(doit, paste0(rundir, "doitall.skj"))
shell("wsl ./doitall.skj >stdout.log", wait = FALSE)


# Extract likelihood info from run 27
diag_liks <- read.MFCLLikelihood(paste0(diagn_dir, "/test_plot_output"))
run01_drop_SSAP$liks <- read.MFCLLikelihood(paste0(run01_drop_SSAP$dir, "/test_plot_output"))
run02_drop_earlysizes$liks <- read.MFCLLikelihood(paste0(run02_drop_earlysizes$dir, "/test_plot_output"))
run03_SSAP_sizes$liks <- read.MFCLLikelihood(paste0(run03_SSAP_sizes$dir, "/test_plot_output"))
run04_add_latesizes$liks <- read.MFCLLikelihood(paste0(run04_add_latesizes$dir, "/test_plot_output"))
run05_add_creep$liks <- read.MFCLLikelihood(paste0(run05_add_creep$dir, "/test_plot_output"))
run06_add_cpuecv$liks <- read.MFCLLikelihood(paste0(run06_add_cpuecv$dir, "/test_plot_output"))
# Load all likelihoods
step_liks <- read.MFCLLikelihood(paste0(basedir, "/run48_fixed_biomass_scale/","test_plot_output"))

slotNames(diag_liks)
effort_dev_penalty(run27liks)- effort_dev_penalty(diag_liks)
slotNames(diag_liks)
total_length_fish(diag_liks)
lapply(length_fish(diag_liks), sum)

# Extract likelioods frm run with last 2 years 
effort_dev_penalty(run10liks)- effort_dev_penalty(diag_liks) # effect on CPUE of removing last 2 years of CPUE
effort_dev_penalty(run16liks)- effort_dev_penalty(diag_liks) # effect on CPUE of removing last 10 years of LL size data

slotNames(diag_liks2)
diag_liks2

modnames <- c("run01_drop_SSAP", "run02_drop_earlysizes","run03_SSAP_sizes", "run04_add_latesizes", "run05_add_creep", "run07_only_RTTP", "run08_only_PTTP", "run09_only_JPTP", "run10_only_SSAP")

# Effect of phase 9 in diagnostic model
# phase9 / phase8
0.2444/0.2510 # FMSY is 2.6% less
2.472e6/2.321e6  # BMSY increases by 6.5%
1.073 / 1  #SBMSY increases by 7%
2.948 / 2.971 # SB/SBMSY reduced by 0.8%
2.861 / 2.800 # Fmult increased by 2.2%
# Natural mortality changes significantly
c(0.4921,0.4623,0.4021,0.3476,0.3207,0.3272,0.3567,0.3987,0.4385,0.4629,0.4652,0.4416,0.3946,0.3424,0.3033,0.2880)/
c(0.4578,0.4381,0.3975,0.3607,0.3449,0.3562,0.3858,0.4229,0.4531,0.4656,0.4574,0.4286,0.3833,0.3363,0.3017,0.2883)


###-----------------------------------------------------------------------------
# Plot spawning potential from different runs
all_liks <- list()
for (i in 1:9) {
  all_liks[[i]] <- read.MFCLLikelihood(paste0(basedir, model=modnames[i], '/test_plot_output'))
}

testrep2 <- list()
for (i in 1:9) {
  testrep2[[i]] <- read.MFCLRep(paste0(basedir, modnames[i], '/plot-08.par.rep'))
}

testrep1 <- list()
for (i in 1:9) {
  testrep1[[i]] <- read.rep(paste0(basedir, modnames[i], '/plot-08.par.rep'))
}
testrep09 <- list()
for (i in c(1,4,5,6,7,8,9)) {
  testrep09[[i]] <- read.rep(paste0(basedir, modnames[i], '/plot-09.par.rep'))
}
testlen <- list()
for (i in 1:8) {
  testlen[[i]] <- read.fit(paste0(basedir, modnames[i], '/length.fit'))
}
slotNames(steprep2[[8]])
adultBiomass(steprep2[[8]])
adultBiomass_nofish(steprep2[[8]])
range(steprep2[[8]])["maxyear"]
eq_rec(steprep2[[8]])
eq_rec_obs(steprep2[[8]])
rec_region(steprep2[[8]])

str(steplen)
steplen[[1]]$dates

names(steprep1[[i]])
steprep2[[i]]$nTimes
steprep2[[i]]$Year1
steprep2[[i]]$yrs
steprep2[[i]]$alltimes
steprep2[[i]]$Recruitment

# Change this to R4MFCL? 
windows()
plot(1972:2021, 1972:2021, ylim = c(0, 2e7), xlab = "Year", ylab="Adult Biomass", type = "n")
for(i in 1:9){ 
  lines(testrep1[[i]]$yrs, apply(testrep1[[i]]$AdultBiomass,1,sum), col = i, lwd=2)
}

windows()
plot(1972:2021, 1972:2021, ylim = c(0, 1), xlab = "Year", ylab="Adult Biomass", type = "n")
for(i in 1:8){ 
  lines(testrep1[[i]]$yrs, apply(testrep1[[i]]$AdultBiomass,1,sum) / apply(testrep1[[i]]$AdultBiomass.nofish,1,sum), col = i, lwd=2)
}
apply(testrep1[[6]]$AdultBiomass,1,sum) / apply(testrep1[[6]]$AdultBiomass.nofish,1,sum)



windows()
plot(1972:2021, 1972:2021, ylim = c(0, 2e7), xlab = "Year", ylab="Adult Biomass", type = "n")
for(i in 1:8){ 
  a <- adultBiomass(steprep2[[i]])
  abom <- apply(a[,,,4,], 2, sum)
  yrs <- as.numeric(colnames(abom))
  lines(yrs, abom, col = i, lwd=2)
}

windows(10,10)
plot(1972:2021, 1972:2021, ylim = c(0, 1), xlab = "Year", ylab="Depletion", type = "n")
for(i in 1:8){ 
  a <- adultBiomass(steprep[[i]])
  abom <- apply(a[,,,4,], 2, sum)
  a <- adultBiomass_nofish(steprep[[i]])
  abnf <- apply(a[,,,4,], 2, sum)
  yrs <- as.numeric(colnames(abom))
  lines(yrs, abom/abnf, col = i, lwd=2)
}
legend("bottomleft", modnames, col=1:8, lwd=2, lty=1)
savePlot(paste0(skj22_dir,"/figures/","stepwise_depletion1.png"), type = "png")

windows(10,10)
plot(2000:2021, 2000:2021, ylim = c(0, 5e7), xlab = "Year", ylab="Recruitment", type = "n")
for(i in 1:8){ 
  if (i==6) x=2 else x=2
  if (i==6) ltx=2 else ltx=1
  lines(testrep1[[i]]$yrs[-(1:160)], apply(testrep1[[i]]$Recruitment,1,sum)[-(1:160)], col = i, lwd=x, lty=ltx)
}
legend("topleft", modnames, col=1:8, lwd=2, lty=c(1,1,1,1,1,2,1,1))

windows(10,10); par(mfrow=c(3,3))
for(i in 1:8){ 
  yl <- max(apply(testrep1[[i]]$Recruitment,1,sum)[-(1:160)])
  plot(2000:2021, 2000:2021, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = modnames[i])
  lines(testrep1[[i]]$yrs[-(1:160)], apply(testrep1[[i]]$Recruitment,1,sum)[-(1:160)], col = i, lwd=1, lty=1)
}
savePlot(paste0(skj22_dir,"/figures/","stepwise_recruitment.png"), type = "png")

apply(allrep2[[1]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[2]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[3]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[4]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[5]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[6]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[7]]$Recruitment,1,sum)[-(1:160)]
apply(allrep2[[8]]$Recruitment,1,sum)[-(1:160)]

# Plot CPUE
frq21 <- data.frame(read.frq(paste0(diagn_dir, "skj.frq"))$mat[,1:8])
frq12 <- data.frame(read.frq(paste0(skj22_dir, "/2012_skj_basecase/skj.frq"))$mat[,1:8])
frq15 <- data.frame(read.frq(paste0(skj22_dir, "/2015_skj_basecase/skj2015/skj.frq"))$mat[,1:8])
frq18 <- data.frame(read.frq(paste0(skj22_dir, "/SPA2018Web/GridML_GeoCPUE_GrowthEstim_SizeWgt50_Steep0.8_M0.3_diag/skj.frq"))$mat[,1:8])
frq12$f2 <- floor((frq12$fishery - .5)/4)+1
tmp <- frq12$f2[frq12$se > 0.1]
tmp <- frq12$fishery[frq12$se > 0.1]
sort(unique(tmp))

windows(9,9); par(mfrow = c(2,2))
for (ff in 1:4) {
  a <- filter(frq21, fishery==ff+21)
  a18 <- filter(frq18, fishery==ff+16 & effort >0)
  a$cpue <- a$catch / a$effort; a$cpue <- a$cpue / mean(a$cpue)
  a18$cpue <- a18$catch / a18$effort; a18$cpue <- a18$cpue / mean(a18$cpue)
  a$yq <- a$year + a$qtr/12; a18$yq <- a18$year + a18$qtr/12
  plot(a$yq, a$cpue, xlim=c(1960,2020), ylim = c(0,5), xlab = "Year", ylab = "CPUE", main = paste("Region", ff), cex=0.7)
  points(a18$yq, a18$cpue, col=2, pch=2, cex=0.7)
}
legend("topright", legend=c(2021, 2018), col=1:2, pch=1:2)
savePlot(paste0(skj22_dir,"/figures/","CPUE_diagnostic_by_assmt.png"), type = "png")

windows(9,9); par(mfrow = c(2,2))
yl <- c(0, max(frq21$se))
for (ff in 22:25) {
  a <- filter(frq21, fishery==ff)
  a$yq <- a$year + a$qtr/12
  plot(a$yq, a$se, xlim=c(1960,2020), ylim = yl, xlab = "Year", ylab = "CPUE Penalties", main = paste("Region", ff-21))
}
savePlot(paste0(skj22_dir,"/figures/","CPUE_penalties.png"), type = "png")

windows(9,9); par(mfrow = c(2,2))
yl <- c(0, max(frq21$se))
for (ff in 1:4) {
  a <- filter(frq21, fishery==ff+21)
  a18 <- filter(frq18, fishery==ff+16 & effort >0)
  a15 <- filter(frq15, fishery==ff & effort >0)
  a12 <- filter(frq12, f2==c(3,10,17,23)[ff] & effort >0)
  a$yq <- a$year + a$qtr/12; a18$yq <- a18$year + a18$qtr/12
  a15$yq <- a15$year + a15$qtr/12; a12$yq <- a12$year + a12$qtr/12
  plot(a$yq, a$se, xlim=c(1960,2020), ylim = yl, xlab = "Year", ylab = "CPUE Penalties", main = paste("Region", ff), cex=0.7)
  points(a18$yq, a18$se, col=2, pch=2, cex=0.7)
  points(a15$yq, a15$se, col=3, pch=3, cex=0.7)
#  points(a12$yq, a12$se, col=4, pch=4, cex=0.7)
}
legend("topleft", legend=c(2021, 2018, 2015), col=1:3, pch=1:3)
savePlot(paste0(skj22_dir,"/figures/","CPUE_penalties_by_assmt.png"), type = "png")


frq21 <- read.frq(paste0(diagn_dir, "skj.frq"))$mat
frq12 <- read.frq(paste0(skj22_dir, "/2012_skj_basecase/skj.frq"))$mat
frq15 <- read.frq(paste0(skj22_dir, "/2015_skj_basecase/skj2015/skj.frq"))$mat
frq18 <- read.frq(paste0(skj22_dir, "/SPA2018Web/GridML_GeoCPUE_GrowthEstim_SizeWgt50_Steep0.8_M0.3_diag/skj.frq"))$mat
frq21 <- as.data.frame(frq21[frq21[,8] > -1,])
frq18 <- as.data.frame(frq18[frq18[,8] > -1,])
frq15 <- as.data.frame(frq15[frq15[,8] > -1,])
frq12 <- as.data.frame(frq12[frq12[,8] > -1,])
frq21 <- filter(frq21, !fishery %in% c(18:21,25)) # Remove eastern
frq12 <- filter(frq12, fishery < 140) # Remove eastern
# a <- table(frq12$fishery,frq12$qtr)
# a[-(1:100),]
frq21$nfreq <- apply(frq21[,8:107],1,sum); frq21$nfreq[frq21$nfreq > 1000] <- 1000
frq18$nfreq <- apply(frq18[,8:107],1,sum); frq18$nfreq[frq18$nfreq > 1000] <- 1000
frq15$nfreq <- apply(frq15[,8:107],1,sum); frq15$nfreq[frq15$nfreq > 1000] <- 1000
frq12$nfreq <- apply(frq12[,8:107],1,sum); frq12$nfreq[frq12$nfreq > 1000] <- 1000
nf21 <- aggregate(nfreq ~ year, data = frq21, sum)
nf18 <- aggregate(nfreq ~ year, data = frq18, sum)
nf15 <- aggregate(nfreq ~ year, data = frq15, sum)
nf12 <- aggregate(nfreq ~ year, data = frq12, sum)

windows()
plot(nf21$year, nf21$nfreq, xlab="Year",  pch=3, col=3, xlim = c(1960, 2020), ylim = c(0, max(nf12$nfreq)))
points(nf18$year, nf18$nfreq, pch = 2, col=2)
points(nf15$year, nf15$nfreq, pch = 1, col=1)
#points(nf12$year, nf12$nfreq, pch = 1, col=1)
legend("topleft", legend=c(2015, 2018, 2021), col = 1:3, pch=1:3)
savePlot(paste0(skj22_dir,"/figures/","LF samples.png"), type = "png")

xfrq21 <- filter(frq21, fishery %in% c(15, 16))
xfrq18 <- filter(frq18, fishery %in% c(13,14))
xfrq15 <- filter(frq15, fishery %in% 9:11)
xfrq12 <- filter(frq12, fishery %in% 137:138)
nf21 <- aggregate(nfreq ~ year, data = xfrq21, sum)
nf18 <- aggregate(nfreq ~ year, data = xfrq18, sum)
nf15 <- aggregate(nfreq ~ year, data = xfrq15, sum)
nf12 <- aggregate(nfreq ~ year, data = xfrq12, sum)
windows() 
plot(nf21$year, nf21$nfreq, xlab="Year", ylab = "# LF samples", 
     pch=3, col=3, xlim = c(1960, 2020), main = "Troll fisheries")
points(nf18$year, nf18$nfreq, pch = 2, col=2)
points(nf15$year, nf15$nfreq, pch = 1, col=1)
#points(nf12$year, nf12$nfreq, pch = 1, col=1)
legend("topleft", legend=c(2015, 2018, 2021), col = 1:3, pch=1:3)
savePlot(paste0(skj22_dir,"/figures/","LF samples troll.png"), type = "png")

xfrq21 <- filter(frq21, !fishery %in% c(15:18))
xfrq18 <- filter(frq18, !fishery %in% c(13:16))
xfrq15 <- filter(frq15, !fishery %in% 9:14)
xfrq12 <- filter(frq12, !fishery %in% 137:140)
nf21 <- aggregate(nfreq ~ year, data = xfrq21, sum)
nf18 <- aggregate(nfreq ~ year, data = xfrq18, sum)
nf15 <- aggregate(nfreq ~ year, data = xfrq15, sum)
nf12 <- aggregate(nfreq ~ year, data = xfrq12, sum)

windows(height=12, width=10); par(mfrow = c(2,1))
plot(nf21$year, nf21$nfreq, xlab="Year", ylab = "# LF samples", log="y", ylim = c(50, max(nf12$nfreq)),
     pch=3, col=3, xlim = c(1960, 2020), main = "Longline fisheries (log scale)")
points(nf18$year, nf18$nfreq, pch = 2, col=2)
points(nf15$year, nf15$nfreq, pch = 1, col=1)
#points(nf12$year, nf12$nfreq, pch = 1, col=1)
legend("bottomright", legend=c( 2015, 2018, 2021), col = 1:3, pch=1:3)
#savePlot(paste0(skj22_dir,"/figures/","LF samples longline log.png"), type = "png")

#windows()
plot(nf21$year, nf21$nfreq, xlab="Year", ylab = "# LF samples", ylim = c(0, max(nf12$nfreq)),
     pch=3, col=3, xlim = c(1960, 2020), main = "Longline fisheries")
points(nf18$year, nf18$nfreq, pch = 2, col=2)
points(nf15$year, nf15$nfreq, pch = 1, col=1)
#points(nf12$year, nf12$nfreq, pch = 1, col=1)
legend("topleft", legend=c(2015, 2018, 2021), col = 1:3, pch=1:3)
savePlot(paste0(skj22_dir,"/figures/","LF samples longline.png"), type = "png")

# Apply ff49 scalar
par21 <- read.par(paste0(diagn_dir, "13.par"))$ffl[,49]
par12 <- read.par(paste0(skj22_dir, "/2012_skj_basecase/08.par"))$ffl[,49]
par15 <- read.par(paste0(skj22_dir, "/2015_skj_basecase/skj2015/12.par"))$ffl[,49]
par18 <- read.par(paste0(skj22_dir, "/SPA2018Web/GridML_GeoCPUE_GrowthEstim_SizeWgt50_Steep0.8_M0.3_diag/13.par"))$ffl[,49]
frq21$nfreq2 <- frq21$nfreq / par21[frq21$fishery]
frq18$nfreq2 <- frq18$nfreq / par18[frq18$fishery]
frq15$nfreq2 <- frq15$nfreq / par15[frq15$fishery]
frq12$nfreq2 <- frq12$nfreq / par12[frq12$fishery]

xfrq21 <- filter(frq21, !fishery %in% c(15:18))
xfrq18 <- filter(frq18, !fishery %in% c(13:16))
xfrq15 <- filter(frq15, !fishery %in% 9:14)
xfrq12 <- filter(frq12, !fishery %in% 137:140)
nf21 <- aggregate(nfreq2 ~ year, data = xfrq21, sum)
nf18 <- aggregate(nfreq2 ~ year, data = xfrq18, sum)
nf15 <- aggregate(nfreq2 ~ year, data = xfrq15, sum)
nf12 <- aggregate(nfreq2 ~ year, data = xfrq12, sum)

windows(height=12, width=10); par(mfrow = c(2,1))
plot(nf21$year, nf21$nfreq2, xlab="Year", ylab = "# LF samples", log="y", ylim = c(1, max(nf21$nfreq2)),
     pch=3, col=3, xlim = c(1960, 2020), main = "Longline fisheries (log scale)")
points(nf18$year, nf18$nfreq2, pch = 2, col=2)
points(nf15$year, nf15$nfreq2, pch = 1, col=1)
#points(nf12$year, nf12$nfreq2, pch = 1, col=1)
legend("bottomright", legend=c(2015, 2018, 2021), col = 1:3, pch=1:3)
#savePlot(paste0(skj22_dir,"/figures/","LF samples longline log.png"), type = "png")

#windows()
plot(nf21$year, nf21$nfreq2, xlab="Year", ylab = "# LF samples", ylim = c(0, max(nf21$nfreq2)),
     pch=3, col=3, xlim = c(1960, 2020), main = "Longline fisheries")
points(nf18$year, nf18$nfreq2, pch = 2, col=2)
points(nf15$year, nf15$nfreq2, pch = 1, col=1)
#points(nf12$year, nf12$nfreq2, pch = 1, col=1)
legend("topleft", legend=c(2015, 2018, 2021), col = 1:3, pch=1:3)
savePlot(paste0(skj22_dir,"/figures/","LF samples longline adjusted.png"), type = "png")


##########----------------------------------------------------------------
#Plot all the runs
##########----------------------------------------------------------------
# Load all the runs
setwd(skj22_dir)
#testrep1 <- testlens <- list()
alldirs <- dir(path="./test/", full.names=TRUE)[-c(1)]
nnm <- dir(path="./test/")[-c(1)]
for (i in 58:length(alldirs)) {
  dd <- alldirs[i]
  # print(dd)
  # print(list.files(path=dd, pattern="par.rep"))
  a <- list.files(path=dd, pattern="par.rep")
  a <- a[-grep("q0",a)]
  a <- a[length(a)]
  if(file.exists(paste0(dd, "/", a)) & length(a) != 0) { 
    res <- read.rep(paste0(dd, "/", a))
    testrep1[[i]] <- res
  }
  if(file.exists(paste0(dd, "/length.fit"))) { 
    lens <- read.fit(paste0(dd, "/length.fit"))
    testlens[[i]] <- lens
  }
  if(file.exists(paste0(dd, "/test_plot_output"))) {
    testrep1[[i]]$lenLiks <- read.MFCLLikelihood(paste0(dd, "/test_plot_output"))
  }
}
for (i in 1:length(alldirs)) testrep1[[i]]$dirname <- nnm[i]


for(i in 1:length(testrep1)){ 
  if(i %in% c(10, 19, 28, 37, 46, 55)){
    savePlot(paste0(skj22_dir,"/figures/allruns2_recruitment_",i-1,".png"), type = "png")
  }
  if(i %in% c(1, 10, 19, 28, 37, 46, 55)) {
    windows(10,10); par(mfrow=c(3,3))}
  if(length(testrep1[[i]]) > 0) {
    yrs <- testrep1[[i]]$yrs
    doyrs <- (max(yrs) - 25):max(yrs)
    gety <- (length(yrs)-100): length(yrs)
    yl <- max(apply(testrep1[[i]]$Recruitment,1,sum)[-(1:160)])
    plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = nnm[i])
    lines(testrep1[[i]]$yrs[gety], apply(testrep1[[i]]$Recruitment,1,sum)[gety], col = i, lwd=1, lty=1)
  }}
savePlot(paste0(skj22_dir,"/figures/allruns2_recruitment_",i,".png"), type = "png")

for(i in 1:length(testrep1)){ 
  if(i %in% c(10, 19, 28, 37, 46, 55)){
    savePlot(paste0(skj22_dir,"/figures/allruns2_recruit2_",i-1,".png"), type = "png")
  }
  if(i %in% c(1, 10, 19, 28, 37, 46, 55)) {
    windows(10,10); par(mfrow=c(3,3))}
  if(length(testrep1[[i]]) > 0) {
    yrs <- testrep1[[i]]$yrs
    doyrs <- testrep1[[i]]$yrs
    gety <- 1:length(yrs)
    yl <- max(apply(testrep1[[i]]$Recruitment,1,sum))
    plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = nnm[i])
    lines(testrep1[[i]]$yrs[gety], apply(testrep1[[i]]$Recruitment,1,sum)[gety], col = i, lwd=1, lty=1)
  }}
savePlot(paste0(skj22_dir,"/figures/allruns2_recruit2_",i,".png"), type = "png")

for(i in 1:length(testrep1)){ 
  if(i %in% c(10, 19, 28, 37, 46, 55)){
    savePlot(paste0(skj22_dir,"/figures/allruns2_F_",i-1,".png"), type = "png")
  }
  if(i %in% c(1, 10, 19, 28, 37, 46, 55)) {
    windows(10,10); par(mfrow=c(3,3))}
  if(length(testrep1[[i]]) > 0) {
    yrs <- testrep1[[i]]$yrs
    doyrs <- (max(yrs) - 25):max(yrs)
    gety <- (length(yrs)-100): length(yrs)
    yl <- max(testrep1[[i]]$FatYrAgeReg[,,1])
    plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="F", type = "n", main = nnm[i])
    for(ag in 1:48) lines(testrep1[[i]]$yrs[gety], testrep1[[i]]$FatYrAgeReg[,ag,1][gety], col = i, lwd=1, lty=1)
}}
savePlot(paste0(skj22_dir,"/figures/allruns2_F_",i,".png"), type = "png")

for(i in 1:length(testrep1)){ 
  if(i %in% c(10, 19, 28, 37, 46, 55)){
    savePlot(paste0(skj22_dir,"/figures/allruns2_R1_depletion_",i-1,".png"), type = "png")
  }
  if(i %in% c(1, 10, 19, 28, 37, 46, 55)) {
    windows(10,10); par(mfrow=c(3,3))}
  if(length(testrep1[[i]]) > 0) {
    yrs <- testrep1[[i]]$yrs
    doyrs <- (max(yrs) - 25):max(yrs)
    gety <- (length(yrs)-100): length(yrs)
    yl <- 1
    plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Depletion", type = "n", main = nnm[i])
    lines(testrep1[[i]]$yrs[gety], (testrep1[[i]]$AdultBiomass[,1]/testrep1[[i]]$AdultBiomass.nofish[,1])
          [gety], col = i, lwd=1, lty=1)
}}
savePlot(paste0(skj22_dir,"/figures/allruns2_R1_depletion_",i,".png"), type = "png")

for(i in 1:length(testrep1)){ 
  if(i %in% c(10, 19, 28, 37, 46, 55)){
    savePlot(paste0(skj22_dir,"/figures/allruns2_SSB_",i-1,".png"), type = "png")
  }
  if(i %in% c(1, 10, 19, 28, 37, 46, 55)) {
    windows(10,10); par(mfrow=c(3,3))}
  if(length(testrep1[[i]]) > 0) {
    yrs <- testrep1[[i]]$yrs
    doyrs <- (max(yrs) - 25):max(yrs)
    gety <- (length(yrs)-100): length(yrs)
    yl <- max(apply(testrep1[[i]]$AdultBiomass,1,sum))
    plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="SSB", type = "n", main = nnm[i])
    lines(testrep1[[i]]$yrs[gety], (apply(testrep1[[i]]$AdultBiomass,1,sum))[gety], col = i, lwd=1, lty=1)
}}
savePlot(paste0(skj22_dir,"/figures/allruns2_SSB_",i,".png"), type = "png")
graphics.off()
 
library(mgcv)

windows(7,10); par(mfrow=c(2,1))
yrs <- doyrs <- testrep1[[3]]$yrs
gety <- 1:length(yrs)
tmp <- data.frame(yrs=testrep1[[3]]$yrs)
tmp$rec <-  apply(testrep1[[3]]$Recruitment,1,sum)
yl <- max(tmp$rec)
plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = "Recruitment in diagnostic model")
lines(tmp$yrs, tmp$rec, col = 1, lwd=1, lty=1)
gammod <- gam(rec ~ s(yrs, k=30), data=tmp)
lines(tmp$yrs, predict.gam(gammod, newdata= data.frame(yrs=tmp$yrs), type = "response"), col=2)
yrs <- testrep1[[3]]$yrs
doyrs <- testrep1[[3]]$yrs
gety <- 1:length(yrs)
tmp <- data.frame(yrs=testrep1[[3]]$yrs)
tmp$yrs2 <- floor(tmp$yrs)
tmp$rec <-  apply(testrep1[[3]]$Recruitment,1,sum)
tmp2 <- aggregate(rec ~ yrs2, data = tmp, FUN=sum)
yl <- max(tmp2$rec)
plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = "Annual Recruitment in diagnostic model")
lines(tmp2$yrs2, tmp2$rec, col = 1, lwd=1, lty=1)
gammod <- gam(rec ~ s(yrs2, k=10), data=tmp2)
lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
savePlot(paste0(skj22_dir,"/figures/diag_recruitment.png"), type = "png")
#savePlot(paste0(skj22_dir,"/figures/diag_recruitment_annual.png"), type = "png")

windows(7,10); par(mfrow=c(2,1))
plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = "Annual Recruitment in diagnostic model")
lines(tmp2$yrs2, tmp2$rec, col = 1, lwd=1, lty=1)
lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
tmp2$dec <- 5*floor(tmp2$yrs2/5)
tmp3 <- aggregate(rec ~ dec, data=tmp2, FUN=mean)
barplot(tmp3$rec, names.arg=tmp3$dec)
savePlot(paste0(skj22_dir,"/figures/diag_recruit_bplt_2021.png"), type = "png")



res2012 <- read.rep("~/../../OneDrive - NIWA/International/WCPFC/skj_2021/2012_skj_basecase/plot-08.par.rep")
res2015 <- read.rep("~/../../OneDrive - NIWA/International/WCPFC/skj_2021/2015_skj_basecase/skj2015/plot-12.par.rep")
res2018 <- read.rep("~/../../OneDrive - NIWA/International/WCPFC/skj_2021/SPA2018Web/GridML_GeoCPUE_GrowthEstim_SizeWgt50_Steep0.8_M0.3_diag/plot-13.par.rep")
res2021 <- read.rep("~/../../OneDrive - NIWA/International/WCPFC/skj_2021/stepwiseskj21/8_S2M1D2R1G1_diagcase/plot-13.par.rep")

windows(7,10); par(mfrow=c(2,1))
yrs <- doyrs <- res2018$yrs
gety <- 1:length(yrs)
tmp <- data.frame(yrs=res2018$yrs)
tmp$rec <-  apply(res2018$Recruitment,1,sum)
yl <- max(tmp$rec)
plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = "Recruitment in diagnostic model")
lines(tmp$yrs, tmp$rec, col = 1, lwd=1, lty=1)
gammod <- gam(rec ~ s(yrs, k=30), data=tmp)
lines(tmp$yrs, predict.gam(gammod, newdata= data.frame(yrs=tmp$yrs), type = "response"), col=2)
yrs <- res2018$yrs
doyrs <- res2018$yrs
gety <- 1:length(yrs)
tmp <- data.frame(yrs=res2018$yrs)
tmp$yrs2 <- floor(tmp$yrs)
tmp$rec <-  apply(res2018$Recruitment,1,sum)
tmp2 <- aggregate(rec ~ yrs2, data = tmp, FUN=sum)
yl <- max(tmp2$rec)
plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = "Annual Recruitment in diagnostic model")
lines(tmp2$yrs2, tmp2$rec, col = 1, lwd=1, lty=1)
gammod <- gam(rec ~ s(yrs2, k=30), data=tmp2)
lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
savePlot(paste0(skj22_dir,"/figures/diag_recruitment_2018.png"), type = "png")

windows(7,10); par(mfrow=c(2,1))
plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = "Annual Recruitment in diagnostic model")
lines(tmp2$yrs2, tmp2$rec, col = 1, lwd=1, lty=1)
lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
tmp2$dec <- 5*floor(tmp2$yrs2/5)
tmp3 <- aggregate(rec ~ dec, data=tmp2, FUN=mean)
barplot(tmp3$rec, names.arg=tmp3$dec)
savePlot(paste0(skj22_dir,"/figures/diag_recruit_bplt_2018.png"), type = "png")

for (resn in c(2012, 2015, 2018, 2021)){
  ddd <- get(paste0("res", resn))
  yrs <- doyrs <- ddd$yrs
  gety <- 1:length(yrs)
  tmp <- data.frame(yrs=ddd$yrs)
  tmp$rec <-  apply(ddd$Recruitment,1,sum)
  yl <- max(tmp$rec)
  plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = paste("Recruitment in diagnostic model", resn))
  lines(tmp$yrs, tmp$rec, col = 1, lwd=1, lty=1)
  gammod <- gam(rec ~ s(yrs, k=30), data=tmp)
  lines(tmp$yrs, predict.gam(gammod, newdata= data.frame(yrs=tmp$yrs), type = "response"), col=2)
}
windows(10,10); par(mfrow=c(2,2))
for (resn in c(2012, 2015, 2018, 2021)){
  ddd <- get(paste0("res", resn))
  yrs <- doyrs <- ddd$yrs
  gety <- 1:length(yrs)
  tmp <- data.frame(yrs=ddd$yrs)
  tmp$yrs2 <- floor(tmp$yrs)
  tmp$rec <-  apply(ddd$Recruitment,1,sum)
  tmp2 <- aggregate(rec ~ yrs2, data = tmp, FUN=sum)
  yl <- max(tmp2$rec)
  plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = paste("Annual Rec diagnostic model", resn))
  lines(tmp2$yrs2, tmp2$rec, col = 1, lwd=1, lty=1)
  gammod <- gam(rec ~ s(yrs2, k=10), data=tmp2)
  lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
}
savePlot(paste0(skj22_dir,"/figures/diag_recruitment_annual_2012_2021.png"), type = "png")

frq2021 <- read.frq(paste0(diagn_dir, "skj.frq"))
frq2012 <- read.frq(paste0(skj22_dir, "/2012_skj_basecase/skj.frq"))
frq2015 <- read.frq(paste0(skj22_dir, "/2015_skj_basecase/skj2015/skj.frq"))
frq2018 <- read.frq(paste0(skj22_dir, "/SPA2018Web/GridML_GeoCPUE_GrowthEstim_SizeWgt50_Steep0.8_M0.3_diag/skj.frq"))
frq2021$fish$ctype

windows(10,10); 
par(mfrow=c(2,2), mar = c(4,4,4,5))
for (resn in c(2012, 2015, 2018, 2021)){
  ddd <- get(paste0("res", resn))
  frq <- get(paste0("frq", resn))
  yrs <- doyrs <- ddd$yrs
  gety <- 1:length(yrs)
  tmp <- data.frame(yrs=ddd$yrs)
  tmp$yrs2 <- floor(tmp$yrs)
  tmp$rec <-  apply(ddd$Recruitment,1,sum)
  tmp2 <- aggregate(rec ~ yrs2, data = tmp, FUN=sum, na.rm=T)
  yl <- max(tmp2$rec)
  plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = paste("Diagnostic model", resn))
  lines(tmp2$yrs2, tmp2$rec, col = 1, lwd=1, lty=1)
  # gammod <- gam(rec ~ s(yrs2, k=10), data=tmp2)
  # lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
  par(new=TRUE)
  a <- data.frame(frq$mat[,1:7])
  a <- aggregate(catch ~ year, data = a[!frq$fish$ctype[a$fishery],], FUN=sum)
  plot(a$year, a$catch, col=2, pch=19, axes=FALSE, xlab="", ylab="",ylim = c(0, max(a$catch)))
  axis(side=4, col = 2)
  mtext("Longline catch (number)", side=4, col=2, line=3)
}
savePlot(paste0(skj22_dir,"/figures/diag_rec_catch_2012_2021.png"), type = "png")

windows(10,10); 
par(mfrow=c(2,2), mar = c(4,4,4,5))
for (resn in c(2012, 2015, 2018, 2021)){
  ddd <- get(paste0("res", resn))
  frq <- get(paste0("frq", resn))
  yrs <- doyrs <- ddd$yrs
  gety <- 1:length(yrs)
  tmp <- data.frame(yrs=ddd$yrs)
  tmp$yrs2 <- floor(tmp$yrs)
  tmp$rec <-  apply(ddd$Recruitment,1,sum)
  tmp2 <- aggregate(rec ~ yrs2, data = tmp, FUN=sum, na.rm=T)
  yl <- max(tmp$rec)
  plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = paste("Diagnostic model", resn))
  lines(tmp$yrs, tmp$rec, col = 1, lwd=1, lty=1)
  # gammod <- gam(rec ~ s(yrs2, k=10), data=tmp2)
  # lines(tmp2$yrs2, predict.gam(gammod, newdata= data.frame(yrs2=tmp2$yrs2), type = "response"), col=2)
  par(new=TRUE)
  a <- data.frame(frq$mat[,1:7])
  a <- aggregate(catch ~ year, data = a[!frq$fish$ctype[a$fishery],], FUN=sum)
  plot(a$year, a$catch, col=2, pch=19, axes=FALSE, xlab="", ylab="",ylim = c(0, max(a$catch)))
  axis(side=4, col = 2)
  mtext("Longline catch (number)", side=4, col=2, line=3)
}
savePlot(paste0(skj22_dir,"/figures/diag_recq_catch_2012_2021.png"), type = "png")

plot(apply(res2012$ObsCatch[!frq2012$fish$ctype,], 1, sum, na.rm=T))
apply(res2012$ObsCatch[!frq2012$fish$ctype,], 1, sum, na.rm=T)

a <- data.frame(frq2021$mat[,1:7])
a <- aggregate(catch ~ year, data = a[!frq2021$fish$ctype[a$fishery],], FUN=sum)
plot(a[,1], a[,2])

names(res2021)
head(res2021$ObsCatch)
apply(res2021$ObsCatch, 2, sum, na.rm=T)[1:240]

enso <- read.csv("enso.csv")
windows(height=8, width=10)
plot(enso$yr2, enso$enso, type = "l", ylim=c(-3,3), xlab="Year", ylab="ENSO")
ddd <- res2021
yrs <- doyrs <- ddd$yrs
gety <- 1:length(yrs)
tmp <- data.frame(yrs=ddd$yrs)
tmp$yrs2 <- floor(tmp$yrs)
tmp$rec <-  apply(ddd$Recruitment,1,sum)
tmp2 <- aggregate(rec ~ yrs2, data = tmp, FUN=sum, na.rm=T)
par(new=TRUE)
#plot(yrs,tmp$rec/mean(tmp$rec)-1, col=2, axes=F, type = "l", xlim=c(1979, 2023), ylim=c(-2,2))
plot(tmp2$yrs,tmp2$rec/mean(tmp2$rec)-1, col=2, axes=F, type = "b", xlim=c(1979, 2023), ylim=c(-.8,.8), lwd=2, xlab="", ylab="")
axis(side=4, col=2)
abline(h=0, lty=2)
legend("topright", legend = c("ENSO", "Recruitment"), col=1:2, lwd=c(1,2))
savePlot(paste0(skj22_dir,"/figures/ENSO versus recruitment.png"), type = "png")

enso_annual <- aggregate(enso ~ floor(yr2), data = enso, FUN=mean, na.rm=T)
names(enso_annual)[1] <- "yr2"
cor(enso_annual$enso[enso_annual$yr2 %in% c(1979:2021)], tmp2$rec[tmp2$yrs2 %in% c(1979:2021)])
cor(enso_annual$enso[enso_annual$yr2 %in% c(1979:2013)], tmp2$rec[tmp2$yrs2 %in% c(1979:2013)])
enso_annual$rec <- NA
enso_annual$rec[enso_annual$yr2 %in% c(1979:2021)] <- tmp2$rec[tmp2$yrs2 %in% c(1979:2021)]

library(mgcv)
mod <-  gam(rec ~ enso, data=enso_annual[enso_annual$yr2 %in% c(1979:2013),])
summary(mod) # Relationship is not significant
# correlation is -0.0532879

# ----------------------------------------------------------

windows()
testrep1[[27]]

barplot(total_length_fish(run27liks)-total_length_fish(diag_liks), names.arg=1:25, main = "run27_diag_last24recs")
savePlot(paste0(skj22_dir,"/figures/","run27_sizelik barplot.png"), type = "png")
barplot(total_length_fish(run48liks)-total_length_fish(diag_liks), names.arg=1:25, main = "run48_fixed_biomass_scale")
savePlot(paste0(skj22_dir,"/figures/","run48_sizelik barplot.png"), type = "png")
barplot(total_length_fish(testrep1[[29]]$lenLiks)-total_length_fish(diag_liks), names.arg=1:25, main = "runfix200_diag_last24recs")
savePlot(paste0(skj22_dir,"/figures/","runfix200_sizelik barplot.png"), type = "png")

diag_lenfit <- read.fit(paste0(source_basedir, model=modnames[8], '/length.fit'))
for (i in 1:25) {
  a <- diag_lenfit$dates[[i]]
  diag_lenfit$dates[[i]]$yq <- a$Year + (a$Month-0.5)/12
}
windows(15,15); par(mfrow=c(5,5), mar = c(1,2,4,1))
for (i in 1:25) {
  plot(diag_lenfit$dates[[i]]$yq,length_fish(run27liks)[[i]]-length_fish(diag_liks)[[i]], xlim = c(1960, 2020), main = i)
}
savePlot(paste0(skj22_dir,"/figures/","run27_sizelik plot.png"), type = "png")
windows(15,15); par(mfrow=c(5,5), mar = c(1,2,4,1))
for (i in 1:25) {
  plot(diag_lenfit$dates[[i]]$yq,length_fish(run48liks)[[i]]-length_fish(diag_liks)[[i]], xlim = c(1960, 2020), main = i)
}
savePlot(paste0(skj22_dir,"/figures/","run48_sizelik plot.png"), type = "png")
windows(15,15); par(mfrow=c(5,5), mar = c(1,2,4,1))
for (i in 1:25) plot(diag_lenfit$dates[[i]]$yq,length_fish(testrep1[[29]]$lenLiks)[[i]]-length_fish(diag_liks)[[i]], xlim = c(1960, 2020), main = i)
savePlot(paste0(skj22_dir,"/figures/","runfix200_sizelik plot.png"), type = "png")




allretro61 <- allretro62 <- allretro71 <- allretro72 <- list()
retrodir <- "./retro/nan_retro"
i=1
for(yy in 2009:2021) {
  repfile6 <- paste0(retrodir, "/", "plot-06_", yy, ".par.rep")
#  repfile7 <- paste0(retrodir, "/", "plot-07_", yy, ".par.rep")
  res6 <- read.rep(repfile6)
  allretro61[[i]] <- res6
#  allretro62[[i]] <- read.MFCLRep(repfile6)
  i <- i+1
}
rm(allretro1, allretro2,allretro71,allretro72)

windows(10,10); par(mfrow=c(3,4))
for(i in 1:11){ 
  yrs <- allretro61[[i]]$yrs
  doyrs <- (max(yrs) - 25):max(yrs)
  gety <- (length(yrs)-100): length(yrs)
  yl <- max(apply(allretro61[[i]]$Recruitment,1,sum)[-(1:160)])
  plot(doyrs, doyrs, ylim = c(0, yl), xlab = "Year", ylab="Recruitment", type = "n", main = (2009:2021)[i])
  lines(allretro61[[i]]$yrs[gety], apply(allretro61[[i]]$Recruitment,1,sum)[gety], col = i, lwd=1, lty=1)
}
savePlot(paste0(skj22_dir,"/figures/","retro_recruitment.png"), type = "png")

proj_dir = "D:/HOME/SAP/2023_SC/skj_investigation/"

# load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(FLR4MFCL)

# read in data for diagnostic model
diagnostic_rep_files = paste0(retrodir, "/plot-06_",2021:2009,".par.rep")
rep_list = lapply(diagnostic_rep_files,read.MFCLRep)

n_dt.list = as.list(rep(NA,length(rep_list)))
for(i in 1:(length(n_dt.list)))
{
  n_dt.list[[i]] = as.data.table(popN(rep_list[[i]])) %>%
    .[,peel:=i-1] %>%
    .[,model:="S2M1D2R1G1"] %>%
    .[,age:=as.numeric(age)/4] %>%
    .[,yy:=as.numeric(year)] %>%
    .[,qq:=as.numeric(season)] %>%
    .[,area:=as.numeric(area)] %>%
    .[,yyqq:=yy+(qq-1)/4] %>%
    .[,.(model,peel,yy,qq,yyqq,area,age,value)]
  
}

n_dt = rbindlist(n_dt.list)

windows()
n_dt %>%
  .[age==0.25] %>%
  .[,ts_plot:=yyqq-max(yyqq),by=peel] %>%
  .[,n_plot:=value/tail(value,n=1),by=.(area,peel)] %>%
  .[ts_plot>=-12] %>%
  ggplot() +
  xlab("Year") +
  ylab("N") +
  facet_wrap(~area,scales="free") +
  geom_path(aes(x=ts_plot,y=n_plot,color=peel,group=peel)) +
  geom_smooth(aes(x=ts_plot,y=n_plot),span=0.1,color="black") +
  viridis::scale_color_viridis("Peel",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
        panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
        panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
        strip.background =element_rect(fill="white"),
        legend.key = element_rect(fill = "white"))
savePlot(paste0(skj22_dir,"/figures/","retro_rec_byreg_smoothed_2021.png"), type = "png")
