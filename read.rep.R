read.rep <- function (rep.file, verbose = FALSE, DEBUG = FALSE) 
{
  cat("Starting read.rep ; ")
  if (0) {
    datfromstr <- function(datstring, transpose = FALSE) {
      out <- if (length(datstring) > 1) {
        datstring %>% sapply("trimws") %>% strsplit(split = "[[:blank:]]+") %>% 
          lapply(., "as.numeric") %>% {
            if (all(sapply(., length) == max(sapply(., 
                                                    length)))) 
              do.call("rbind", .)
            else matrix(unlist(lapply(., "[", 1:max(sapply(., 
                                                           length)))), nrow = length(.), byrow = TRUE)
          } %>% `colnames<-`(NULL) %>% `rownames<-`(NULL) %>% 
          {
            if (transpose) 
              t(.)
            else .
          }
      }
      else {
        datstring %>% trimws() %>% strsplit(split = "[[:blank:]]+") %>% 
          1[[]] %>% as.numeric()
      }
      return(out)
    }
  }
  # charVec2numVec <- function(x) {
  #   if (length(x) > 0) 
  #     x %>% trimws() %>% strsplit(., split = " +") %>% 
  #     1[[]] %>% as.numeric()
  #   else NA
  # }
  charVec2numVec <- function(x) {
    if (length(x) > 0) 
      x %>% trimws() %>% strsplit(split = " +") %>% unlist() %>% as.numeric()
    else NA
  }
  # getNumVector <- function(keyword, a, offset = 1) {
  #   pos1 <- grep(pattern = keyword, x = a)
  #   if (length(pos1) == 1) {
  #     a[pos1 + offset] %>% trimws() %>% strsplit(., split = "[[:blank:]]+") %>% 
  #       1[[]] %>% as.numeric()
  #   }
  #   else if (length(pos1) == 0) {
  #     NA
  #   }
  #   else {
  #     stop("grep(", keyword, ")=", pos1)
  #   }
  # }
  getNumVector <- function(keyword, a, offset = 1) {
    pos1 <- grep(pattern = keyword, x = a)
    if (length(pos1) == 1) {
      a[pos1 + offset] %>% trimws() %>% strsplit(., split = "[[:blank:]]+") %>% 
        unlist() %>% as.numeric()
    }
    else if (length(pos1) == 0) {
      NA
    }
    else {
      stop("grep(", keyword, ")=", pos1)
    }
  }
  a <- readLines(rep.file)
  viewerVer <- strsplit(a[2], split = " +")[[1]][2]
  viewerVer.num <- as.numeric(viewerVer)
  frqfilename <- strsplit(a[3], split = " +")[[1]][5]
  inputParfilename <- strsplit(a[4], split = " +")[[1]][6]
  outputParfilename <- strsplit(a[6], split = " +")[[1]][6]
  MFCL.version.string <- strsplit(a[5], split = " +")[[1]][5]
  MFCL.version <- list()
  MFCL.version$Major <- as.numeric(strsplit(MFCL.version.string, 
                                            ".")[[1]][1])
  MFCL.version$Minor <- as.numeric(strsplit(MFCL.version.string, 
                                            ".")[[1]][2])
  MFCL.version$Build <- as.numeric(strsplit(MFCL.version.string, 
                                            ".")[[1]][3])
  MFCL.version$Revision <- as.numeric(strsplit(MFCL.version.string, 
                                               ".")[[1]][4])
  pos1 <- grep("# Number of time periods", a)
  nTimes <- datfromstr(a[pos1 + 1])
  pos1 <- grep("# Year 1", a)
  Year1 <- datfromstr(a[pos1 + 1])
  pos1 <- grep("# Number of regions", a)
  nReg <- datfromstr(a[pos1 + 1])
  if (viewerVer.num >= 3) {
    pos1 <- grep("# Number of species", a)
    nSp <- 1
    nSp <- datfromstr(a[pos1 + 1])
  }
  else {
    nSp <- 1
  }
  if (viewerVer.num >= 4) {
    pos1 <- grep("# Multi species pointer", a)
    spPtr <- NA
    spPtr <- datfromstr(a[pos1 + 1])
  }
  else {
    spPtr <- NA
  }
  if (viewerVer.num >= 4) {
    pos1 <- grep("# Species sex pointer", a)
    spSexPtr <- NA
    spSexPtr <- datfromstr(a[pos1 + 1])
  }
  else {
    spSexPtr <- NA
  }
  pos1 <- grep("# Number of age classes", a)
  nAges <- datfromstr(a[pos1 + 1])
  if (viewerVer.num >= 3) {
    pos1 <- grep("# Regions species pointer", a)
    regSpPtr <- NA
    regSpPtr <- datfromstr(a[pos1 + 1])
  }
  else {
    regSpPtr <- NA
  }
  pos1 <- grep("# Number of recruitments per year", a)
  nRecs.yr <- datfromstr(a[pos1 + 1])
  pos1 <- grep("# Number of fisheries", a)
  nFisheries <- datfromstr(a[pos1 + 1])
  if (verbose) 
    cat("L79 ; ")
  if (viewerVer.num >= 3) {
    pos1 <- grep("# Fishery selectivity seasons", a)
    SelexSeasons <- datfromstr(a[pos1 + 1])
  }
  else {
    SelexSeasons <- rep(1, nFisheries)
  }
  if (viewerVer.num >= 3) {
    pos1 <- grep("# Fishery selectivity time-blocks", a)
    SelexTblocks <- datfromstr(a[pos1 + 1])
  }
  else {
    SelexTblocks <- rep(1, nFisheries)
  }
  pos1 <- grep("# Number of realizations per fishery", a)
  nRlz.fsh <- datfromstr(a[pos1 + 1])
  if (verbose) 
    cat("L83 ; ")
  pos1 <- grep("# Region for each fishery", a)
  Region.fsh <- datfromstr(a[pos1 + 1])
  pos1 <- grep("# Time of each realization by fishery", a)
  Rlz.t.fsh <- matrix(nrow = nFisheries, ncol = max(nRlz.fsh))
  for (i in 1:nFisheries) {
    Rlz.t.fsh[i, 1:nRlz.fsh[i]] <- datfromstr(a[pos1 + i])
  }
  if (verbose) 
    cat("L89 ; ")
  if (DEBUG) 
    browser()
  pos1 <- grep("# Mean lengths at age", a)
  mean.LatAge <- datfromstr(a[pos1 + 1])
  pos1 <- grep("# SD of length at age", a)
  sd.LatAge <- datfromstr(a[pos1 + 1])
  pos1 <- grep("# Mean weights at age", a)
  mean.WatAge <- datfromstr(a[pos1 + 1])
  if (nSp > 1) {
    colnames(mean.LatAge) <- colnames(sd.LatAge) <- colnames(mean.WatAge) <- paste0(1:nAges[1])
    rownames(mean.LatAge) <- rownames(sd.LatAge) <- rownames(mean.WatAge) <- c("Male", 
                                                                               "Female")
  }
  else {
    mean.LatAge <- as.data.frame(t(mean.LatAge))
    sd.LatAge <- as.data.frame(t(sd.LatAge))
    mean.WatAge <- as.data.frame(t(mean.WatAge))
    colnames(mean.LatAge) <- colnames(sd.LatAge) <- colnames(mean.WatAge) <- paste0(1:nAges[1])
    rownames(mean.LatAge) <- rownames(sd.LatAge) <- rownames(mean.WatAge) <- c("Both")
  }
  if (verbose) 
    cat("L103 ; ")
  if (DEBUG) 
    browser()
  pos1 <- grep("# Natural mortality at age", a)
  MatAge <- datfromstr(a[pos1 + 1:nSp])
  if (verbose) 
    cat("L107 ; ")
  pos1 <- grep("# Selectivity by age class ", a)
  SelAtAge <- datfromstr(a[(pos1 + 1):(pos1 + sum(SelexTblocks))])
  pos1 <- grep("# Implicit catchability by realization ", a)
  qAtAge <- matrix(nrow = nFisheries, ncol = max(nRlz.fsh))
  for (i in 1:nFisheries) {
    qAtAge[i, 1:nRlz.fsh[i]] <- datfromstr(a[pos1 + i])
  }
  if (verbose) 
    cat("L113 ; ")
  if (DEBUG) 
    browser()
  # pos1 <- grep("# Catchability\\+effort dev\\. by realization ", 
  #              a)
  # qEdevAtAge <- matrix(nrow = nFisheries, ncol = max(nRlz.fsh))
  # for (i in 1:nFisheries) {
  #   qEdevAtAge[i, 1:nRlz.fsh[i]] <- datfromstr(a[pos1 + 
  #                                                  i])
  # }
  pos1 <- grep("# Fishing mortality by age class \\(across\\) and year \\(down\\)", 
               a)
  if (viewerVer.num < 3) {
    FbyAgeYr <- datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
  }
  else {
    FbyAgeYr <- if (nSp == 1) {
      datfromstr(a[(pos1 + 2):(pos1 + 1 + nTimes)])
    }
    else {
      nAges1 <- nAges[1]
      xxx <- scan(text = a[pos1 + 2:(nTimes * nSp + 100)], 
                  n = nTimes * nAges1 * nSp, what = 0, comment.char = "#", 
                  quiet = !verbose)
      dim(xxx) <- c(nAges1, nTimes, nSp)
      aperm(xxx, c(2, 1, 3))
    }
  }
  if (verbose) 
    cat("L134 ; ")
  yrs <- Year1 + (0:(nTimes - 1))/nRecs.yr + ifelse(nRecs.yr == 
                                                      1, 0, 1/(2 * nRecs.yr))
  alltimes <- sort(unique(as.vector(Rlz.t.fsh)))
  pos1 <- grep("# Fishing mortality by age class \\(across\\), year \\(down\\) and region \\(block\\)", 
               a)
  nAges1 <- if (nSp > 1) {
    nAges[1]
  }
  else {
    nAges
  }
  FatYrAgeReg <- array(dim = c(nTimes, nAges1, nReg))
  xxx <- scan(text = a[pos1 + 1:(nTimes * nReg + 100)], n = nTimes * 
                nAges1 * nReg, what = 0, comment.char = "#", quiet = !verbose)
  dim(xxx) <- c(nAges1, nTimes, nReg)
  FatYrAgeReg <- aperm(xxx, c(2, 1, 3))
  if (verbose) 
    cat("L154 ; ")
  pos1 <- grep("# Population Number by age \\(across\\), year \\(down\\) and region", 
               a, ignore.case = TRUE)
  nAges1 <- if (nSp > 1) {
    nAges[1]
  }
  else {
    nAges
  }
  NatYrAgeReg <- array(dim = c(nTimes, nAges1, nReg))
  if (verbose) 
    cat("L158 ; ")
  xxx <- scan(text = a[pos1 + 1:(nTimes * nReg + 100)], n = nTimes * 
                nAges1 * nReg, what = 0, comment.char = "#", quiet = !verbose)
  dim(xxx) <- c(nAges1, nTimes, nReg)
  NatYrAgeReg <- aperm(xxx, c(2, 1, 3))
  if (verbose) 
    cat("L163 ; ")
  pos1 <- grep("# Exploitable population biomass by fishery (down) and by year-season  (across)", 
               a, fixed = T)
  NexpbyYrFsh <- if (length(pos1) != 0) {
    xxx <- scan(text = a[pos1 + 1:nFisheries], n = nTimes * 
                  nFisheries, what = 0, comment.char = "#", quiet = !verbose)
    dim(xxx) <- c(nTimes, nFisheries)
    t(xxx)
  }
  else {
    NA
  }
  pos1 <- grep(a, pattern = "# Exploitable population in same units as catch by fishery (down) and year-season (across)", 
               fixed = T)
  ExPopCUnitsbyYrFsh <- if (length(pos1) != 0) {
    xxx <- scan(text = a[pos1 + 1:nFisheries], n = nTimes * 
                  nFisheries, what = 0, comment.char = "#", quiet = !verbose)
    dim(xxx) <- c(nTimes, nFisheries)
    t(xxx)
  }
  else {
    NA
  }
  pos1 <- grep(x = a, pattern = "# Recruitment")
  Recruitment <- datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
  if (verbose) 
    cat("L181 ; ")
  pos1 <- grep("# Total biomass$", a)
  TotBiomass <- datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
  pos1 <- grep("# Adult biomass$", a)
  AdultBiomass <- datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
  pos1 <- grep("# Relative biomass by region \\(across\\) and year \\(down\\)$", 
               a)
  RelBiomass <- datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
  pos1 <- grep("# Observed catch by fishery \\(down\\) and time \\(across\\)", 
               a)
  ObsCatch <- datfromstr(a[pos1 + 1:nFisheries])
  pos1 <- grep("# Predicted catch by fishery \\(down\\) and time \\(across\\)", 
               a)
  PredCatch <- datfromstr(a[pos1 + 1:nFisheries])
  pos1 <- grep("# Observed CPUE by fishery \\(down\\) and time \\(across\\)", 
               a)
  ObsCPUE <- datfromstr(a[pos1 + 1:nFisheries])
  pos1 <- grep("# Predicted CPUE by fishery \\(down\\) and time \\(across\\)", 
               a)
  PredCPUE <- datfromstr(a[pos1 + 1:nFisheries])
  BH.ar <- NULL
  pos1 <- grep("# Yield analysis option", a)
  YieldOption <- as.numeric(a[pos1 + 1])
  if (YieldOption == 1) {
    if (length(grep("# BH autocorrelation", a)) > 0) {
      BH.ar <- as.numeric(strsplit(a[grep("# BH autocorrelation", 
                                          a)], split = " +")[[1]][6])
      pos.offset <- 4
    }
    else {
      BH.ar <- NULL
      pos.offset <- 3
    }
    tmp <- strsplit(a[grep("Beverton-Holt stock-recruitment relationship report", 
                           a) + 1], split = " +")[[1]]
    columns.vec <- if (length(tmp) == 13) {
      c(4, 7, 10, 13)
    }
    else {
      c(4, 7, 10)
    }
    names(columns.vec) <- if (length(tmp) == 13) {
      c("alpha", "beta", "steepness", "varRdev")
    }
    else {
      c("alpha", "beta", "steepness")
    }
    tmp <- as.numeric(tmp[columns.vec])
    names(tmp) <- names(columns.vec)
    SRR <- tmp
    SPR0 <- 4 * SRR[2] * SRR[3]/(SRR[1] * (1 - SRR[3]))
    names(SPR0) <- "SPR0"
    SSB0 <- SRR[1] * SPR0 - SRR[2]
    names(SSB0) <- "SSB0"
    R0 <- SRR[1] - SRR[2]/SPR0
    names(R0) <- R0
    names(R0) <- "R0"
  }
  else {
    SRR <- NA
    SPR0 <- NULL
    SSB0 <- NULL
    R0 <- NULL
  }
  if (verbose) 
    cat("L225 ; ")
  pos1 <- grep("# Observed spawning Biomass", a)
  Obs.SB <- if (length(pos1) > 0) {
    datfromstr(a[pos1 + 1])
  }
  else {
    NULL
  }
  pos1 <- grep("# Observed recruitment", a)
  Obs.R <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("# Spawning Biomass", a)
  Pred.SB <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("# Predicted recruitment", a)
  Pred.R <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("steepness =", a)
  steep <- as.numeric(unlist(strsplit(a[pos1], split = "[[:blank:]]+"))[10])
  pos1 <- grep("# MSY", a)
  MSY <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# F multiplier at MSY", a)
  Fmult <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# F at MSY", a)
  Fmsy <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# Total biomass at MSY", a)
  Bmsy <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# Adult biomass at MSY", a)
  SBmsy <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# current Total Biomass to Total biomass at MSY", 
               a)
  B.Bmsy <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# current Adult Biomass to Adult Biomass at MSY", 
               a)
  SB.SBmsy <- as.numeric(a[pos1 + 1])
  pos1 <- grep("# Effort multiplier", a)[1]
  Effmult <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("# Equilibrium yield", a)
  Eq.yield <- charVec2numVec(a[pos1 + 1])
  YFcurr <- Eq.yield[Effmult == 1]
  if (length(YFcurr) == 0) 
    YFcurr <- 0
  pos1 <- grep("# Equilibrium adult biomass", a)
  Eq.SB <- charVec2numVec(a[pos1 + 1])
  SBFcurr <- Eq.SB[Effmult == 1]
  SB0 <- Eq.SB[Effmult == 0]
  if (length(SBFcurr) == 0) 
    SBFcurr <- 0
  Eq.B <- getNumVector(keyword = "# Equilibrium total biomass", 
                       a)
  BFcurr <- Eq.B[Effmult == 1]
  if (length(BFcurr) == 0) 
    BFcurr <- 0
  B0 <- Eq.B[Effmult == 0]
  Eq.SB.SBmsy <- getNumVector(keyword = "# Adult biomass over adult biomass at MSY", 
                              a)
  Eq.B.Bmsy <- getNumVector(keyword = "# Total biomass over total biomass at MSY", 
                            a)
  pos1 <- grep("# Aggregate F over F at MSY", a)
  Eq.F.Fmsy <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("# Aggregate F$", a)
  Eq.aggF <- charVec2numVec(a[pos1 + 1])
  if (nSp == 1) {
    pos1 <- grep("# Effort multiplier", a)[2]
    YPR.effmult <- charVec2numVec(a[pos1 + 1])
    pos1 <- grep("# Yield per recruit", a)[2]
    YPR <- charVec2numVec(a[pos1 + 1])
  }
  else {
    YPR.effmult <- NULL
    YPR <- NULL
  }
  TagGrps <- RepRates <- nTagRetPds <- TagRetPds <- ObsTagReturns <- PredTagReturns <- MaxLiberty <- ObsvPredbyLib <- 0
  pos1 <- grep("# Grouping indicator \\(0", a)
  if (length(pos1) > 0) {
    TagGrps <- charVec2numVec(a[pos1 + 1])
    if (max(TagGrps) == 0) {
      TagGrps <- 1:length(TagGrps)
    }
    pos1 <- grep("# Reporting rates by fishery \\(no time", 
                 a)
    xxx <- grep("# Reporting rates by by fishery by tag group \\(no time", 
                a)
    if (length(pos1) > 0) {
      RepRates <- as.numeric(a[(pos1 + 1):(pos1 + nFisheries)])
    }
    if (length(xxx) > 0) {
      pos1 <- xxx
      rm(xxx)
      RepRates <- sapply(a[(pos1 + 1):(pos1 + nFisheries)], 
                         datfromstr, USE.NAMES = F)
    }
    pos1 <- grep("# No. of time periods associated with tag returns", 
                 a)
    nTagRetPds <- charVec2numVec(a[pos1 + 1])
    pos1 <- grep("# Time periods associated with ", a)
    TagRetPds <- charVec2numVec(a[pos1 + 1])
    pos1 <- grep("# Observed tag returns by time period", 
                 a)
    ObsTagReturns <- sapply(a[(pos1 + 1):(pos1 + max(TagGrps))], 
                            datfromstr, USE.NAMES = F)
    pos1 <- grep("# Predicted tag returns by time period", 
                 a)
    PredTagReturns <- sapply(a[(pos1 + 1):(pos1 + max(TagGrps))], 
                             datfromstr, USE.NAMES = F)
    pos1 <- grep("# Maximum time at liberty", a)
    MaxLiberty <- charVec2numVec(a[pos1 + 1])
    pos1 <- grep("# Observed vs predicted tag returns by time at liberty", 
                 a)
    ObsvPredbyLib <- t(sapply(a[(pos1 + 1):(pos1 + MaxLiberty)], 
                              datfromstr, USE.NAMES = F))
  }
  if (nReg > 1) {
    pos1 <- grep("# Movement analysis", a)
    if (length(pos1) > 0) {
      nMPars <- grep("# Region 2", a[pos1:length(a)])[1] - 
        3
      MoveRates <- array(dim = c(nMPars, nReg, nReg))
      for (i in 1:nReg) {
        MoveRates[, , i] <- t(sapply(a[(pos1 + 2):(pos1 + 
                                                     nMPars + 1)], datfromstr, USE.NAMES = F))
        pos1 <- pos1 + nMPars + 1
      }
    }
    else {
      MoveRates <- NA
    }
  }
  else {
    MoveRates <- NA
  }
  if (verbose) 
    cat("L299 ; ")
  pos1 <- grep("# length-sample components of likelihood by fishery", 
               a)
  lenLiks <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("# weight-sample components of likelihood by fishery", 
               a)
  wtLiks <- charVec2numVec(a[pos1 + 1])
  pos1 <- grep("# Total biomass in absence of fishing", a)
  if (length(pos1) > 0) {
    TotalBiomass.nofish <- datfromstr(a[(pos1 + 1):(pos1 + 
                                                      nTimes)])
    pos1 <- grep("# Adult biomass in absence of fishing", 
                 a)[1]
    AdultBiomass.nofish <- if (length(pos1) > 0) {
      datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
    }
    else {
      NA
    }
    pos1 <- grep("# Exploitable populations in absence of fishing", 
                 a)[1]
    ExplBiomass.nofish <- if (length(pos1) > 0) {
      datfromstr(a[(pos1 + 1):(pos1 + nTimes)])
    }
    else {
      NA
    }
    pos1 <- grep("# Predicted catch for interaction analysis", 
                 a)[1]
    if (length(pos1) > 0) {
      PredCatch.interact <- matrix(nrow = nFisheries, 
                                   ncol = max(nRlz.fsh))
      for (i in 1:nFisheries) {
        PredCatch.interact[i, 1:nRlz.fsh[i]] <- datfromstr(a[pos1 + 
                                                               i])
      }
    }
    else PredCatch.interact <- NA
  }
  else {
    TotalBiomass.nofish <- AdultBiomass.nofish <- ExplBiomass.nofish <- PredCatch.interact <- NULL
  }
  tSelAtAge <- t(SelAtAge)
  if (verbose) 
    cat("L318 ; ")
  if (DEBUG) 
    browser()
  tblocks <- !all(SelexTblocks == 1)
  nfish <- nFisheries
  nfishWTblocks <- ifelse(tblocks, sum(SelexTblocks)/nSp, 
                          nFisheries/nSp)
  if (verbose) 
    cat("L322 ; ")
  if (DEBUG) 
    browser()
  FL0 <- unlist(sapply(1:nFisheries, function(i) {
    nblk <- SelexTblocks[i]
    tmp <- if (nblk == 1) {
      paste(i)
    }
    else {
      paste(i, 1:nblk, sep = "_")
    }
    if (i < 10) {
      paste0("0", tmp)
    }
    else {
      tmp
    }
  }))
  if (verbose) 
    cat("L328 ; ")
  if (DEBUG) 
    browser()
  dimnames(tSelAtAge)[[1]] <- paste0(1:dim(tSelAtAge)[1])
  dimnames(tSelAtAge)[[2]] <- if (nSp > 1) {
    paste0("FL", FL0, c(rep("Male", nfishWTblocks), rep("Female", 
                                                        nfishWTblocks)))
  }
  else {
    paste0("FL", FL0)
  }
  tSelAtAge <- apply(tSelAtAge, 1:2, as.numeric)
  SelAtAge.wide <- data.table::as.data.table(t(tSelAtAge))
  SelAtAge.wide$Gender <- if (nSp > 1) {
    c(rep("Male", nfishWTblocks), rep("Female", nfishWTblocks))
  }
  else {
    rep("Both", nfishWTblocks)
  }
  fishlab <- paste(rep(FL0, 2), paste0("FL", rep(substr(FL0, 
                                                        1, 2), 2)), sep = "_")
  if (verbose) 
    cat("L341 ; ")
  SelAtAge.wide$Fishery <- if (nSp == 1) {
    fishlab[1:nfishWTblocks]
  }
  else {
    rep(fishlab[1:nfishWTblocks], 2)
  }
  if (verbose) 
    cat("L348 ;")
  mean.LatAge.long <- mean.LatAge %>% as.data.frame() %>% 
    mutate(Gender = rownames(.)) %>% gather(key = AgeClass, 
                                            value = mean.LatAge, -!!sym("Gender"))
  mean.WatAge.long <- mean.WatAge %>% as.data.frame() %>% 
    mutate(Gender = rownames(.)) %>% gather(key = AgeClass, 
                                            value = mean.WatAge, -!!sym("Gender"))
  sd.LatAge.long <- sd.LatAge %>% as.data.frame() %>% mutate(Gender = rownames(.)) %>% 
    gather(key = AgeClass, value = sd.LatAge, -!!sym("Gender"))
  SelAtAge.long <- SelAtAge.wide %>% unite(col = "Fishery_Gender", 
                                           !!!syms(c("Fishery", "Gender")), sep = "-") %>% tidyr::gather(key = "AgeClass", 
                                                                                                         value = "selex", -!!sym("Fishery_Gender")) %>% separate(col = "Fishery_Gender", 
                                                                                                                                                                 into = c("Fishery", "Gender"), sep = "-") %>% mutate(Age = as.numeric(!!sym("AgeClass")))
  if (verbose) 
    cat("L356 ; ")
  if (DEBUG) 
    browser()
  SelAtAge.long %<>% inner_join(., mean.LatAge.long, by = c("Gender", 
                                                            "AgeClass")) %>% inner_join(., mean.WatAge.long, by = c("Gender", 
                                                                                                                    "AgeClass")) %>% inner_join(., sd.LatAge.long, by = c("Gender", 
                                                                                                                                                                          "AgeClass"))
  if (verbose) 
    cat("L360 ; ")
  if (DEBUG) 
    browser()
  rep.obj <- list(nTimes = nTimes, Year1 = Year1, nReg = nReg, 
                  nAges = nAges, nRecs.yr = nRecs.yr, yrs = yrs, alltimes = alltimes, 
                  nFisheries = nFisheries, nRlz.fsh = nRlz.fsh, Region.fsh = Region.fsh, 
                  Rlz.t.fsh = Rlz.t.fsh, mean.LatAge = mean.LatAge, sd.LatAge = sd.LatAge, 
                  mean.WatAge = mean.WatAge, MatAge = MatAge, SelAtAge = SelAtAge, 
                  qAtAge = qAtAge, FbyAgeYr = FbyAgeYr, 
                  FatYrAgeReg = FatYrAgeReg, NatYrAgeReg = NatYrAgeReg, 
                  NexpbyYrFsh = NexpbyYrFsh, ExPopCUnitsbyYrFsh = ExPopCUnitsbyYrFsh, 
                  Recruitment = Recruitment, TotBiomass = TotBiomass, 
                  AdultBiomass = AdultBiomass, RelBiomass = RelBiomass, 
                  ObsCatch = ObsCatch, PredCatch = PredCatch, ObsCPUE = ObsCPUE, 
                  PredCPUE = PredCPUE, YieldOption = YieldOption, SRR = SRR, 
                  Obs.SB = Obs.SB, Obs.R = Obs.R, Pred.SB = Pred.SB, Pred.R = Pred.R, 
                  steep = steep, MSY = MSY, Fmult = Fmult, Fmsy = Fmsy, 
                  Bmsy = Bmsy, SBmsy = SBmsy, B.Bmsy = B.Bmsy, SB.SBmsy = SB.SBmsy, 
                  Effmult = Effmult, Eq.yield = Eq.yield, YFcurr = YFcurr, 
                  Eq.SB = Eq.SB, SBFcurr = SBFcurr, SB0 = SB0, Eq.B = Eq.B, 
                  B0 = B0, BFcurr = BFcurr, Eq.SB.SBmsy = Eq.SB.SBmsy, 
                  Eq.B.Bmsy = Eq.B.Bmsy, Eq.F.Fmsy = Eq.F.Fmsy, Eq.aggF = Eq.aggF, 
                  YPR.effmult = YPR.effmult, YPR = YPR, TagGrps = TagGrps, 
                  RepRates = RepRates, nTagRetPds = nTagRetPds, TagRetPds = TagRetPds, 
                  ObsTagReturns = ObsTagReturns, PredTagReturns = PredTagReturns, 
                  MaxLiberty = MaxLiberty, ObsvPredbyLib = ObsvPredbyLib, 
                  MoveRates = MoveRates, lenLiks = lenLiks, wtLiks = wtLiks, 
                  TotalBiomass.nofish = TotalBiomass.nofish, AdultBiomass.nofish = AdultBiomass.nofish, 
                  ExplBiomass.nofish = ExplBiomass.nofish, PredCatch.interact = PredCatch.interact, 
                  version = viewerVer.num, nSp = nSp, spPtr = spPtr, spSexPtr = spSexPtr, 
                  regSpPtr = regSpPtr, filename = rep.file, frqfilename = frqfilename, 
                  SelexSeasons = SelexSeasons, SelexTblocks = SelexTblocks, 
                  SSB0 = SSB0, SPR0 = SPR0, R0 = R0, BH.ar = BH.ar, SelAtAge.wide = SelAtAge.wide, 
                  SelAtAge.long = SelAtAge.long, inputParfilename = inputParfilename, 
                  outputParfilename = outputParfilename, MFCL.version.string = MFCL.version.string, 
                  MFCL.version = MFCL.version)
  return(invisible(rep.obj))
}
