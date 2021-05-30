library(MEDIPS)
library(BSgenome.Mmusculus.UCSC.mm10)

dips <- function(x) {
  mm10 = "BSgenome.Mmusculus.UCSC.mm10"
  dip = MEDIPS.createSet(
    file = x,
    BSgenome = mm10,
    extend = 300,
    shift = 0,
    window_size = 50,
    paired = T
  )
  CS = MEDIPS.couplingVector(pattern = "CG",
                             refObj = dip)
  MEDIPS.exportWIG(
    Set = dip,
    CSet = CS,
    file = paste(
      sub(pattern = ".bam", replacement = "", x),
      ".wig",
      sep = ""
    ),
    format = "rms"
  )
}

dips("1KO_5mC_G1.bam")
dips("1KO_5mC_G2.bam")
dips("1KO_5hmC_G1.bam")
dips("1KO_5hmC_G2.bam")