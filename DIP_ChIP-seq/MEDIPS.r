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
      sub(pattern = "trimmed.bam", replacement = "", x),
      "MEDIPS_bs-50.wig",
      sep = ""
    ),
    format = "rms"
  )
}

dips("1KO_medip-G1.bam")
dips("1KO_medip-G2.bam")
dips("1KO_hmedip-G1.bam")
dips("1KO_hmedip-G2.bam")