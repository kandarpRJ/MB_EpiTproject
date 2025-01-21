library(karyoploteR)
m6a <- import.bed("/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_m6A.cov4_pmod10_mod2_omod0.bed")
m5c <- import.bed("/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_m5C.cov4_pmod10_mod2_omod0.bed")
ino <- import.bed("/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_I.cov4_pmod10_mod2_omod0.bed")
psU <- import.bed("/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_psU.cov4_pmod10_mod2_omod0.bed")
kp <- plotKaryotype(genome = "hs1")
kpPlotRegions(kp, data=m6a, r0=0, r1=0.5, data.panel = 1, col = "darkred")
kpPlotRegions(kp, data=m5c, r0=0.5, r1=1, data.panel = 1, col = "darkblue")
kpPlotRegions(kp, data=I, r0=0, r1=0.5, data.panel = 2, col = "darkgreen")
kpPlotRegions(kp, data=psU, r0=0.5, r1=1, data.panel = 2, col = "darkorange")


png("m6a_chr_density.png", res = 600, width = 4000, height = 3800)
kp<-plotKaryotype(genome = "hs1")
kpPlotDensity(kp, data=m6a, col="darkred")
dev.off()

png("m5c_chr_density.png", res = 600, width = 4000, height = 3800)
kp<-plotKaryotype(genome = "hs1")
kpPlotDensity(kp, data=m5c, col="darkblue")
dev.off()

png("inosine_chr_density.png", res = 600, width = 4000, height = 3800)
kp<-plotKaryotype(genome = "hs1")
kpPlotDensity(kp, data=ino, col="darkgreen")
dev.off()

png("psU_chr_density.png", res = 600, width = 4000, height = 3800)
kp<-plotKaryotype(genome = "hs1")
kpPlotDensity(kp, data=psU, col="darkorange")
dev.off()






