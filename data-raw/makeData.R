###
### mIHC Package Data
###

setwd("/Users/hortowe/my_tool_repos/spatialMIHC")
library(colorspace)
library(data.table)

###
### AP DP20 Colors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

ap_dp20Colors_dt <- data.table("Population" = c("Mono/macs", "Mono/macs CSF1R+", 
                                                "Mono/macs CSF1R-", "Dendritic cells", "T cells", "CD8 T cells", 
                                                "CD3+ CD8- T cells", "B cells", "PanCK+", "PanCK+ CC1+ (pyroptosis)", 
                                                "PanCK+ Gasdermin+ (pyroptosis)", "PanCK+ CC8+ CC3- (early apoptosis)", 
                                                "PanCK+ CC3+ (apoptosis)", "PanCK+ H2AX+ (apoptosis)", "PanCK+ CD71+ (ferroptosis)", 
                                                "PanCK+ S-phase proliferation", "PanCK+ broad proliferation", 
                                                "Total CC1", "Total Gasdermin", "Total CC3", "Total CC8", "Total H2AX", 
                                                "Total CD71", "CSF1R+", "CC1+", "GASDERMIN+", "CC3+", "CC8+", 
                                                "H2AX+", "CD71+", "BRDU+", "KI67+"),
                               "Hex" = c("#CA3BD4", "#6D0C6A", "#CB88E7", "#F8EB74", "#abe1ed", "#DE2D31", "#B4C6E7", "#F99539", 
                                         "#833C0C", "#FF64AA", "#87dec2", "#03B050", "#81F1F7", "#7F71FF", "#FDE834", "#F3CFC4",
                                         "#f097ac", "#87dec2", "#FF64AA", "#7031A0", "#03B050", "#135E8E", "#F99539", "#7F71FF",
                                         "#87dec2", "#FF64AA", "#7031A0", "#03B050", "#135E8E", "#F99539", "#CB88E7", "#FDE834"))

ap_dp20Colors_v <- ap_dp20Colors_dt$Hex
names(ap_dp20Colors_v) <- ap_dp20Colors_dt$Population

###
### RCN Colors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

rcnColors_dt <- data.table("RCN" = 1:21, 
                           "Hex" = c(RColorBrewer::brewer.pal(7, "Set2"), 
                                     RColorBrewer::brewer.pal(9, "Set1"), 
                                     RColorBrewer::brewer.pal(5, "Dark2")))

rcnColors_v <- rcnColors_dt$Hex
names(rcnColors_v) <- rcnColors_dt$RCN

###
### Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

usethis::use_data(ap_dp20Colors_v, overwrite = T)
usethis::use_data(rcnColors_v, overwrite = T)