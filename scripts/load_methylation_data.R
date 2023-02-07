#!/usr/bin/env R

# https://github.com/santoesha/catnon/blob/master/scripts/R/catnon_methylation.R

# example CSV:
# "seqID","BatchID","Slide","Array","MethylationBatch","sex","whord","cortrd","Trt3","age","mgmtold","pfs","ss","tss","tpfs","biomicsID","IDH1mutation","IDH1status","IDH2mutation","IDH2status","IDH","MethylationStatus"
# 2,"203197470197_R08C01",203197470197,"R08C01","7","female","0","stable/decreasing dose for at least 2 weeks","RT",40.380561259,"methylated","No","No",3976,3976,"I19-1218-01","IDH1 R132H","mutant","wildtype","wildtype","mutant","methylated"
# 3,"202259490003_R08C01",202259490003,"R08C01","3","male","0","no","TMZ/RT->TMZ",33.535934292,"unmethylated","Yes","Yes",3550,2613,"I18-1076-01","IDH1 R132H","mutant","wildtype","wildtype","mutant","methylated"
# 5,"201870610161_R08C01",201870610161,"R08C01","2","female","0","stable/decreasing dose for at least 2 weeks","RT",51.605749487,"unmethylated","Yes","Yes",391,154,"I18-1076-02","wildtype","wildtype","wildtype","wildtype","wildtype","unmethylated"

a = read.csv("data/GLASS_OD/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7, stringsAsFactors = F)

