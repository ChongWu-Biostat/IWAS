
setwd("/Users/chong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/IWAS_source")
library(data.table)
library(matlib)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(MASS)

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
suppressMessages(library(aSPU2))

source("dist_support.R")
# Some codes are directly copied from TWAS source codes
option_list = list(
make_option("--sumstats", action="store", default=NA, type='character',
help="summary statistics (rds file and must have SNP and Z column headers) [required]"),
make_option("--out", action="store", default=NA, type='character',
help="Path to output file [required]"),
make_option("--weights", action="store", default=NA, type='character',
help="File listing molecular weight (rds files and must have columns ID,CHR,P0,P1, Weights) [required]"),
make_option("--ref_ld", action="store", default=NA, type='character',
help="Reference LD files in binary PLINK format [required]"),
make_option("--gene_list", action="store", default=NA, type='character',
help="Gene sets we want to analyze, currently only gene sets from a single chromosome are supported [required]"),
make_option("--test_type", action="store", default="aSPU", type='character',
help="Test we want to perform, the default is aSPU for single weights. If we want to combine mulitple weights, we can set it to daSPU [default: aSPU]"),
make_option("--weight_type", action="store", default="ST31TA", type='character',
help="Weight we want to use. [default: ST31TA]")

)

opt = parse_args(OptionParser(option_list=option_list))




opt = list(
sumstats="./Example/IGAP_chr22.rds",
out= "./Example/example_res.rds",
weights = "./WEIGHTS/ADNI1_wgt.rds",
ref_ld = "./LDREF/hapmap_CEU_r23a_hg19",
gene_list ="./Example/gene_list.txt",
test_type ="aSPU",
weight_type = "ST31CV")

# This function is from TWAS (http://gusevlab.org/projects/fusion/#typical-analysis-and-output)
allele.qc = function(a1,a2,ref1,ref2) {
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;
    
    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
    return(snp)
}

Voldmn <- c(
"ST31CV", "ST90CV",  # Inferior parietal
"ST32CV", "ST91CV", # Inferior temporal
"ST39CV", "ST98CV",  # Medial orbitofrontal
"ST44CV", "ST103CV", # Parahippocampal
"ST52CV", "ST111CV", # Precuneus
"ST50CV", "ST109CV", # Posterior cingulate
"ST29SV", "ST88SV") # Hippocampus



gene.list = read.table(opt$gene_list)
gene.list[,1] = as.character(gene.list[,1])
gene.list = gene.list[,1]

weights = readRDS(opt$weights)

weights = weights[weights[,1] %in% gene.list,]


hg19 = read.table("glist-hg19.txt")

hg19[,1] = as.numeric(levels(hg19[,1]))[hg19[,1]]
hg19[,2] = as.numeric(hg19[,2])
hg19[,3] = as.numeric(hg19[,3])


colnames(hg19) = c("CHR","startbp","endbp","gene")

expand = 1e6

hg19$startbp = hg19$startbp - expand; hg19$startbp[hg19$startbp<0] <- 0
hg19$endbp = hg19$endbp + expand

wgtlist0 = hg19[hg19[,4] %in% gene.list,] # gene list
wgtlist0[,4] = as.character(wgtlist0[,4])


snp.inf = unique(weights[,"SNP"])


sumstat.orgin = readRDS(opt$sumstats)
m = match(snp.inf, sumstat.orgin$SNP_map)
sumstat.orgin = sumstat.orgin[m,]

if(opt$test_type == "aSPU") {
    out.res = as.data.frame(matrix(NA,nrow(wgtlist0),15))
    colnames(out.res) = c("gene","CHR","P0","P1","#nonzero_SNPs","TWAS_asy","SSU_asy","SPU(1)","SPU(2)","SPU(3)","SPU(4)","SPU(5)","SPU(6)","SPU(Inf)","aSPU")
} else {
    out.res = as.data.frame(matrix(NA,nrow(wgtlist0),23))
    colnames(out.res) = c( "gene","CHR","P0","P1","#SNPs","#Weights","SPU(1,1)","SPU(2,1)","SPU(3,1)","SPU(Inf,1)","SPU(1,2)","SPU(2,2)","SPU(3,2)","SPU(Inf,2)","SPU(1,3)","SPU(2,3)","SPU(3,3)","SPU(Inf,3)","SPU(1,Inf)","SPU(2,Inf)","SPU(3,Inf)","SPU(Inf,Inf)","daSPU"  )
}

for ( w in 1:nrow(wgtlist0) ) {
    tryCatch({
        # Load in summary stats
        genename = wgtlist0$gene[w]
        genename = as.character(genename)
        sumstat = sumstat.orgin
        
        wgt.dat = weights[weights[,1] == genename,]
        
        mysnps = wgt.dat[,"SNP"]
        mysnps = unique(mysnps)
        
        write.table(mysnps,paste0("mysnps_",genename,".txt"),row.names=F,col.names=F,quote=F,append=F)
        
        system(paste0("./plink --bfile ",opt$ref_ld," --extract mysnps_", genename,".txt --maf 0.01 --make-bed --out test_",genename))
        
        # Load in reference data
        genos = read_plink(paste0("test_",genename),impute="avg")
        system(paste0("rm test_",genename,"*"))
        system(paste0("rm mysnps_",genename,"*"))
        
        # Match summary data to input, record NA where summary data is missing
        m = match( genos$bim[,2] , sumstat$SNP_map )
        sum.missing = is.na(m)
        sumstat = sumstat[m,]
        sumstat$SNP = genos$bim[,2]
        sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
        sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]
        
        # QC / allele-flip the input and output
        qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )
        
        # Flip Z-scores for mismatching alleles
        sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
        sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
        sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]
        
        # Remove strand ambiguous SNPs (if any)
        if ( sum(!qc$keep) > 0 ) {
            genos$bim = genos$bim[qc$keep,]
            genos$bed = genos$bed[,qc$keep]
            sumstat = sumstat[qc$keep,]
        }
        
        snps = mysnps
        m = match( snps , genos$bim[,2] )
        m.keep = !is.na(m)
        snps = snps[m.keep]
        wgt.matrix = wgt.dat[m[m.keep],]
        
        
        ###############
        ### revise it
        ###
        
        wgt.matrix = wgt.matrix[,colSums(is.na(wgt.matrix)) < 10]
        
        
        R2 = wgt.matrix[1,colnames(wgt.matrix) %in% paste0("CV_Rsq_",1:14)]
        R2.name = names(R2)
        R2 = as.numeric(R2)
        names(R2) = R2.name
        R2  = R2[R2 >0.01]
        
        
        tmp.name = names(R2)
        
        tmp.index = as.numeric(gsub("CV_Rsq_","",tmp.name))
        tmp.index = tmp.index[!is.na(tmp.index)]
        
        wgt.matrix = wgt.matrix[,c("gene","chr","start","end","Nsnp","SNP","ref","alt","alpha",Voldmn[tmp.index])]
        
        genos$bed = genos$bed[,m[m.keep]]
        genos$bim = genos$bim[m[m.keep],]
        
        # Flip WEIGHTS for mismatching alleles
        qc = allele.qc( wgt.matrix[,"ref"] , wgt.matrix[,"alt"] ,  genos$bim [,5] ,  genos$bim[,6] )
        
        wgt.matrix2 = wgt.matrix[,10:dim(wgt.matrix)[2]]
        wgt.matrix2 = apply(wgt.matrix2,2,as.numeric)
        wgt.matrix2[qc$flip,]  = -1 * wgt.matrix2[qc$flip,]
        wgt.matrix[,10:dim(wgt.matrix)[2]] = wgt.matrix2
        
        if ( sum(!qc$keep) > 0 ) {
            genos$bim = genos$bim[qc$keep,]
            genos$bed = genos$bed[,qc$keep]
            sumstat = sumstat[qc$keep,]
            wgt.matrix = wgt.matrix[qc$keep,]
        }
        
        cur.genos = scale(genos$bed)
        cur.bim = genos$bim
        
        cur.FAIL = FALSE
        
        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , sumstat$SNP_map)
        cur.Z = sumstat$Z[m]
        
        # Compute LD matrix
        cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
        cur.miss = is.na(cur.Z)
        
        # Impute missing Z-scores
        if ( sum(cur.miss) != 0 ) {
            if ( sum(!cur.miss) == 0 ) {
                cat( "WARNING : " , unlist(wgtlist0[w,]) , " had no overlapping GWAS Z-scores\n")
                cur.FAIL = TRUE
            } else {
                cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
            }
        }
        
        if ( !cur.FAIL ) {
            
            
            U = cur.Z
            U = as.matrix(U)
            
            V = cur.LD
            weight = wgt.matrix[,10:dim(wgt.matrix)[2]]
            
            
            weight.name = colnames(weight)
            weight = apply(weight,2,as.numeric)
            
            # standardize the weights
            for(i in 1:dim(weight)[2]) {
                weight[,i] = weight[,i] /sum(abs(weight[,i]))
            }
            
            colnames(weight) = weight.name
            
            name = rownames(V)
            rownames(U) = name
            weight = as.matrix(weight)
            weight = weight[,colSums(abs(weight))>0]
            weight.name = colnames(weight)
            
            
            if(opt$test_type =="aSPU") {
                
                if(sum(weight.name %in% opt$weight_type) <1 ) {
                     cat( "Note : " , unlist(wgtlist0[w,]) , " had no pre-specificed weight type.")
                    next
                }
                weight = weight[,opt$weight_type]
                
                
                non_zero_ind <- (weight != 0)
                
                diag_element <- as.vector(weight[non_zero_ind])
                diag_sd <- diag_element/sum(abs(diag_element))
                weight_diag <- diag(diag_sd,nrow = length(diag_sd))
                
                Zstat.w <- weight_diag %*% cur.Z[non_zero_ind]
                corSNP.w <- weight_diag %*% cur.LD[non_zero_ind, non_zero_ind] %*% t(weight_diag)
                
                pSum = Sum(U=Zstat.w, CovS=corSNP.w)
                pSSU = SumSqU(U=Zstat.w, CovS=corSNP.w)
                
                res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e3)
                
                if(min(res$pvs) < 5e-3) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e4)
                }
                
                if(min(res$pvs) < 5e-4) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e5)
                }
                
                if(min(res$pvs) < 5e-5) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e6)
                }
                
                if(min(res$pvs) < 5e-6) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e7)
                }
                
                out.res[w,] = c(as.character(wgtlist0[w,"gene"]),wgtlist0[w,"CHR"],wgtlist0[w,"startbp"]- expand,wgtlist0[w,"endbp"] - expand,sum(weight!=0),pSum,pSSU,res$pvs)
            
                
            } else {
                daSPU.p = daSPU(U,V,weight,pow1=c(1,2,3,Inf),pow2=c(1,2,3,Inf),n.perm = 1e3)
                
                if(min(daSPU.p) < 5e-3) {
                    daSPU.p = daSPU(U,V,weight,pow1=c(1,2,3,Inf),pow2=c(1,2,3,Inf),n.perm = 1e4)
                }
                
                
                if(min(daSPU.p) < 5e-4) {
                    daSPU.p = daSPU(U,V,weight,pow1=c(1,2,3,Inf),pow2=c(1,2,3,Inf),n.perm = 1e5)
                }
                
                if(min(daSPU.p) < 5e-5) {
                    daSPU.p = daSPU(U,V,weight,pow1=c(1,2,3,Inf),pow2=c(1,2,3,Inf),n.perm = 1e6)
                }
                
                if(min(daSPU.p) < 5e-6) {
                    daSPU.p = daSPU(U,V,weight,pow1=c(1,2,3,Inf),pow2=c(1,2,3,Inf),n.perm = 1e7)
                }
                
                pvalue = daSPU.p
                
                tmp.res = c(as.character(wgtlist0[w,"gene"]),wgtlist0[w,"CHR"],wgtlist0[w,"startbp"]-expand,wgtlist0[w,"endbp"]-expand,dim(weight)[1],dim(weight)[2],pvalue)
                out.res[w,] = tmp.res
            }
        }
        
        print(w)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

out.res = out.res[!is.na(out.res[,1]),]
saveRDS(out.res,opt$out)

write.table(out.res, "output.txt")
