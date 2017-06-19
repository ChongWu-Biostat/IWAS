suppressMessages(library(glmnet))
suppressMessages(library('plink2R'))

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


## parameter set-up

args = commandArgs(TRUE)
job = (eval(parse(text=args[[1]])))

job = as.numeric(job)

pheno = "12thick"
pfx = "V7"
phase = "AD"
#pheno = (as.character(args[[2]]))
#pfx = (as.character(args[[3]]))
#phase = (as.character(args[[4]]))


outd=paste("/home/panwei/wuxx0845/TWAS/ADNI/",pfx,"_",phase,"_",pheno,sep="")
datadir = paste0("/home/panwei/wuxx0845/TWAS/ADNI/data/temp/",pfx,"_",phase,"_",pheno)

system(paste("mkdir -p ", outd, sep=""))
system(paste("mkdir -p ", datadir, sep=""))

setwd("/home/panwei/wuxx0845/TWAS/ADNI/data/")

dat <- readRDS("ADNI_cross_sectional.rds") #cross-sectional
col.name = colnames(dat)

## thickness, 12 ROIs related to the default mode network
Voldmn <- c(
"ST31TA", "ST90TA",  # Inferior parietal
"ST32TA", "ST91TA",  # Inferior temporal
"ST39TA", "ST98TA",  # Medial orbitofrontal
"ST44TA", "ST103TA", # Parahippocampal
"ST52TA", "ST111TA", # Precuneus
"ST50TA", "ST109TA") # Posterior cingulate


# Other ROIs
#Voldmn2 = col.name[grepl("TA$", col.name)]
#Voldmn = c(Voldmn, Voldmn2[!Voldmn2 %in% Voldmn]) # other ROIs

cov <-dat[dat$t == 0,] # we use baseline information

cov = cov[!duplicated(cov[,"PTID"]),]
rownames(cov) = cov[,"PTID"]

if(phase != "ALL") {
    cov = cov[cov[,"DX_bl"]== "CN" |  cov[,"DX_bl"]== "AD",]
}

hg19 = read.table("glist-hg19.txt")
hg19 = hg19[order(hg19[,1]),]
hg19 = hg19[hg19[,1]!= "X" & hg19[,1]!= "XY" & hg19[,1]!= "Y", ]

hg19[,1] = as.numeric(hg19[,1])
hg19[,2] = as.numeric(hg19[,2])
hg19[,3] = as.numeric(hg19[,3])

colnames(hg19) = c("chr","startbp","endbp","gene")

expand = 1e6
alpha = 0.5 # elastic-net with alpha
hg19$startbp = hg19$startbp - expand; hg19$startbp[hg19$startbp<0] <- 0
hg19$endbp = hg19$endbp + expand

each = ceiling(nrow(hg19)/500)
last.job.num = floor(nrow(hg19)/each)
if(job==last.job.num) geneset = hg19[(job*each+1):nrow(hg19),] else geneset = hg19[(job*each+1):(job*each+each),]


cov$pteducat = as.numeric(as.character(cov$PTEDUCAT))
cov$icv_bl = as.numeric(as.character(cov$ICV_bl))/10^5
cov$age = as.numeric(as.character(cov$AGE))
cov$ptgender = ifelse(cov$PTGENDER=="Female",1,0)
cov$pthand = as.numeric(as.character(cov$PTHAND))-1

Z = as.matrix(cov[,c("pthand","ptgender","pteducat","icv_bl","age")])
Ymat = matrix(nrow=nrow(Z),ncol=length(Voldmn))


z.rowname = rownames(Z)
rownames(Ymat) = z.rowname
colnames(Ymat) = Voldmn

for(i in 1:length(Voldmn)){
    y = as.numeric(as.character(cov[,Voldmn[i]]))
    fit1 = glm(y~Z,family="gaussian")
    yres = y-fitted.values(fit1)
    Ymat[,i] = yres
    yres = NA
}

# count = 0

out.inf = as.data.frame(matrix(NA,nrow(geneset),5))

for(i in 1:nrow(geneset)){
    tryCatch({
        #     i=124
        system(paste0("plink --bfile ","/home/panwei/wuxx0845/Imputed_ADNI1_ADNI2_1000G/ADNI1/chr",geneset[i,1],
        " --chr ",geneset[i,1]," --from-kb ",geneset[i,2]/1e3," --to-kb ",
        geneset[i,3]/1e3," --maf 0.05 --hwe 0.05  --extract Hapmap_SNPlist.txt --make-bed --out ",datadir,"/", geneset[i,4]))
        
        genos = read_plink(paste0(datadir,"/",geneset[i,4]),impute="avg")
        system(paste0("rm ",datadir,"/",geneset[i,4],"*"))
        
        #### extract X0:the un-weighted SNPs and Y
        
        ### remove ambiguous SNPs
        qc = allele.qc(  genos$bim[,5]  ,genos$bim[,6] , genos$bim[,5] , genos$bim[,6] )
        
        if ( sum(!qc$keep) > 0 ) {
            genos$bim = genos$bim[qc$keep,]
            genos$bed = genos$bed[,qc$keep]
        }
        
        tmp.geno = genos$bed #353

        
        tmp = strsplit(rownames(tmp.geno),"_")
        tmp = unlist(tmp)
        tmp.len = length(tmp)
        raw.id2 = paste(tmp[1:tmp.len %%4==2],tmp[1:tmp.len %%4==3],tmp[1:tmp.len %%4==0],sep="_")
        rownames(tmp.geno) = raw.id2
        
        used.subj = intersect(raw.id2,cov[,"PTID"])
        tmp.geno = tmp.geno[used.subj,]
        Z = Z[used.subj,]
        Ymat = Ymat[used.subj,]
        
        X0 = scale(tmp.geno)
        
        refAllele = genos$bim[,5]
        alt = genos$bim[,6]
        
        X0 = as.matrix(X0)
        
        Xres = X0
        
        for(j in 1:ncol(X0)){
            fit2 = glm(X0[,j]~Z,family="gaussian")
            Xres[,j] = X0[,j]-fitted.values(fit2)
        }
        
        Xres = scale(Xres)
        eff.wgt = matrix( 0 , nrow=ncol(Xres), ncol=ncol(Ymat))
        R2 = rep(NA, ncol = ncol(Ymat))
        for(k2 in 1:ncol(Ymat)){
            yres = Ymat[,k2]
            sds = apply( Xres  , 2 , sd )
            keep = sds != 0 & !is.na(sds)
            set.seed(123)
            enet = cv.glmnet( x=Xres[,keep] , y=yres , alpha=alpha , nfold=5 , intercept=T , standardize=F )
            eff.wgt[ keep, k2 ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
            # lasso.fit = glmnet(Xres,yres,alpha=alpha,family="gaussian")
            # cv.out = cv.glmnet(Xres,yres,alpha=alpha,family="gaussian")
            # bestlam = cv.out$lambda.min
            # lasso.coef[,k2] = predict(lasso.fit,type="coefficients",s=bestlam)[-1]
            yhat <- Xres %*% eff.wgt[,k2]
            R2[k2] <- cor(cbind(yres, yhat))[1,2]^2
            print(k2)
        }
        
        apply(eff.wgt, 2, function(x) sum(x!=0))
        R2 = matrix(rep(R2, nrow(eff.wgt)), nrow = nrow(eff.wgt), byrow = T)
        out = cbind(as.character(geneset[i,4]),geneset[i,1],geneset[i,2],geneset[i,3],ncol(X0),
        colnames(X0),refAllele,alt,alpha,eff.wgt,R2)
        
        out.inf[i,] = c(as.character(geneset[i,4]),geneset[i,1],geneset[i,2],geneset[i,3],ncol(X0))
        
        saveRDS(out.inf,paste(outd,"/out_",job,".rds",sep=""))
        write.table(out.inf,paste(outd,"/out_inf_",job,".txt",sep=""),append=T,col.names=F,row.names=F,quote=F)
        
        colnames(out) = c("gene","chr","start","end","Nsnp","SNP","ref","alt","alpha",Voldmn,paste0("R2_",Voldmn))
        
        out = out[apply(eff.wgt, 1, function(x) sum(abs(x)))!=0,]
        ### if only one row left
        if(is.null(nrow(out))) out = matrix(out, nrow=1)
        
        if(sum(abs(eff.wgt))!=0) {
            write.table(out,paste(outd,"/out_",as.character(geneset[i,4]),".txt",sep=""),append=T,col.names=F,row.names=F,quote=F)
            saveRDS(out,paste(outd,"/out_",as.character(geneset[i,4]),".rds",sep=""))
        }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


saveRDS(out.inf,paste(outd,"/out_",job,".rds",sep=""))
write.table(out.inf,paste(outd,"/out_inf_",job,".txt",sep=""),append=T,col.names=F,row.names=F,quote=F)
