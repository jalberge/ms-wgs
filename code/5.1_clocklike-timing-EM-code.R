#Myeloma hyperdiploidy timing
#Dec 2023

options(stringsAsFactors = F)
library(data.table)

#Functions
dbinomtrunc = function(x, size, prob, minx=3) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}
#


annot_all=fread("donor_annot.tsv",data.table=F)
maf=fread("20231019_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level_combined_signatures_annotated.tsv",data.table=F,quote=F)

#annot=annot_all[annot_all$Cohort=="SU2C"&annot_all$HRDTx=="HRD",]
annot=annot_all[annot_all$HRDTx=="HRD",]
rownames(annot)=annot$"Sample ID"
hrd_chrs_poss=c(3,5,7,9,11,15,19,21)

for(donor in annot$"Sample ID"){
  maf_patient=maf[maf$Participant_ID==donor,]
  maf_patient=maf_patient[!is.na(maf_patient$IS_SCNA)&!is.na(maf_patient$SCNA_NB),]
  chr_gain=table(maf_patient$Chromosome[maf_patient$IS_SCNA&maf_patient$SCNA_NB==2&maf_patient$SCNA_NA==0&maf_patient$SCNA_tau>2.1])
  num=which(rownames(annot)==donor)
  
  if(length(chr_gain)>0){
    chr_nongain=table(maf_patient$Chromosome[!(maf_patient$IS_SCNA&maf_patient$SCNA_NB==2&maf_patient$SCNA_NA==0&maf_patient$SCNA_tau>2.1)])
    frac_gained=chr_gain/(chr_gain+chr_nongain[names(chr_gain)])
    
    gained_chrs=names(chr_gain)[!names(chr_gain)%in%names(chr_nongain)|frac_gained>0.9]
    hrd_chrs=gained_chrs[gained_chrs%in%hrd_chrs_poss]
    nonhrd_gains=gained_chrs[!gained_chrs%in%hrd_chrs_poss]
    sig_col=106:123
    maf_patient$ml_sig=apply(maf_patient[,106:123],1,function(x) colnames(maf_patient)[106:123][which.max(x)])
    # NV_vec=maf_patient$t_alt_count[mut_hrd]
    # NR_vec=maf_patient$t_alt_count[mut_hrd]+maf_patient$t_ref_count[mut_hrd]
    # res = binom_mix(NV_vec,NR_vec,nrange=1:3,mode='Full')
    # 
    # Duplicated=2
    # Tot_CN=3
    # purity=maf_patient$purity[1]
    # cutoff=0.05
    # Major_cluster=which(abs(res$p-Duplicated/Tot_CN*purity)<cutoff)
    # Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff)
    # Prop_Duplicated=res$prop[Major_cluster]
    # Prop_NonDuplicated=sum(res$prop[Minor_cluster])
    purity=maf_patient$purity[1]
    # if(donor=="IL01308") purity=0.3
    # if(donor=="IL01864") purity=0.5
     if(donor=="MBp09") purity=0.9
    # if(donor=="pM10598") purity=0.4
    # if(donor=="pM10845") purity=0.4
    # if(donor=="pM5738") purity=0.4
    
    timing_hrd_chrs=timing_hrd_chrs_high=timing_hrd_chrs_low=num_m1=num_m2=c()
    for(chr in hrd_chrs){
      mut_hrd=maf_patient$Chromosome%in%hrd_chrs&maf_patient$clonal==1&!is.na(maf_patient$clonal)&
        maf_patient$ml_sig%in%c("SBS1","SBS5")&maf_patient$Chromosome==chr
      if(sum(mut_hrd)>4){
      Duplicated=2
      Tot_CN=3

      
      simple_gain=all(sort(unique(maf_patient$multiplicity[mut_hrd]))==c(1,2))
      if(!simple_gain){
        hrd_chrs=hrd_chrs[hrd_chrs!=chr]
      }else{
        NV_vec=maf_patient$t_alt_count[mut_hrd]
        NR_vec=maf_patient$t_alt_count[mut_hrd]+maf_patient$t_ref_count[mut_hrd]
        res = binom_mix(NV_vec,NR_vec,nrange=1:3,mode='Truncated',criterion = "AIC")
        
        if(res$n==2){
          Major_cluster=which.max(res$p)
          Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff) 
        }else{
          cutoff=0.08
          Major_cluster=which(abs(res$p-Duplicated/Tot_CN*purity)<cutoff)
          Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff) 
        }
        
        if(length(Major_cluster)>0&length(Minor_cluster)>0){
        if(any(Major_cluster%in%Minor_cluster))Major_cluster=Major_cluster[!Major_cluster%in%Minor_cluster]
        Prop_Duplicated=res$prop[Major_cluster]
        Prop_NonDuplicated=sum(res$prop[Minor_cluster])   
        conf.intv=poisson.test(c(round(sum(mut_hrd)*Prop_NonDuplicated),
                                 round(sum(mut_hrd)*Prop_Duplicated)))$conf.int
        # Prop_NonDuplicated=sum(mut_hrd&maf_patient$multiplicity==1)
        # Prop_Duplicated=sum(mut_hrd&maf_patient$multiplicity==2)
         Time=Tot_CN/(Duplicated+Prop_NonDuplicated/Prop_Duplicated)
        # conf.intv=poisson.test(c(Prop_NonDuplicated,Prop_Duplicated))$conf.int
        Conf.Time1=Tot_CN/(Duplicated+conf.intv[1])
        Conf.Time2=Tot_CN/(Duplicated+conf.intv[2])
        num=which(annot$"Sample ID"==donor)
        timing_hrd_chrs=c(timing_hrd_chrs,Time)
        #num_m1=c(num_m1,sum(mut_hrd&maf_patient$multiplicity==1))
        #num_m2=c(num_m2,sum(mut_hrd&maf_patient$multiplicity==2))
        num_m1=c(num_m1,round(sum(mut_hrd)*Prop_NonDuplicated))
        num_m2=c(num_m2,round(sum(mut_hrd)*Prop_Duplicated))
        timing_hrd_chrs_low=c(timing_hrd_chrs_low,Conf.Time2)
        timing_hrd_chrs_high=c(timing_hrd_chrs_high,Conf.Time1)
        }
      }
      }
    }
    if(length(hrd_chrs)>0){
    annot$HRD_chrs[num]=paste(hrd_chrs,collapse=",")
    
    mut_hrd=maf_patient$Chromosome%in%hrd_chrs&maf_patient$clonal==1&!is.na(maf_patient$clonal)&maf_patient$ml_sig%in%c("SBS1","SBS5")
    if(sum(mut_hrd)>5){
      Duplicated=2
      Tot_CN=3
      NV_vec=maf_patient$t_alt_count[mut_hrd]
      NR_vec=maf_patient$t_alt_count[mut_hrd]+maf_patient$t_ref_count[mut_hrd]
      res = binom_mix(NV_vec,NR_vec,nrange=1:3,mode='Truncated',criterion="AIC")
      if(res$n==2){
        Major_cluster=which.max(res$p)
        Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff) 
      }else{
        cutoff=0.08
        Major_cluster=which(abs(res$p-Duplicated/Tot_CN*purity)<cutoff)
        Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff) 
      }
      if(length(Major_cluster)>0&length(Minor_cluster)>0){
        if(any(Major_cluster%in%Minor_cluster))Major_cluster=Major_cluster[!Major_cluster%in%Minor_cluster]
        Prop_Duplicated=res$prop[Major_cluster]
        Prop_NonDuplicated=sum(res$prop[Minor_cluster])
        conf.intv=poisson.test(c(round(sum(mut_hrd)*Prop_NonDuplicated),
                                 round(sum(mut_hrd)*Prop_Duplicated)))$conf.int
        
        #Prop_NonDuplicated=sum(mut_hrd&maf_patient$multiplicity==1)
        #Prop_Duplicated=sum(mut_hrd&maf_patient$multiplicity==2)
        
        Time=Tot_CN/(Duplicated+Prop_NonDuplicated/Prop_Duplicated)
        # conf.intv=poisson.test(c(round(sum(mut_hrd)*Prop_NonDuplicated),
        #                          round(sum(mut_hrd)*Prop_Duplicated)))$conf.int
        # conf.intv=poisson.test(c(Prop_NonDuplicated,Prop_Duplicated))$conf.int
        
        Conf.Time1=Tot_CN/(Duplicated+conf.intv[1])
        Conf.Time2=Tot_CN/(Duplicated+conf.intv[2])
        annot$HRD_timing[num]=Time
        annot$HRD_timing_low[num]=Conf.Time2
        annot$HRD_timing_high[num]=Conf.Time1
        
        if(!is.null(timing_hrd_chrs)){
          annot$HRD_timing_chrs[num]=paste(round(timing_hrd_chrs,digits = 3),collapse=",")
          annot$HRD_timing_low_chrs[num]=paste(round(timing_hrd_chrs_low,digits = 3),collapse=",")
          annot$HRD_timing_high_chrs[num]=paste(round(timing_hrd_chrs_high,digits = 3),collapse=",")
        }
        #annot$num_m1[num]=paste(num_m1,collapse=",")
        #annot$num_m2[num]=paste(num_m2,collapse=",")
        annot$num_m1[num]=round(sum(mut_hrd)*Prop_NonDuplicated)
        annot$num_m2[num]=round(sum(mut_hrd)*Prop_Duplicated)
        annot$purity[num]=purity
        
      }else{
        annot$HRD_timing[num]="NO CLUSTERING"
        annot[num,c("HRD_timing_low","HRD_timing_high","HRD_chrs",
                    "HRD_timing_chrs","HRD_timing_low_chrs","HRD_timing_high_chrs","num_m1","num_m2","purity")]=NA
      }
    }else{
      annot$HRD_timing[num]="NO MUTS"
      annot[num,c("HRD_timing_low","HRD_timing_high","HRD_chrs",
                  "HRD_timing_chrs","HRD_timing_low_chrs","HRD_timing_high_chrs","num_m1","num_m2","purity")]=NA
    }
    }else{
      annot$HRD_timing[num]="NO SIMPLE GAINS"
      
      annot[num,c("HRD_timing_low","HRD_timing_high","HRD_chrs",
                  "HRD_timing_chrs","HRD_timing_low_chrs","HRD_timing_high_chrs","num_m1","num_m2","purity")]=NA
    }
  }else{
    annot$HRD_timing[num]="NO GAINS"
    annot[num,c("HRD_timing_low","HRD_timing_high","HRD_chrs",
                "HRD_timing_chrs","HRD_timing_low_chrs","HRD_timing_high_chrs","num_m1","num_m2","purity")]=NA
  }
}
write.table(annot,"MM_HRD_timing_all.txt",quote=F,row.names = F,sep="\t")

annot_flt=annot[!grepl("NO",annot$HRD_timing),]
annot_flt$HRD_timing=as.numeric(annot_flt$HRD_timing)
annot_flt$Real_time=annot_flt$HRD_timing*annot_flt$Age
annot_flt$Real_time_low=annot_flt$HRD_timing_low*annot_flt$Age
annot_flt$Real_time_high=annot_flt$HRD_timing_high*annot_flt$Age

annot_flt$Stage=factor(annot_flt$Stage, levels = c("MGUS", "LRSMM", "IRSMM", "HRSMM", "NDMM", "MM"))
annot_flt=annot_flt[order(annot_flt$Stage,annot_flt$Age),]

pdf("purity_vs_timing.pdf")
plot(annot_flt$HRD_timing,annot_flt$purity,xlab="Timing",ylab="Purity",pch=21,bg="steelblue")
dev.off()
#annot_flt=annot[!is.na(annot$Real_time),]
pdf("absolute_timing_2024_02_24.pdf",width=15,height=6,useDingbats = F)
plot(annot_flt$Age,ylim=c(0,90),pch=19,xaxt='n',xlab="",ylab="Age")
points(annot_flt$Real_time,pch=19,col='steelblue')
segments(x0=1:nrow(annot_flt),y0=annot_flt$Real_time_low,y1=annot_flt$Real_time_high, col='steelblue')     
abline(v=cumsum(table(annot_flt$Stage))+0.5,lty='dashed')
text(y = 90, x= c(cumsum(table(annot_flt$Stage))+c(0,cumsum(table(annot_flt$Stage))[-6]))/2+0.5,labels=unique(annot_flt$Stage))
dev.off()

annot$HRD_timing_final=annot$HRD_timing
annot$HRD_timing_low_final=annot$HRD_timing_low
annot$HRD_timing_high_final=annot$HRD_timing_high

for(donor in rownames(annot)[!is.na(annot$HRD_timing)]){
  num=which(rownames(annot)==donor)
  
  chrs=unlist(strsplit(annot$HRD_chrs[num],split=','))
  ord=order(as.numeric(chrs))
  
  dup=as.numeric(unlist(strsplit(annot$num_m2[num],split=',')))
  nondup=as.numeric(unlist(strsplit(annot$num_m1[num],split=',')))
  is_diff=rep(FALSE,length(chrs))
  if(length(chrs)>1){
    for(n in 1:length(chrs)){
      is_diff[n]=poisson.test(c(nondup[n],dup[n]),c(sum(nondup[-n]),sum(dup[-n])))$p.value<0.01
    }
  }

  if(any(is_diff)){
    new_dup=sum(dup[!is_diff])
    new_nondup=sum(nondup[!is_diff])
    new_time=Tot_CN/(Duplicated+new_nondup/new_dup)
  
    conf.intv=poisson.test(c(new_nondup,new_dup))$conf.int
    new_high=Tot_CN/(Duplicated+conf.intv[1])
    new_low=Tot_CN/(Duplicated+conf.intv[2])
    pdf(paste0(donor,"_HRD.pdf"))
    plot(c(annot$HRD_timing[num],new_time,unlist(strsplit(annot$HRD_timing_chrs[num],split=','))[ord]),
         pch=19,col=c('black',"firebrick",rep('steelblue',length(chrs))),ylim=c(-0.05,1.25),main=donor,
         ylab="Relative time",cex=1.5,xaxt='n',xlab="")
    is_diff_txt=rep("",length(chrs))
    is_diff_txt[is_diff]='*'
    
    segments(x0=1:(length(chrs)+2),
             y0=as.numeric(c(annot$HRD_timing_low[num],new_low,unlist(strsplit(annot$HRD_timing_low_chrs[num],split=','))[ord])),
             y1=as.numeric(c(annot$HRD_timing_high[num],new_high,unlist(strsplit(annot$HRD_timing_high_chrs[num],split=','))[ord])), 
             col=c('black','firebrick',rep('steelblue',length(chrs))),lwd=2)   
    abline(h=1,lty='dashed')
    text(x = 1:(length(chrs)+2), y= 0,labels=c("HRD","HRD_corrected",paste0("chr",chrs[ord],is_diff_txt[ord])))
    dev.off()
    annot$HRD_timing_final[num]=new_time
    annot$HRD_timing_low_final[num]=new_low
    annot$HRD_timing_high_final[num]=new_high
    
    
  }else{
    pdf(paste0(donor,"_HRD.pdf"))
    plot(c(annot$HRD_timing[num],unlist(strsplit(annot$HRD_timing_chrs[num],split=','))[ord]),
         pch=19,col=c('black',rep('steelblue',length(chrs))),ylim=c(-0.05,1.25),main=donor,
         ylab="Relative time",cex=1.5,xaxt='n',xlab="")
  
    
    segments(x0=1:(length(chrs)+1),
             y0=as.numeric(c(annot$HRD_timing_low[num],unlist(strsplit(annot$HRD_timing_low_chrs[num],split=','))[ord])),
             y1=as.numeric(c(annot$HRD_timing_high[num],unlist(strsplit(annot$HRD_timing_high_chrs[num],split=','))[ord])), 
             col=c('black',rep('steelblue',length(chrs))),lwd=2)   
    abline(h=1,lty='dashed')
    text(x = 1:(length(chrs)+1), y= 0,labels=c("HRD",paste0("chr",chrs[ord])))
    dev.off()
    
  }
}

annot$Real_time=annot$HRD_timing_final*annot$Age
annot$Real_time_low=annot$HRD_timing_low_final*annot$Age
annot$Real_time_high=annot$HRD_timing_high_final*annot$Age


annot_flt=annot[!is.na(annot$Real_time),]
annot_flt=annot_flt[order(annot_flt$Stage,annot_flt$HRD_timing_final),]


pdf("relative_timing.pdf",width=15,height=6,useDingbats = F)
plot(annot_flt$HRD_timing,ylim=c(0,1.15),pch=19,xaxt='n',xlab="",ylab="Relative time",col='steelblue')
segments(x0=1:nrow(annot_flt),y0=annot_flt$HRD_timing_low_final,y1=annot_flt$HRD_timing_high_final, col='steelblue')     
abline(v=cumsum(table(annot_flt$Stage))+0.5,lty='dashed')
text(y = 0, x= c(cumsum(table(annot_flt$Stage))+c(0,cumsum(table(annot_flt$Stage))[-6]))/2+0.5,labels=unique(annot_flt$Stage))
abline(h=1,col='firebrick',lty='dashed')
text(y = c(0,1.05), x= 0,labels=c("Zygote","Most recent clonal ancestor"),pos=4)
dev.off()

annot_flt=annot_flt[order(annot_flt$Stage,annot_flt$Age),]

pdf("age_w_relative_timing.pdf",width=15,height=6,useDingbats = F)
plot(annot_flt$Age,ylim=c(0,90),pch=19,xaxt='n',xlab="",ylab="Age")
points(annot_flt$Real_time,pch=19,col='steelblue')
segments(x0=1:nrow(annot_flt),y0=annot_flt$Real_time_low,y1=annot_flt$Real_time_high, col='steelblue')     
abline(v=cumsum(table(annot_flt$Stage))+0.5,lty='dashed')
text(y = 5, x= c(cumsum(table(annot_flt$Stage))+c(0,cumsum(table(annot_flt$Stage))[-6]))/2+0.5,labels=unique(annot_flt$Stage))
dev.off()


annot_flt=annot_flt[order(annot_flt$Stage,annot_flt$Real_time),]
# 
# pdf("age_w_relative_timing_order_hrd.pdf",width=15,height=6,useDingbats = F)
# plot(annot_flt$Age,ylim=c(0,90),pch=19,xaxt='n',xlab="",ylab="Age")
# points(annot_flt$Real_time,pch=19,col='steelblue')
# segments(x0=1:nrow(annot_flt),y0=annot_flt$Real_time_low,y1=annot_flt$Real_time_high, col='steelblue')     
# abline(v=cumsum(table(annot_flt$Stage))+0.5,lty='dashed')
# text(y = 5, x= c(cumsum(table(annot_flt$Stage))+c(0,cumsum(table(annot_flt$Stage))[-6]))/2+0.5,labels=unique(annot_flt$Stage))
# dev.off()

plot(pmin(1,annot_flt$HRD_timing))

wilcox.test(annot_flt$Real_time[annot$Stage=="MGUS"],
            annot_flt$Real_time[annot$Stage!="MGUS"])
nonhrd_time=c()
for(chr in nonhrd_gains){
  mut_cnv=maf_patient$Chromosome==chr&maf_patient$clonal==1&!is.na(maf_patient$clonal)&maf_patient$ml_sig=="SBS5"
  
  NV_vec=maf_patient$t_alt_count[mut_cnv]
  NR_vec=maf_patient$t_alt_count[mut_cnv]+maf_patient$t_ref_count[mut_cnv]
  res = binom_mix(NV_vec,NR_vec,nrange=1:3,mode='Truncated')
  
  
  Major_cluster=which(abs(res$p-Duplicated/Tot_CN*purity)<cutoff)
  Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff)
  Prop_Duplicated=res$prop[Major_cluster]
  Prop_NonDuplicated=sum(res$prop[Minor_cluster])
  Time=Tot_CN/(Duplicated+Prop_NonDuplicated/Prop_Duplicated)
  conf.intv=poisson.test(c(round(sum(mut_cnv)*Prop_NonDuplicated),
                           round(sum(mut_cnv)*Prop_Duplicated)))$conf.int
  Conf.Time1=Tot_CN/(Duplicated+conf.intv[1])
  Conf.Time2=Tot_CN/(Duplicated+conf.intv[2])
  nonhrd_time=rbind(nonhrd_time,c(chr,Prop_Duplicated,Prop_NonDuplicated,Time,Conf.Time1,Conf.Time2))
}



p=hist(NV_vec/NR_vec,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
       main=paste0(donor,", (n=",length(NV_vec),")"))
cols=c("red","blue","green","magenta","cyan")

y_coord=max(p$density)-0.5
y_intv=y_coord/5

for (i in 1:res$n){
  depth=rpois(n=5000,lambda=median(NR_vec))
  sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
  sim_VAF=sim_NV/depth
  sim_VAF=sim_VAF[sim_NV>3]
  dens=density(sim_VAF)
  lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
  y_coord=y_coord-y_intv/2
  text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2)))
  segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
}

#------
# Maura's data

caveman=rbind(read.table("Downloads/Project1212Caveman.txt",header=T,sep="\t"),
              read.table("Downloads/ChapmanCaveman.txt",header=T,sep="\t")[,1:64])
purity_df=read.table("Downloads/ccf_samples.txt",header=T,sep="\t")
caveman$NV=0
caveman$NV[caveman$Alt=="A"]=caveman$FAZ.Tum[caveman$Alt=="A"]+caveman$RAZ.Tum[caveman$Alt=="A"]
caveman$NV[caveman$Alt=="C"]=caveman$FCZ.Tum[caveman$Alt=="C"]+caveman$RCZ.Tum[caveman$Alt=="C"]
caveman$NV[caveman$Alt=="G"]=caveman$FGZ.Tum[caveman$Alt=="G"]+caveman$RGZ.Tum[caveman$Alt=="G"]
caveman$NV[caveman$Alt=="T"]=caveman$FTZ.Tum[caveman$Alt=="T"]+caveman$RTZ.Tum[caveman$Alt=="T"]
caveman$NR=rowSums(caveman[,grepl("Z.Tum",colnames(caveman))])
battenberg=read.table("Downloads/CNV_batte.txt",header=T,sep="\t")
hrd_chrs_poss=c(3,5,7,9,11,15,19,21)

battenberg=battenberg[battenberg$nMaj1_A>1&battenberg$Chrom%in%hrd_chrs_poss&battenberg$frac1_A==1,]


samples=purity_df$Sample

battenberg$length=battenberg$end-battenberg$start
annot=data.frame(Sample=samples,
                 HRD_timing=0,
                 HRD_timing_low=0,
                 HRD_timing_high=0,
                 num_muts=0,
                 prop_duplicated=0,
                 prop_nonduplicated=0,
                 purity=0)
for(n in 1:length(samples)){
  sample=samples[n]
  caveman_sub=caveman[caveman$Sample==sample,]
  battenberg_sub=battenberg[battenberg$sample==sample&battenberg$nMaj1_A==2&battenberg$nMin1_A==1,]
  hrd_chr=unique(battenberg_sub$Chrom[battenberg_sub$length>1e7])
  battenberg_sub=battenberg_sub[battenberg_sub$Chrom%in%hrd_chr,]
  if(nrow(battenberg_sub)>0){
    mut_hrd=rep(F,nrow(caveman_sub))
    for(k in 1:nrow(battenberg_sub)){
      mut_hrd=(caveman_sub$Chrom==battenberg_sub$Chrom[k]&
                 caveman_sub$Pos>battenberg_sub$start[k]&
                 caveman_sub$Pos<battenberg_sub$end[k])|mut_hrd
    }
    if(sum(mut_hrd)>0){
      NR_vec=caveman_sub$NR[mut_hrd]
      NV_vec=caveman_sub$NV[mut_hrd]
      
      res = binom_mix(NV_vec,NR_vec,nrange=1:4,mode='Truncated',criterion="AIC")
      
      Duplicated=2
      Tot_CN=3
      purity=purity_df$CCF[purity_df$Sample==sample]
      if(sample=="PD26411a") purity=0.85
      if(sample=="PD26411c") purity=0.7
      if(sample=="PD26414a") purity=0.5
      if(sample=="PD26416d") purity=0.75
      
      cutoff=0.08
      Major_cluster=which(abs(res$p-Duplicated/Tot_CN*purity)<cutoff)
      Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff)
      Prop_Duplicated=res$prop[Major_cluster]
      Prop_NonDuplicated=sum(res$prop[Minor_cluster])
      if(length(Prop_Duplicated)>0&length(Prop_NonDuplicated)>0&Prop_NonDuplicated>0){
        
        Time=Tot_CN/(Duplicated+Prop_NonDuplicated/Prop_Duplicated)
        conf.intv=poisson.test(c(round(sum(mut_hrd)*Prop_NonDuplicated),
                                  round(sum(mut_hrd)*Prop_Duplicated)))$conf.int
        
        Conf.Time1=Tot_CN/(Duplicated+conf.intv[1])
        Conf.Time2=Tot_CN/(Duplicated+conf.intv[2])
        annot$HRD_timing[n]=Time
        annot$HRD_timing_low[n]=Conf.Time2
        annot$HRD_timing_high[n]=Conf.Time1
        annot$prop_nonduplicated[n]=Prop_NonDuplicated
        annot$prop_duplicated[n]=Prop_Duplicated
        annot$num_muts[n]=sum(mut_hrd)
        annot$purity[n]=purity
      }else{
        annot$HRD_timing[n]="NO CLUSTERING"
      }
    }else{
      annot$HRD_timing[n]="NO MUTS"
    }
  }else{
    annot$HRD_timing[n]="NO HRD"
  }
}
table(annot$HRD_timing[grepl("NO",annot$HRD_timing)])
annot[annot$HRD_timing=="NO CLUSTERING",]
annot_flt=annot[!grepl("NO",annot$HRD_timing),]
plot(annot_flt$purity,annot_flt$HRD_timing)

annot_flt=annot[!grepl("NO",annot$HRD_timing),]

write.table(annot,"Francesco_data_MM_timing.txt",quote=F,row.names = F,sep="\t")

plot(annot_flt$HRD_timing,ylim=c(0,1.15),pch=19,xaxt='n',xlab="",ylab="Relative time",col='steelblue')
segments(x0=1:nrow(annot_flt),y0=annot_flt$HRD_timing_low_final,y1=annot_flt$HRD_timing_high_final, col='steelblue')     
abline(v=cumsum(table(annot_flt$Stage))+0.5,lty='dashed')
text(y = 0, x= c(cumsum(table(annot_flt$Stage))+c(0,cumsum(table(annot_flt$Stage))[-6]))/2+0.5,labels=unique(annot_flt$Stage))
abline(h=1,col='firebrick',lty='dashed')
text(y = c(0,1.05), x= 0,labels=c("Zygote","Most recent clonal ancestor"),pos=4)
