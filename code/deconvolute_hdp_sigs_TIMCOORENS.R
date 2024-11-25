#Deconvolute (HDP) mutational signatures into reference signatures and arrive at final set
options(stringsAsFactors = F)
library(hdp)
library(RColorBrewer)
library(lsa)
library(lattice)

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

#Load HDP signatures
hdp_sigs=read.table("hdp_sigs.txt")
colnames(hdp_sigs)=gsub("X","N",colnames(hdp_sigs))

#Load reference signatures, make sure they are ordered correctly and actually sum to 1
ref=read.csv("/lustre/scratch117/casm/team268/tc16/COSMIC_Mutational_Signatures_v3.1.csv")
rownames(ref)=paste0(substr(ref$Subtype,1,1),"[",ref$Type,"]",substr(ref$Subtype,3,3))
ref=ref[,-c(1,2)]
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
full_vec = paste0(rep(c("A","C","G","T"),each=4),"[",rep(sub_vec,each=16),"]",rep(c("A","C","G","T"),times=4))
ref=ref[full_vec,]
ref=apply(ref,2,as.numeric)
ref[is.na(ref)|ref==0]=0.00001
ref=t(t(ref)/colSums(ref))

#Assess cosine similarities for all reference signatures
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, "Cosine_similarities.txt",sep="\t",quote=F)

#plot output
pdf("cosine_similarities.pdf", height=5, width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

## The actuasl decomposition

#First iteration; decomposed hdp sigs into all suspected sigs 

#Pick your candidate signatures
gdsigs=c("SBS1","SBS2","SBS5","SBS6","SBS7a","SBS7b","SBS7d","SBS8","SBS13","SBS16", "SBS17b", "SBS18","SBS22","SBS23","SBS32","SBS35","SBS40","SBS48","SBS52","SBS88") 

add="v1" #add text to output files

signatures=t(ref[,gdsigs])
sample_list=paste0("N",c(0:10)) #The list of signatures to be deconvolved
colnames(hdp_sigs)=paste0("N",c(0:10))
profiles=hdp_sigs[,sample_list]


#The First Iteration
signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
rownames(signature_fraction) = rownames(signatures)
colnames(signature_fraction) = sample_list
maxiter <- 1000

for (j in 1:length(sample_list)) {
  freqs = profiles[,j]
  freqs[is.na(freqs)] = 0
  # EM algowith to estimate the signature contribution
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

#Plot initial deconvolution and save results
pdf(paste0("Deconvolution_hdp_sigs_R1_",add,"_priors.pdf"), height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fraction[nrow(signature_fraction):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(signature_fraction, paste0("hdp_known_sigs_broken_down_into_pcawg_gd_sigs_",add,"_priors.txt"), sep="\t", col.names=T, row.names = T, quote=F)


#Only keep candidate signature above a certain contribution and do a second round of deconvolution
tresh=0.1
sigs_deconv_R2=list()
for(n in 1:length(sample_list)){
  sigs_deconv_R2[[n]]=rownames(signature_fraction)[signature_fraction[,n]>tresh]
}
names(sigs_deconv_R2)=colnames(signature_fraction)

#Some signatures only have one major contributor - leaving those out
sigs_to_deconv=names(sigs_deconv_R2)[unlist(lapply(sigs_deconv_R2,length))>1]

ref_sigs_R2=sort(unique(unlist(sigs_deconv_R2)))
signature_fractionR2=matrix(NA,ncol=length(sigs_to_deconv),nrow=length(ref_sigs_R2))
rownames(signature_fractionR2)=ref_sigs_R2
colnames(signature_fractionR2)=sigs_to_deconv
#repeat the deconvolution with the identified constitutive signatures
n=1
for(s in sigs_to_deconv){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])
  
  signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
  rownames(signature_fraction) = rownames(signatures)
  colnames(signature_fraction) = sample_list
  maxiter <- 1000
  
  freqs = profiles[,s]
  freqs[is.na(freqs)] = 0
  
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  signature_fractionR2[gdsigs,n] = alpha
  n=n+1
  reconsbs <- rep(0,96)
  for (g in gdsigs) {
    reconsbs=reconsbs+(ref[,g]*alpha[g])
  }
  cosine_reconst=cosine(x=reconsbs, y=hdp_sigs[,s])
  print(paste0(s,": ",cosine_reconst))
  pdf(paste0("HDP_",s,"_reconstitution_",add,"_priors.pdf"))
  par(mfrow=c(length(alpha)+2,1))
  par(mar=c(1,2,4,1))
  barplot(hdp_sigs[,s], col=mut.cols, main=paste0("HDP ",s),names.arg="")
  barplot(reconsbs, col=mut.cols, main=paste0("Reconstituted ",s," cosine similarity to original: ", round(cosine_reconst, digits=2)))
  for (g in gdsigs) {
    barplot(ref[,g], col=mut.cols, main=paste0("PCAWG ", g, " accounts for ", round(alpha[g], digits=2)))
  }
  dev.off()
}

#If there are any novel signatures in the set, note them here with a name you like, such as:
sigs_deconv_R2$N14="SBS91"

saveRDS(sigs_deconv_R2,"hdp2refsigs.Rdata")


#Combine hdp signatures that did not get deconvolved and reference signatures into final table
final_sigs=cbind(hdp_sigs[,!colnames(hdp_sigs)%in%names(sigs_deconv_R2)],ref[,ref_sigs_R2])
#Rename the HDP components that didn't get deconvoluted, if necessary
colnames(final_sigs)[1:4]=paste0("SBS10",LETTERS[2:5])
colnames(final_sigs)[1]="SBS35like"
#Order them - optional
final_sigs=final_sigs[,paste0("SBS",c(1,5,"7a","7b","10a","10B","10C","10D","10E","17b",28,"35like","A","B"))]
final_sigs=ref[,sort(unique(unlist(sigs_deconv_R2)))]
write.table(final_sigs,"final_sigs_2020_08_31.txt")
