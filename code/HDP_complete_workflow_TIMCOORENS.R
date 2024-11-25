## HDP Flow
# Tim Coorens, 2020

# PART 1: Generate input matrix
# This part of the code generates a matrix of trinucleotide mutation counts
# Can be skipped if already available

options(stringsAsFactors = F)
library(data.table)

patients=read.table("patients.txt")[,1]
patients=patients[!patients%in%c("PD40097","PD44583","PD44594")]
all_muts=c()
for (p in patients){
  if(file.exists(paste0(p,"/",p,".snp_assigned_to_branches.txt"))){
    all_muts=rbind(all_muts,fread(paste0(p,"/",p,".snp_assigned_to_branches.txt"),header=T,sep="\t",data.table=F))
    print(p)
  }
}

mutlist_to_96_contexts=function(mutlist,version){
  if(version==37) genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
  if(version==38) genomeFile = "/lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  samples=unique(mutlist$SampleID)
  trinuc_mut_mat=matrix(0,ncol=96,nrow=length(samples))
  for (n in 1:length(samples)){
    s=samples[n]
    mutations=as.data.frame(mutlist[mutlist$SampleID==s,c("Chr","Pos","Ref","Alt")])
    colnames(mutations) = c("chr","pos","ref","mut")
    mutations$pos=as.numeric(mutations$pos)
    if(version==37) mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
    if(version==38) utations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% paste0("chr",c(1:22,"X","Y")),]
    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                       as.numeric(mutations$pos)+1))))
    ntcomp = c(T="A",G="C",C="G",A="T")
    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
    mutations$trinuc_ref_py = mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A","G")) { # Purine base
        mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
        mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    trinuc_mut_mat[n,]=freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat)=full_vec
  rownames(trinuc_mut_mat)=samples
  return(trinuc_mut_mat)
}

# Start with a list of mutations per sample (or per branch, if phylogeny)
# make sure colnames are Chr Pos Ref Alt SampleID (order shouldn't matter)

## If required/desired, subset max. mutations per sample here
freqs=table(all_muts$SampleID)
samples=names(freqs)
max_mut_num=3000
mut_list_subsampled=c()
for(s in samples){
  mut_list_sub=all_muts[all_muts$SampleID==s,]
  if(nrow(mut_list_sub)<=max_mut_num){
    mut_list_subsampled=rbind(mut_list_subsampled,mut_list_sub)
  }else{
    random_select=sample(size=max_mut_num,x=1:nrow(mut_list_sub),replace=F)
    mut_list_subsampled=rbind(mut_list_subsampled,mut_list_sub[random_select,])
  }
  print(s)
}
##

trinuc_mut_mat=mutlist_to_96_contexts(all_muts[,c(1:4,7)])
#if subsampled:
trinuc_mut_mat=mutlist_to_96_contexts(mut_list_subsampled[,c(1:4,7)],version=38)


samples=rownames(trinuc_mut_mat)
key_table=data.frame(Sample=samples,
                     Patient=substr(samples,1,7))
write.table(trinuc_mut_mat,"trinuc_mut_mat_full.txt")
write.table(key_table,"key_table_full.txt")

## PART 2: Run multiple HDP chains

# I run this using the following simple bash script:

#for n in $(seq 1 20);
#do
#bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R'select[mem>10000] rusage[mem=10000]' -M10000 -J $n Rscript hdp_single_chain.R $n
#done


###
# hdp_single_chain.R

options(stringsAsFactors = F)
library(hdp)
lower_threshold=100

n=as.numeric(commandArgs(T)[1])
mutations=read.table("trinuc_mut_mat.txt")
key_table=read.table("key_table.txt")

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

#Hierarchy is set per patient, can change if wanted
freq=table(key_table$Patient)

hdp_mut <- hdp_init(ppindex = c(0, rep(1,length(freq)),rep(2:(length(freq)+1), times=freq)), # index of parental node
                    cpindex = c(1, rep(2,length(freq)),rep(3:(length(freq)+2), times=freq)), # index of the CP to use
                    hh = rep(1, 96), # prior is uniform over 96 categories
                    alphaa = rep(1,length(freq)+2), # shape hyperparameters for 2 CPs
                    alphab = rep(1,length(freq)+2))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut, 
                       dpindex = (length(freq)+2):numdp(hdp_mut), # index of nodes to add data to
                       mutations)

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10,seed=n*300)

chain=hdp_posterior(hdp_activated,
                    burnin=20000,
                    n=100,
                    seed=n*1000,
                    space=200,
                    cpiter=3)
saveRDS(chain,paste0("hdp_chain_",n,".Rdata"))

### PART 3: combine the results
options(stringsAsFactors = F)
library(hdp)

chlist <- vector("list", 20)
for (i in 1:20){
  if(file.exists(paste0("hdp_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

mut_example_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi) #This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi,"HDP_multi_chain.Rdata")

pdf("muts_attributed.pdf")
plot_comp_size(mut_example_multi, bty="L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

#dev.new(width=12,height=4)
#par(mfrow=c(3,4))


for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_",i,".pdf"),width=12,height=4)
  
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"),
                      incl_nonsig=TRUE, cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')

pdf("signature_attribution.pdf",width=10,height=8)

plot_dp_comp_exposure(mut_example_multi, dpindices=(length(mut_example_multi@comp_dp_counts)-nrow(counts)+1):length(mut_example_multi@comp_dp_counts), incl_nonsig = T,ylab_exp = 'Signature exposure', leg.title = 'Signature',
                      col=c(RColorBrewer::brewer.pal(12, "Set3"),"magenta","firebrick",RColorBrewer::brewer.pal(8, "Dark2")))
dev.off()

mutations=read.table("../trinuc_mut_mat_full.txt")
key_table=read.table("../key_table_full.txt")

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

freq=table(key_table$Patient)

dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
# FIXME JB who is mutations here? what if length(freq)+1+1 greater than nrow(mutations)? 
# Why not similar to above?
# length(mut_example_multi@comp_dp_counts)-nrow(counts)+1):length(mut_example_multi@comp_dp_counts)

mean_assignment <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),,drop=FALSE])
mean_assignment <- t(dp_distn["mean"][length(freq)+1+1:nrow(mutations),,drop=FALSE])
write.table(mean_assignment,"mean_assignment_hdp.txt")

mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
lower=mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))

write.table(mean_sigs,"hdp_sigs.txt")

#PART 4: Deconvolute HDP signatures into reference signatures and arrive at final set
options(stringsAsFactors = F)
library(hdp)
library(RColorBrewer)
library(lsa)
library(lattice)

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

#Load HDP signatures
hdp_sigs=read.table("hdp_sigs.txt")
#Load reference signatures
ref=read.csv("/lustre/scratch117/casm/team268/tc16/COSMIC_Mutational_Signatures_v3.1.csv")
ref=read.table("../annot/cosmic/COSMIC_v3.3.1_SBS_GRCh37.txt", header=TRUE)
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

colnames(hdp_sigs)=gsub("X","N",colnames(hdp_sigs))

#First iteration; decomposed hdp sigs into all suspected sigs 

#Make selection of priors sigs - alternatively, use all COSMIC sigs, but this leads to overfitting
gdsigs=c("SBS1","SBS2","SBS4","SBS5","SBS7a","SBS7b","SBS13","SBS16", "SBS17b", "SBS18","SBS22","SBS32","SBS35","SBS40","SBS41","SBS88")

add="" #add something to titles to differentiate multiple runs [OPTIONAL]

signatures=t(ref[,gdsigs])
sample_list=paste0("N",c(0:(ncol(hdp_sigs)-1))) 
colnames(hdp_sigs)=paste0("N",c(0:(ncol(hdp_sigs)-1))) 
profiles=hdp_sigs[,sample_list]

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

sigs_deconv_R2=list()
for(n in 1:length(sample_list)){
  sigs_deconv_R2[[n]]=rownames(signature_fraction)[signature_fraction[,n]>0.15]
}
names(sigs_deconv_R2)=colnames(signature_fraction)


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

sigs_deconv_R2$N11="N11"
sigs_deconv_R2$N14="SBS91"

saveRDS(sigs_deconv_R2,"hdp2refsigs.Rdata")

#Combine hdp signatures that did not get deconvolved and reference signatures into final table
final_sigs=cbind(hdp_sigs[,"N14"],ref[,ref_sigs_R2])
#Rename the HDP components that didn't get deconvoluted
write.table(final_sigs,"final_sigs_2020_08_31.txt")
  

#PART 5: Fit signatures to observed counts
options(stringsAsFactors = F)
library(sigfit)
library(RColorBrewer)
#library(ape)
#library(ggtree)
data("cosmic_signatures_v3")

final_sigs=t(read.table("final_sigs_2020_08_31.txt"))
all_counts=read.csv("/lustre/scratch119/casm/team267ms/lm14/panbody_map/signatures/HDP/001_Input/mut_count.panbodyTestes_trees_11082020.csv")
rownames(all_counts)=all_counts$X
all_counts=all_counts[,-1]
hdp_exposures=read.table("mean_assignment_hdp.txt")
hdp_sigs=read.table("hdp_sigs.txt")
colnames(hdp_sigs)=gsub("X","N",colnames(hdp_sigs))
colnames(hdp_exposures)=gsub("X","N",colnames(hdp_exposures))

key_table=read.table("/lustre/scratch117/casm/team268/tc16/PanBody/Full_cohort_key.txt",header=T)

colnames(hdp_exposures)=paste0("N",0:14)
hdp2ref=readRDS("hdp2refsigs.Rdata")
hdp2ref$N11="N11"
hdp2ref$N14="N14"

data("cosmic_signatures_v3")

colnames(all_counts)=colnames(final_sigs)=colnames(cosmic_signatures_v3)
final_sigs=rbind(final_sigs,t(hdp_sigs[,c("N11","N14")]))
tissues=unique(key_table$TissueType4)

patients=unique(substr(rownames(all_counts),1,7))

conv=c()
for(p in c("PD28690","PD43850","PD43851")){
  data=read.table(paste0("/lustre/scratch117/casm/team268/tc16/PanBody/",p,"/",p,"_snp_assigned_to_branches_full_IDs.txt"),header=T)
  conv=rbind(conv,unique(data[,c("SampleID","Sample.ID")]))
}
full=conv$Sample.ID
names(full)=conv$SampleID
rownames(all_counts)[rownames(all_counts)%in%conv$SampleID]=full[rownames(all_counts)[rownames(all_counts)%in%conv$SampleID]]
rownames(hdp_exposures)[rownames(hdp_exposures)%in%conv$SampleID]=full[rownames(hdp_exposures)[rownames(hdp_exposures)%in%conv$SampleID]]

conv=c()
for(p in patients[!patients%in%c("PD28690","PD43850","PD43851")]){
  data=read.table(paste0("/lustre/scratch117/casm/team267/tc16/PanBody/Testes_calls/",p,"_snp_assigned_to_branches_full_IDs.txt"),header=T)
  conv=rbind(conv,unique(data[,c("SampleID","Sample.ID")]))
}
full=conv$Sample.ID
names(full)=conv$SampleID
rownames(all_counts)[rownames(all_counts)%in%conv$SampleID]=full[rownames(all_counts)[rownames(all_counts)%in%conv$SampleID]]
rownames(hdp_exposures)[rownames(hdp_exposures)%in%conv$SampleID]=full[rownames(hdp_exposures)[rownames(hdp_exposures)%in%conv$SampleID]]

sigs_per_tissue_matrix=matrix(0,ncol=length(unique(unlist(hdp2ref))),nrow=57)
colnames(sigs_per_tissue_matrix)=unique(unlist(hdp2ref))
n=0
rownames_mat=c()
sigs_per_tissue=list()

exp_thresh=0.075
min_sample_num=2
for(patient in patients){
  tissues=unique(key_table$TissueType4[key_table$DonorID==patient])
  sigs_per_tissue[[patient]]=list()
  for(tissue in tissues){
    samples_tissue=key_table$SampleID[key_table$TissueType4==tissue&key_table$DonorID==patient]
    select=grepl(paste(samples_tissue,collapse="|"),rownames(hdp_exposures))
    if(sum(select)){
      n=n+1
      hdp_sigs_present=colnames(hdp_exposures)[colSums(hdp_exposures[select,]>exp_thresh)>min(min_sample_num,sum(select)-1)]
      ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present])))
      sigs_per_tissue[[patient]][[tissue]]=ref_sigs_present
      sigs_per_tissue_matrix[n,ref_sigs_present]=1
      rownames_mat=c(rownames_mat,paste0(patient,"_",tissue))
    }
  }
}

rownames(sigs_per_tissue_matrix)=rownames_mat
pdf(paste0("heatmap_sigs_",exp_thresh*100,"pc_",min_sample_num,".pdf"),width=12,height=10)
heatmap(sigs_per_tissue_matrix,scale='none',cexRow = 0.6,margins = c(5,13))
dev.off()
write.table(sigs_per_tissue_matrix,paste0("sigs_",exp_thresh*100,"pc_",min_sample_num,".txt"))


saveRDS(sigs_per_tissue,"sigs_per_tissue.Rdata")

sample_counts=read.table("/lustre/scratch117/casm/team268/tc16/PanBody/All_Matched.SBS96.all",header=T)
rownames(sample_counts)=sample_counts$MutationType
sample_counts=sample_counts[,-1]
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
full_vec = paste0(rep(c("A","C","G","T"),each=4),"[",rep(sub_vec,each=16),"]",rep(c("A","C","G","T"),times=4))
sample_counts=sample_counts[full_vec,]

sigs_fitting=unique(unlist(sigs_per_tissue))
sample_signatures=sample_signatures_lower=sample_signatures_upper=matrix(0,ncol=length(sigs_fitting),nrow=ncol(sample_counts))
rownames(sample_signatures)=rownames(sample_signatures_lower)=rownames(sample_signatures_upper)=colnames(sample_counts)
colnames(sample_signatures)=colnames(sample_signatures_lower)=colnames(sample_signatures_upper)=sigs_fitting

for(patient in patients){
  print(patient)
  tissues=unique(key_table$TissueType4[key_table$DonorID==patient])
  for(tissue in tissues){
    print(tissue)
    ref_sigs_present=sigs_per_tissue[[patient]][[tissue]]
    samples_tissue=key_table$SampleID[key_table$TissueType4==tissue&key_table$DonorID==patient]
    if(!is.null(ref_sigs_present)&file.exists(paste0("pars_objects/",patient,"_",tissue,".Rdata"))&"SBS1"%in%ref_sigs_present){
      pars=readRDS(paste0("pars_objects/",patient,"_",tissue,".Rdata"))
      rownames(pars$mean)=rownames(pars$lower)=rownames(pars$upper)=samples_tissue
      sample_signatures[rownames(pars$mean),colnames(pars$mean)]=as.matrix(pars$mean)
      sample_signatures_lower[rownames(pars$lower),colnames(pars$lower)]=as.matrix(pars$lower)
      sample_signatures_upper[rownames(pars$upper),colnames(pars$upper)]=as.matrix(pars$upper)
    }
    if(!is.null(ref_sigs_present)&file.exists(paste0("pars_objects/",patient,"_",tissue,"_SBS1_added.Rdata"))){
      pars=readRDS(paste0("pars_objects/",patient,"_",tissue,"_SBS1_added.Rdata"))
      rownames(pars$mean)=rownames(pars$lower)=rownames(pars$upper)=samples_tissue
      sample_signatures[rownames(pars$mean),colnames(pars$mean)]=as.matrix(pars$mean)
      sample_signatures_lower[rownames(pars$lower),colnames(pars$lower)]=as.matrix(pars$lower)
      sample_signatures_upper[rownames(pars$upper),colnames(pars$upper)]=as.matrix(pars$upper)
    }
    if(!is.null(ref_sigs_present)&!file.exists(paste0("pars_objects/",patient,"_",tissue,"_SBS1_added.Rdata"))&!"SBS1"%in%ref_sigs_present){
      ref_sigs_present=c(ref_sigs_present,"SBS1")
      counts=t(sample_counts[,samples_tissue])
      colnames(counts)=colnames(cosmic_signatures_v3)
      #ref_sigs_present=names(sort(sig_order[ref_sigs_present]))
      
      fit=fit_signatures(counts=counts,signatures = final_sigs[ref_sigs_present,],
                         iter = 20000,warmup = 10000,model="poisson",chains = 2)
      pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
      rownames(pars$mean)=rownames(pars$lower)=rownames(pars$upper)=samples_tissue
      saveRDS(pars,paste0("pars_objects/",patient,"_",tissue,"_SBS1_added.Rdata"))
      sample_signatures[rownames(pars$mean),colnames(pars$mean)]=as.matrix(pars$mean)
      sample_signatures_lower[rownames(pars$lower),colnames(pars$lower)]=as.matrix(pars$lower)
      sample_signatures_upper[rownames(pars$upper),colnames(pars$upper)]=as.matrix(pars$upper)
    }
  }
}

for(patient in patients){
  print(patient)
  tissues=unique(key_table$TissueType4[key_table$DonorID==patient])
  for(tissue in tissues){
    print(tissue)
    if(!is.null(ref_sigs_present)&file.exists(paste0("pars_objects/",patient,"_",tissue,".Rdata"))){
      pars=readRDS(paste0("pars_objects/",patient,"_",tissue,"_incl_1540.Rdata"))
      sample_signatures[rownames(pars$mean),colnames(pars$mean)]=as.matrix(pars$mean)
      sample_signatures_lower[rownames(pars$lower),colnames(pars$lower)]=as.matrix(pars$lower)
      sample_signatures_upper[rownames(pars$upper),colnames(pars$upper)]=as.matrix(pars$upper)
    }
  }
}

sigs_fitting=unique(unlist(sigs_per_tissue))
sample_signatures=sample_signatures_lower=sample_signatures_upper=matrix(0,ncol=length(sigs_fitting),nrow=ncol(sample_counts))
rownames(sample_signatures)=rownames(sample_signatures_lower)=rownames(sample_signatures_upper)=colnames(sample_counts)
colnames(sample_signatures)=colnames(sample_signatures_lower)=colnames(sample_signatures_upper)=sigs_fitting
for(patient in patients){
  print(patient)
  tissues=unique(key_table$TissueType4[key_table$DonorID==patient])
  for(tissue in tissues){
    print(tissue)
    ref_sigs_present=sigs_per_tissue[[patient]][[tissue]]
    samples_tissue=key_table$SampleID[key_table$TissueType4==tissue&key_table$DonorID==patient]
    if(!is.null(ref_sigs_present)&!file.exists(paste0("pars_objects/",patient,"_",tissue,"_1540_2020_11_08 .Rdata"))){
      counts=t(sample_counts[,samples_tissue])
      colnames(counts)=colnames(cosmic_signatures_v3)
      #ref_sigs_present=names(sort(sig_order[ref_sigs_present]))
      ref_sigs_present=unique(c(ref_sigs_present,"SBS1","SBS5","SBS40"))
      
      fit=fit_signatures(counts=counts,signatures = final_sigs[ref_sigs_present,],
                         iter = 20000,warmup = 10000,model="poisson",chains = 2)
      pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
      rownames(pars$mean)=rownames(pars$lower)=rownames(pars$upper)=samples_tissue
      saveRDS(pars,paste0("pars_objects/",patient,"_",tissue,"_1540_2020_10_25.Rdata"))
      sample_signatures[rownames(pars$mean),colnames(pars$mean)]=as.matrix(pars$mean)
      sample_signatures_lower[rownames(pars$lower),colnames(pars$lower)]=as.matrix(pars$lower)
      sample_signatures_upper[rownames(pars$upper),colnames(pars$upper)]=as.matrix(pars$upper)
    }
  }
}
write.table(sample_signatures,"sample_signatures_SBS1540_2020_11_09.txt")
write.table(sample_signatures_lower,"sample_signatures_lower_SBS1540_2020_11_09.txt")
write.table(sample_signatures_upper,"sample_signatures_upper_SBS1540_2020_11_09.txt")

