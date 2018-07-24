#CDR3 AA generation probability estimation.

library(parallel)
library(Biostrings)

load("setP.rda")

#immmunoseq -> IMGT V and J segments naming convention. 
Vv<-structure(c("TRBV10-1", "TRBV10-2", "TRBV10-3", "TRBV11-1", "TRBV11-2", 
                "TRBV11-3", "TRBV12-3", "TRBV12-4", "TRBV12-5", "TRBV13", "TRBV14", 
                "TRBV15", "TRBV16", "TRBV18", "TRBV19", "TRBV2", "TRBV20-1", 
                "TRBV21-1", "TRBV23-1", "TRBV24-1", "TRBV25-1", "TRBV27", "TRBV28", 
                "TRBV29-1", "TRBV3-1", "TRBV30", "TRBV4-1", "TRBV4-2", "TRBV4-3", 
                "TRBV5-1", "TRBV5-4", "TRBV5-5", "TRBV5-6", "TRBV5-8", "TRBV6-1", 
                "TRBV6-2", "TRBV6-3", "TRBV6-4", "TRBV6-5", "TRBV6-6", "TRBV6-7", 
                "TRBV7-1", "TRBV7-2", "TRBV7-3", "TRBV7-4", "TRBV7-6", "TRBV7-7", 
                "TRBV7-8", "TRBV7-9", "TRBV9"), .Names = c("TCRBV10-01", "TCRBV10-02", 
                                                           "TCRBV10-03", "TCRBV11-01", "TCRBV11-02", "TCRBV11-03", "TCRBV12-03", 
                                                           "TCRBV12-04", "TCRBV12-05", "TCRBV13-01", "TCRBV14-01", "TCRBV15-01", 
                                                           "TCRBV16-01", "TCRBV18-01", "TCRBV19-01", "TCRBV02-01", "TCRBV20-01", 
                                                           "TCRBV21-01", "TCRBV23-01", "TCRBV24-01", "TCRBV25-01", "TCRBV27-01", 
                                                           "TCRBV28-01", "TCRBV29-01", "TCRBV03-01", "TCRBV30-01", "TCRBV04-01", 
                                                           "TCRBV04-02", "TCRBV04-03", "TCRBV05-01", "TCRBV05-04", "TCRBV05-05", 
                                                           "TCRBV05-06", "TCRBV05-08", "TCRBV06-01", "TCRBV06-02", "TCRBV06-03", 
                                                           "TCRBV06-04", "TCRBV06-05", "TCRBV06-06", "TCRBV06-07", "TCRBV07-01", 
                                                           "TCRBV07-02", "TCRBV07-03", "TCRBV07-04", "TCRBV07-06", "TCRBV07-07", 
                                                           "TCRBV07-08", "TCRBV07-09", "TCRBV09-01"))
Jv<-structure(c("TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", 
            "TRBJ1-6", "TRBJ2-1", "TRBJ2-2", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", 
            "TRBJ2-6", "TRBJ2-7"), .Names = c("TCRBJ01-01", "TCRBJ01-02", 
                                              "TCRBJ01-03", "TCRBJ01-04", "TCRBJ01-05", "TCRBJ01-06", "TCRBJ02-01", 
                                              "TCRBJ02-02", "TCRBJ02-03", "TCRBJ02-04", "TCRBJ02-05", "TCRBJ02-06", 
                                              "TCRBJ02-07"))





estimate_pgen_aa<-function(data,V="TRBV9",J="TRBJ2-3",iter=3,nrec=5e5,cores=1,colname="CDR3.amino.acid.sequence",pr=beta.prob,segments=segmentsc2){
  res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    for (i in 1:iter)
                    {r<-r+table(factor(gen_beta(n = nrec,V,J,pr=pr,segments = segments),levels = data[[colname]]))}
                    r
                    },mc.cores = cores)  
  sim_num<-rowSums(do.call(cbind,res)) 
  data<-cbind(sim_num,data)
  data
}

estimate_pgen_aa_alpha<-function(data,V="TRBV9",J="TRBJ2-3",iter=3,nrec=5e5,cores=1,colname="CDR3.amino.acid.sequence"){
  res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    for (i in 1:iter)
                    {r<-r+table(factor(gen_alpha(n = nrec,V,J),levels = data[[colname]]))}
                    r
                  },mc.cores = cores)  
  sim_num<-rowSums(do.call(cbind,res)) 
  data<-cbind(sim_num,data)
  data
}

estimate_pgen_aa_alpha_repertoire<-function(data,VJ_usage=VJ_usage_alpha_tw,iter=3,nrec=5e5,cores=1,colname="CDR3.amino.acid.sequence"){
  res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    for (i in 1:iter)
                    {r<-r+table(factor(gen_alpha_repertoire_CDR3(n = nrec,VJ_usage = VJ_usage),levels = data[[colname]]))}
                    r
                  },mc.cores = cores)  
  sim_num<-rowSums(do.call(cbind,res)) 
  data<-cbind(sim_num,data)
  data
}

estimate_pgen_aa_beta_repertoire<-function(data,VJ_usage=get_vj("../Pattern/VDJDB/public-epitope-master/rearr_model/hip_nonfunc_vj_usage.txt"),iter=3,nrec=5e5,cores=1,colname="CDR3.amino.acid.sequence"){
  res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    for (i in 1:iter)
                    {r<-r+table(factor(gen_beta_repertoire_CDR3(n = nrec,VJ_usage = VJ_usage),levels = data[[colname]]))}
                    r
                  },mc.cores = cores)  
  sim_num<-rowSums(do.call(cbind,res)) 
  data<-cbind(sim_num,data)
  data
}

estimate_pgen_aa_VJ_alpha_repertoire<-function(data,VJ_usage=VJ_usage_alpha_tw,iter=3,nrec=5e5,cores=1,colname="CDR3VJ"){
  res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    for (i in 1:iter)
                    {r<-r+table(factor(gen_alpha_repertoire_short_VJ(n = nrec,VJ_usage = VJ_usage),levels = data[[colname]]))}
                    r
                  },mc.cores = cores)  
  sim_num<-rowSums(do.call(cbind,res)) 
  data<-cbind(sim_num,data)
  data
}

estimate_pgen_aa_VJ_beta_repertoire<-function(data,VJ_usage=get_vj("../Pattern/VDJDB/public-epitope-master/rearr_model/hip_nonfunc_vj_usage.txt"),iter=3,nrec=5e5,cores=1,colname="CDR3VJ"){
  res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    for (i in 1:iter)
                    {r<-r+table(factor(gen_beta_repertoire_short_VJ(n = nrec,VJ_usage = VJ_usage),levels = data[[colname]]))}
                    r
                  },mc.cores = cores)  
  sim_num<-rowSums(do.call(cbind,res)) 
  data<-cbind(sim_num,data)
  data
}

estimate_pgen_aa_div<-function(data,V="TRBV9",J="TRBJ2-3",iter=3,nrec=5e5,cores=1,colname="CDR3.amino.acid.sequence",targets=""){
    res<-mclapply(1:cores,FUN = 
                  function(x){
                    r<-rep(0,nrow(data));
                    vars<-character()
                    for (i in 1:iter)
                      {
                      gen<-gen_beta(n = nrec,V,J,translate=F)
                      genaa<-translate_cdr(gen)
                      r<-r+table(factor(gen_beta(n = nrec,V,J),levels = data[[colname]]))
                      vars<-c(vars,gen[genaa%in%targets])
                      }
                      
                      list(r,vars)
                    },mc.cores = cores) 
  sim_num<-rowSums(do.call(cbind,lapply(res,"[[",1))) 
  data<-cbind(sim_num,data)
  vars=unlist(lapply(res,"[[",2))
  list(data=data,variants=vars)
}

translate_cdr<-function(.seq){
  as.character(translate(DNAStringSet(.seq)))}

markov_chain<-function(x,probsm){
  res<-rep("",length(x))
  res[x=="A"]<-sample(size = length(x[x=="A"]),x = c("A","C","G","T"),replace = T,prob = probsm["A",1:4])
  res[x=="C"]<-sample(size = length(x[x=="C"]),x = c("A","C","G","T"),replace = T,prob = probsm["C",1:4])
  res[x=="G"]<-sample(size = length(x[x=="G"]),x = c("A","C","G","T"),replace = T,prob = probsm["G",1:4])
  res[x=="T"]<-sample(size = length(x[x=="T"]),x = c("A","C","G","T"),replace = T,prob = probsm["T",1:4])
  res
}

mchain1<-function(xstart,ls,probsm){
  res<-xstart
  for (i in 1:max(ls)) 
  {res[ls>=i]<-paste0(res[ls>=i],markov_chain(substr(res[ls>=i],nchar(res[ls>=i]),nchar(res[ls>=i])),probsm))}
  substr(res,2,nchar(res))
}

mchain2<-function(xstart,ls,probsm){
  res<-xstart
  for (i in 1:max(ls)) 
  {res[ls>=i]<-paste0(markov_chain(substr(res[ls>=i],nchar(res[ls>=i]),nchar(res[ls>=i])),probsm),res[ls>=i])}
  substr(res,1,nchar(res)-1)
}

#universal low_level VDJ recombination function. 
gen_VDJ<-function(n,Vseq,Jseq,Dseq,delV,delJ,delD,insLenVD,insLenDJ,insNuclVD,insNuclDJ,sep="",inframe=T,translate=T){
  if (n==0) return()
  #trim Vs according to delV vec
  Vpart<-substr(rep(Vseq,times = n),1,nchar(Vseq)-sample(0:(length(delV)-1),size = n,replace = T,prob = delV))
  #trim Js according to delJ vec
  Jpart<-substr(rep(Jseq,times = n),1+sample(0:(length(delJ)-1),size = n,replace = T,prob = delJ),nchar(Jseq))
  #generate VD insertion lengths according insLenVD
  insVD.num.vec<-sample(0:(length(insLenVD)-1),size = n,replace = T, prob=insLenVD)
  #generate VD insertion sequences with markov chain according to insNuclVD
  VDins<-mchain1(xstart = substr(Vpart,nchar(Vpart),nchar(Vpart)),ls = insVD.num.vec,probsm = insNuclVD)
  #generate DJ insertion lengths according insLenVD
  insDJ.num.vec<-sample(0:(length(insLenVD)-1),size = n,replace = T, prob=insLenDJ)
  #generate DJ insertion sequences with markov chain according to insNuclVD
  DJins<-mchain2(xstart = substr(Jpart,1,1),ls = insDJ.num.vec,probsm = insNuclDJ)
  #preprocess D deletion matrix for efficient subsetting.
  colnames(delD)<-1:ncol(delD)
  row.names(delD)<-1:nrow(delD)
  indmat<-as.data.frame(as.table(delD))
  indmat<-cbind(as.numeric(indmat[,1]),as.numeric(indmat[,2]))
  #generate table for D substings
  Dmat<-t(matrix(indmat[sample(1:length(delD),n,prob = delD,replace = T),c(1,2)],ncol=2))
  Dtrimmed<-substr(rep(Dseq,times = n),Dmat[1,],nchar(Dseq)-Dmat[2,]+1)
  #filter result : we now have problems here!
  res<-(paste(Vpart,VDins,Dtrimmed,DJins,Jpart,sep = sep))
  res
}

getins<-function(n)paste0(sample(c("A","T","G","C"),replace = T,size=n),collapse="")

#universal low_level VJ recombination function. 
gen_VJ<-function(n,Vseq,Jseq,delV,delJ,insLen,sep="",inframe=T,translate=T){
  if (n==0) return()
  #trim Vs according to delV vec
  Vpart<-substr(rep(Vseq,times = n),1,nchar(Vseq)-sample(0:(length(delV)-1),size = n,replace = T,prob = delV))
  #trim Js according to delJ vec
  Jpart<-substr(rep(Jseq,times = n),1+sample(0:(length(delJ)-1),size = n,replace = T,prob = delJ),nchar(Jseq))
  #generate VD insertion lengths according insLenVD
  ins.num.vec<-sample(0:(length(insLen)-1),size = n,replace = T, prob=insLen)
  VJins<-sapply(ins.num.vec,getins)
  #generate VD insertion sequences with markov chain according to insNuclVD
  #VDins<-mchain1(xstart = substr(Vpart,nchar(Vpart),nchar(Vpart)),ls = insVD.num.vec,probsm = insNuclVD)
  #generate DJ insertion lengths according insLenVD
  #insDJ.num.vec<-sample(0:(length(insLenVD)-1),size = n,replace = T, prob=insLenDJ)
  #generate DJ insertion sequences with markov chain according to insNuclVD
  #DJins<-mchain2(xstart = substr(Jpart,1,1),ls = insDJ.num.vec,probsm = insNuclDJ)
  #preprocess D deletion matrix for efficient subsetting.
  #colnames(delD)<-1:ncol(delD)
  #row.names(delD)<-1:nrow(delD)
  #indmat<-as.data.frame(as.table(delD))
  #indmat<-cbind(as.numeric(indmat[,1]),as.numeric(indmat[,2]))
  #generate table for D substings
  #Dmat<-t(matrix(indmat[sample(1:length(delD),n,prob = delD,replace = T),c(1,2)],ncol=2))
  #Dtrimmed<-substr(rep(Dseq,times = n),Dmat[1,],nchar(Dseq)-Dmat[2,]+1)
  #filter result : we now have problems here!
  res<-(paste(Vpart,VJins,Jpart,sep = sep))
  res
}
segmentsc2<-segments
#function to generate beta chain sequence
gen_beta<-function(n,V,J,pr=beta.prob,translate=T,inframe_only=T,segments=segmentsc2,sep="")
{
  tmp<-pr$P.ins.nucl[,1]
  pr$P.ins.nucl<-as.matrix(pr$P.ins.nucl[,-1])
  pr$P.ins.len<-pr$P.ins.len[1:15,]
  row.names(pr$P.ins.nucl)<-tmp
  Ddist<-round(prop.table(pr$P.J.D[row.names(pr$P.J.D)==J,])*n) 
  Vi<-which(colnames(pr$P.del.V)==V)
  Ji<-which(colnames(pr$P.del.J)==J)
  if((V%in%segments$TRBV$V.alleles)&(J%in%segments$TRBJ$J.alleles))
    {
  Vs<-segments$TRBV$Nucleotide.sequence[segments$TRBV$V.alleles==V]
  Js<-segments$TRBJ$Nucleotide.sequence[segments$TRBJ$J.alleles==J]
  res<-c(gen_VDJ(n=Ddist[1],
          Vseq = Vs,
          Jseq = Js,
          Dseq = segments$TRBD[1,]$Nucleotide.sequence,
          delV = pr$P.del.V[,Vi],
          delJ = pr$P.del.J[,Ji],
          delD = pr$P.del.D1,
          insLenVD =pr$P.ins.len[,2], 
          insLenDJ = pr$P.ins.len[,3],
          insNuclVD = pr$P.ins.nucl[,1:4],
          insNuclDJ = pr$P.ins.nucl[,5:8],
          sep = sep ),
    gen_VDJ(n=Ddist[2],
            Vseq = Vs,
            Jseq = Js,
            Dseq = segments$TRBD[2,]$Nucleotide.sequence,
            delV = pr$P.del.V[,Vi],
            delJ = pr$P.del.J[,Ji],
            delD = pr$P.del.D2,
            insLenVD =pr$P.ins.len[,2], 
            insLenDJ = pr$P.ins.len[,3],
            insNuclVD = pr$P.ins.nucl[,1:4],
            insNuclDJ = pr$P.ins.nucl[,5:8],
            sep = sep ))
  if(inframe_only|translate)res<-res[nchar(res)%%3==0]
  if(translate){translate_cdr(res)}
  else{res}}
}

gen_alpha<-function(n,V,J,pr=alpha_twins_average,translate=T,inframe_only=T,segments=segmentsc)
{
  #tmp<-pr$P.ins.nucl[,1]
  #pr$P.ins.nucl<-as.matrix(pr$P.ins.nucl[,-1])
  pr$P.ins.len<-pr$P.ins.len[1:15,]
  #row.names(pr$P.ins.nucl)<-tmp
  #Ddist<-round(prop.table(pr$P.J.D[row.names(pr$P.J.D)==J,])*n) 
  Vi<-which(colnames(pr$P.del.V)==V)
  Ji<-which(colnames(pr$P.del.J)==J)
  if((V%in%segments$TRAV$V.alleles)&(J%in%segments$TRAJ$J.alleles))
  {
    Vs<-segments$TRAV$Nucleotide.sequence[segments$TRAV$V.alleles==V]
    Js<-segments$TRAJ$Nucleotide.sequence[segments$TRAJ$J.alleles==J]
    res<-gen_VJ(n= n,
                   Vseq = Vs,
                   Jseq = Js,
                   delV = pr$P.del.V[,Vi],
                   delJ = pr$P.del.J[,Ji],
                   insLen =pr$P.ins.len[,2], 
                   sep = "" )
    if(inframe_only|translate)res<-res[nchar(res)%%3==0]
    if(translate){translate_cdr(res)}
    else{res}}
}


get_vj<-function(file){
  VJ_aging<-fread(file)
  names(VJ_aging)[3]<-"clonotypes"
  VJ_aging[,p:=prop.table(clonotypes),]
  VJ_aging[,VJ:=paste0(v,"_",j),]
  setkey(VJ_aging,VJ)
  VJ_aging
}

gen_beta_repertoire<-function(n,VJ_usage,pr=beta.prob,translate=T,inframe_only=T,sep=""){
  #sample n into VJ_classes. 
  #VJ usage is table V J prob. 
  VJ_usage$sample<-rmultinom(1,size=n,prob = VJ_usage$clonotypes)
  resl<-list()
  for (i in 1:nrow(VJ_usage)){
    if(VJ_usage$sample[i]!=0){
    gen_seq<-gen_beta(n=VJ_usage$sample[i],V=VJ_usage$v[i],J=VJ_usage$j[i],translate=translate, inframe_only = inframe_only,sep=sep)
    if(length(gen_seq)!=0)
    resl[[i]]<-data.frame(seq=gen_seq,V=VJ_usage$v[i],J=VJ_usage$j[i],stringsAsFactors = F)}
  }
  do.call(rbind,resl)
}

gen_beta_repertoire_short<-function(n,VJ_usage,pr=beta.prob,translate=T,inframe_only=T){
  #sample n into VJ_classes. 
  #VJ usage is table V J prob. 
  VJ_usage$sample<-rmultinom(1,size=n,prob = VJ_usage$clonotypes)
  resl<-list()
  for (i in 1:nrow(VJ_usage)){
    if(VJ_usage$sample[i]!=0){
      gen_seq<-gen_beta(n=VJ_usage$sample[i],V=VJ_usage$v[i],J=VJ_usage$j[i])
      if(length(gen_seq)!=0)
        resl[[i]]<-data.frame(seq=gen_seq,V=VJ_usage$v[i],J=VJ_usage$j[i],stringsAsFactors = F)}
  }
  tmp<-do.call(rbind,resl)
  paste0(tmp$seq,"_",tmp$V)
}

gen_beta_repertoire_short_VJ<-function(n,VJ_usage,pr=beta.prob,translate=T,inframe_only=T){
  #sample n into VJ_classes. 
  #VJ usage is table V J prob. 
  VJ_usage$sample<-rmultinom(1,size=n,prob = VJ_usage$clonotypes)
  resl<-list()
  for (i in 1:nrow(VJ_usage)){
    if(VJ_usage$sample[i]!=0){
      gen_seq<-gen_beta(n=VJ_usage$sample[i],V=VJ_usage$v[i],J=VJ_usage$j[i])
      if(length(gen_seq)!=0)
        resl[[i]]<-data.frame(seq=gen_seq,V=VJ_usage$v[i],J=VJ_usage$j[i],stringsAsFactors = F)}
  }
  tmp<-do.call(rbind,resl)
  paste0(tmp$seq,"_",tmp$V,"_",tmp$J)
}

gen_alpha_repertoire_short_VJ<-function(n,VJ_usage,pr=beta.prob,translate=T,inframe_only=T){
  #sample n into VJ_classes. 
  #VJ usage is table V J prob. 
  VJ_usage$sample<-rmultinom(1,size=n,prob = VJ_usage$clonotypes)
  resl<-list()
  for (i in 1:nrow(VJ_usage)){
    if(VJ_usage$sample[i]!=0){
      gen_seq<-gen_alpha(n=VJ_usage$sample[i],V=VJ_usage$v[i],J=VJ_usage$j[i])
      if(length(gen_seq)!=0)
        resl[[i]]<-data.frame(seq=gen_seq,V=VJ_usage$v[i],J=VJ_usage$j[i],stringsAsFactors = F)}
  }
  tmp<-do.call(rbind,resl)
  paste0(tmp$seq,"_",tmp$V,"_",tmp$J)
}

gen_beta_repertoire_CDR3<-function(n,VJ_usage,pr=beta.prob,translate=T,inframe_only=T){
  #sample n into VJ_classes. 
  #VJ usage is table V J prob. 
  VJ_usage$sample<-rmultinom(1,size=n,prob = VJ_usage$clonotypes)
  resl<-list()
  for (i in 1:nrow(VJ_usage)){
    if(VJ_usage$sample[i]!=0){
      gen_seq<-gen_beta(n=VJ_usage$sample[i],V=VJ_usage$v[i],J=VJ_usage$j[i])
      if(length(gen_seq)!=0)
        resl[[i]]<-gen_seq}
  }
  do.call(c,resl)
}

gen_alpha_repertoire_CDR3<-function(n,VJ_usage,pr=alpha.prob,translate=T,inframe_only=T){
  #sample n into VJ_classes. 
  #VJ usage is table V J prob. 
  VJ_usage$sample<-rmultinom(1,size=n,prob = VJ_usage$clonotypes)
  resl<-list()
  for (i in 1:nrow(VJ_usage)){
    if(VJ_usage$sample[i]!=0){
      gen_seq<-gen_alpha(n=VJ_usage$sample[i],V=VJ_usage$v[i],J=VJ_usage$j[i])
      if(length(gen_seq)!=0)
        resl[[i]]<-gen_seq}
  }
  do.call(c,resl)
}
