library(data.table)
library(igraph)
library(stringdist)
load("VDJT.rda")
source("generation.R")
#load generation

add_space<-function(df,hugedf,volume=66e6){#add space occupied by sequence neighbours. hugedf contains information of CDR3aa sequence and its generative probability. 
  huge<-hugedf$sim_num
  names(huge)<-hugedf$CDR3.amino.acid.sequence
  if(nrow(df)==0){return(0)}
  else{
    space<-numeric(nrow(df));
    tstl<-lapply(df$CDR3.amino.acid.sequence,all_other_variants_one_mismatch)
    for (i in 1:nrow(df))
    {
      space[i]=sum(huge[tstl[[i]]])
    }
    space_n=space/volume
    df$space=space
    df$space_n=space_n
    df
  }
}

add_p_val<-function(df,total,correct=9.41){#adds column with p_value  
  if(!is.null(nrow(df))){
    if(correct=="auto"){
      tmp=df$space_n*total
      correct<-coef(lm(df$D~tmp+0))[1]
    }  
    df$p_val<-pbinom(q=df$D,size = total,prob = correct*df$space_n,lower.tail = F)
    df}
  else{return(df)}
}

igraph_from_seqs<-function(seqs,max_errs=1) {
  graph<-graph.empty(n = length(seqs), directed=F)
  tmp<-stringdistmatrix(seqs,seqs,method="hamming")
  graph<-add.edges(graph, t(which(tmp<=max_errs,arr.ind=T)))
  graph<-igraph::simplify(graph)
  graph<-set.vertex.attribute(graph, 'label', V(graph), seqs)
  graph
}

filter_data_dt<-function(DT){
  DT[,leftgr:=.GRP,.(substr(CDR3.amino.acid.sequence,1,floor(nchar(CDR3.amino.acid.sequence)/2)),bestVGene,bestJGene)]
  DT[,rightgr:=.GRP,.(substr(CDR3.amino.acid.sequence,floor(nchar(CDR3.amino.acid.sequence)/2)+1,nchar(CDR3.amino.acid.sequence)),bestVGene,bestJGene)]
  DT[,D_left:=filter_data_no_igraph(.SD),.(leftgr)]
  DT[,D_right:=filter_data_no_igraph(.SD),.(rightgr)]
  DT[,D_id:=.N,.(CDR3.amino.acid.sequence)]
  DT[,D:=(D_left+D_right-D_id-1),]
  DT
}

filter_data<-function(df)
{
  gr<-igraph_from_seqs(df$CDR3.amino.acid.sequence)
  df$D=degree(gr)
  df$cl=clusters(gr)$membership
  df$total_n=nrow(df)
  #df
  df[df$D>0,]
}  

filter_data_no_igraph<-function(df)
{
  if (nrow(df)>1){
  tmp <- stringdistmatrix(df$CDR3.amino.acid.sequence,df$CDR3.amino.acid.sequence,method="hamming")
    rowSums(tmp<=1,na.rm = T)}
  else{1}
}  

select_sign<-function(sign_list,D_thres=2,P_thres=0.001,cor_method="BH"){ #performs multiple testing correction and returns list of significant results.
  lapply(sign_list,function(x)x[D>D_thres&space!=0,,][p.adjust(p_val,method=cor_method)<P_thres,,]) 
}

all_other_letters<-function(str,ind=8){
  aa<-c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", 
        "P", "Q", "R", "S", "T", "V", "W", "Y")
  paste0(substr(str,1,ind-1),aa,substr(str,ind+1,nchar(str)))
}

all_other_variants_one_mismatch<-function(str){
  unique(as.vector(sapply(2:(nchar(str)-1),all_other_letters,str=str)))
}

all_other_variants_one_mismatch_regexp<-function(str){
  unique(as.vector(sapply(2:(nchar(str)-1),function(x){tmp<-str;substr(tmp,x,x)<-"X";tmp})))
}

convert_comblist_to_df<-function(comblist)
{
  newl<-list()
  for (i in 1:length(comblist[[1]])){
    newl[[i]]<-lapply(comblist,"[[",i)
    names(newl)[i]<-names(comblist[[1]])[i]
    newl[[i]]<-do.call(rbind,newl[[i]][sapply(newl[[i]],is.list)])
  }
  newl
}

#pipeline functions
make_rda_folder<-function(DTlist,folder="",prefix="",VJDT=VDJT){
  dir.create(folder, showWarnings = FALSE)
  VJDT<-as.data.table(VJDT)
  VJDT[,bestVGene:=V,]
  VJDT[,bestJGene:=J,]
  for (i in 1:nrow(VJDT)){
    all_short_i<-lapply(DTlist,function(x)x[bestVGene==VJDT$bestVGene[i]&bestJGene==VJDT$bestJGene[i]&Read.count>1,,]) 
    all_short_int<-lapply(all_short_i,filter_data)
    all_short_int2<-lapply(all_short_int,function(x)x[D>2,])
    hugel<-unlist(lapply(unique(unlist(lapply(all_short_int2,function(x){if(nrow(x)>0)x[,CDR3.amino.acid.sequence,]}))),all_other_variants_one_mismatch))
    shrep<-data.frame(CDR3.amino.acid.sequence=unique(hugel))
    fname=paste0(prefix,VJDT[i,as.character(bestVGene),],"_",VJDT[i,as.character(bestJGene),],".rda",collapse="")
    if (nrow(shrep)!=0)
      save(shrep,file=paste0(folder,"/",fname,collapse = ""))
  }
}

compute_pgen_rda_folder<-function(folder,prefix="",iter=50,cores=8,nrec=5e5,silent=T){
  fnames<-list.files(folder,,full.names = T)
  fnames_s<-list.files(folder,full.names = F)
  fnames<-grep(pattern = "res_",fnames,invert = T,value = T)
  fnames_s<-grep(pattern = "res_",fnames_s,invert = T,value = T)
  fnames_s<-gsub(".rda","",fnames_s)
  fnames_s<-gsub(prefix,"",fnames_s)
  VJlist<-do.call(rbind,strsplit(fnames_s,"_"))
  for (i in 1:nrow(VJlist))
    if (VJlist[i,1]%in%segments$TRBV$V.alleles&VJlist[i,2]%in%segments$TRBJ$J.alleles)#test if present
    {
      if(!silent)print(fnames_s[i]) 
      if(!silent)print(format(Sys.time(), "%a %b %d %X %Y"))
      load(fnames[i])
      res<-data.frame()
      if(nrow(shrep)!=0)
        res<-estimate_pgen_aa(data = shrep,iter = iter,cores=cores,nrec=nrec,V=VJlist[i,1],J=VJlist[i,2])
      if(!silent)print("all iterations done")
      save(res,file=paste0(folder,"/","res_",prefix,fnames_s[i],".rda",collapse = ""))
      if(!silent)print(format(Sys.time(), "%a %b %d %X %Y"))
      if(!silent)print("result saved")
      rm(res)
    }
}

parse_rda_folder<-function(DTlist,folder,prefix="",Q=9.41,volume=66e6,silent=T){# gets folder, returns space and space_n, and add significance also.  
  fnames<-list.files(folder,pattern = "res_",full.names = T)
  fnames_s<-list.files(folder,pattern = "res_",full.names = F)
  fnames_s<-gsub("res_","",fnames_s)
  fnames_s<-gsub(".rda","",fnames_s)
  fnames_s<-gsub(prefix,"",fnames_s)
  VJlist<-do.call(rbind,strsplit(fnames_s,"_"))
  resl<-list()
  for (i in 1:nrow(VJlist)){
    if(!silent)print(i)
    all_short_i<-lapply(DTlist,function(x)x[bestVGene==VJlist[i,1]&bestJGene==VJlist[i,2]&Read.count>1,,]) 
    all_short_int<-lapply(all_short_i,filter_data)
    all_short_int2<-lapply(all_short_int,function(x)x[D>2,])
    load(fnames[i])
    all_short_int2_space<-lapply(all_short_int2,add_space,hugedf = res,volume=volume)
    for (j in 1:length(all_short_int2_space))
    { 
      all_short_int2_space[[j]]<-add_p_val(all_short_int2_space[[j]],total = nrow(all_short_i[[j]][Read.count>1,,]),correct=Q)
    }
    resl[[i]]<-all_short_int2_space
  }
  names(resl)<-fnames_s
  resl
}

olga_parallel_wrapper_beta<-function(DT,cores=1){
  DT<-DT[!grepl(CDR3.amino.acid.sequence,pattern = "*",fixed = T)&((nchar(CDR3.nucleotide.sequence)%%3)==0)]
  DT<-DT[bestVGene%in%row.names(OLGAVJ)&bestJGene%in%colnames(OLGAVJ)]
#  DT[,ind:=1:.N,] this is wrong for tmp2
  fn<-paste0("tmp",1:cores,".tsv")
  fn2<-paste0("tmp_out",1:cores,".tsv")
  
  for (f in c(fn,fn2))if (file.exists(f)) file.remove(f)
  
  DTl<-split(DT, sort((1:nrow(DT)-1)%%cores+1))
  for (i in 1:length(DTl))
  write.table(as.data.frame(DTl[[i]][,.(CDR3.amino.acid.sequence,bestVGene,bestJGene,ind),]),quote=F,row.names = F,sep = "\t",file = fn[i])
  olga_commands<-paste0("olga-compute_pgen --humanTRB --display_off --seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i tmp",1:cores,".tsv -o tmp_out",1:cores,".tsv")
  system(paste0(olga_commands,collapse=" & "))
  readline(prompt="Press [enter] to continue")
  fnt<-do.call(rbind,lapply(fn2,fread))
  DT$Pgen<-fnt$V2
  DT
}
olga_parallel_wrapper_alpha<-function(DT,cores=1,withoutVJ=F){#test_it
  fn<-paste0("tmp",cores,".tsv")
  for (f in fn)if (file.exists(f)) file.remove(f)
  DTl<-split(DT, sort((1:nrow(DT)-1)%%cores+1))
  for (i in 1:length(DTl))
    write.table(as.data.frame(DTl[[i]][,.(CDR3.amino.acid.sequence,bestVGene,bestJGene,ind),]),quote=F,row.names = F,sep = "\t",file = fn[i])
  
  olga_commands<-paste0("olga-compute_pgen --humanTRA --display_off --seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i tmp",1:cores,".tsv -o tmp_out",1:cores,".tsv")
  if (withoutVJ)   olga_commands<-paste0("olga-compute_pgen --humanTRA --display_off --seq_in 0  --lines_to_skip 1 -d 'tab' -i tmp",1:cores,".tsv -o tmp_out",1:cores,".tsv")

  system(paste0(olga_commands,collapse=" & "))
  fn<-paste0("tmp_out",cores,".tsv")
  fnt<-do.call(rbind,lapply(fn,fread))
  DT$Pgen<-fnt$V2
  DT
}

output_olga_DT<-function(DT,path="",Q=9.41)
{
  load("OLGA_V_J_hum_beta.rda")
  DT<-DT[!grepl(CDR3.amino.acid.sequence,pattern = "*",fixed = T)&((nchar(CDR3.nucleotide.sequence)%%3)==0)]
  DT<-DT[bestVGene%in%row.names(OLGAVJ)&bestJGene%in%colnames(OLGAVJ)]
  
  tmp<-filter_data_dt(DT)
  tmp[,n_total:=.N,.(bestVGene,bestJGene)]
  tmp<-tmp[D>2][,ind:=1:.N,] #sequences themselves
  tmp2<-tmp[,.(bestVGene,bestJGene,CDR3.amino.acid.sequence=all_other_variants_one_mismatch_regexp(CDR3.amino.acid.sequence)),ind] #all one mismatch variants with X (regexp!)
  write.table(as.data.frame(tmp[,.(CDR3.amino.acid.sequence,bestVGene,bestJGene,ind),]),quote=F,row.names = F,sep = "\t",file = "tmp.tsv")
  write.table(as.data.frame(tmp2[,.(CDR3.amino.acid.sequence,bestVGene,bestJGene,ind),]),quote=F,row.names = F,sep = "\t",file = "tmp2.tsv")
  fn <- "tmp_out.tsv"
  if (file.exists(fn)) file.remove(fn)
  fn <- "tmp2_out.tsv"
  if (file.exists(fn)) file.remove(fn)
  system("olga-compute_pgen --humanTRB --display_off --seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i tmp.tsv -o tmp_out.tsv")
  system("olga-compute_pgen --humanTRB --display_off --seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i tmp2.tsv -o tmp2_out.tsv")
  probstmp<-fread("tmp_out.tsv")
  tmp$Pgen<-probstmp$V2
  probstmp2<-fread("tmp2_out.tsv")
  tmp2$Pgen<-probstmp2$V2
  tmp$Pgen1<-tmp2[,sum(Pgen),ind]$V1 #add all nums... add p-values...
  tmp[,Pgen3:=Pgen1-Pgen*(nchar(CDR3.amino.acid.sequence)-2),]
  tmp[,p_value:=ppois(D,lambda = 3*Q*n_total*Pgen3/(OLGAVJ[cbind(bestVGene,bestJGene)]),lower.tail = F),]# 007 is VJ_comb coeff!!!
  tmp
}
output_olga_DT_alpha<-function(DT,path="",Q=9.41)
{
  load("OLGA_A_hum_alpha.rda")
  DT<-DT[!grepl(CDR3.amino.acid.sequence,pattern = "*",fixed = T)&((nchar(CDR3.nucleotide.sequence)%%3)==0)]
  DT<-DT[bestVGene%in%row.names(OLGAVJA)&bestJGene%in%colnames(OLGAVJA)]
  tmp<-filter_data_dt(DT)
  tmp[,n_total:=.N,.(bestVGene,bestJGene)]
  tmp<-tmp[D>2][,ind:=1:.N,] #sequences themselves
  tmp2<-tmp[,.(bestVGene,bestJGene,CDR3.amino.acid.sequence=all_other_variants_one_mismatch_regexp(CDR3.amino.acid.sequence)),ind] #all one mismatch variants with X (regexp!)
  write.table(as.data.frame(tmp[,.(CDR3.amino.acid.sequence,bestVGene,bestJGene,ind),]),quote=F,row.names = F,sep = "\t",file = "tmp.tsv")
  write.table(as.data.frame(tmp2[,.(CDR3.amino.acid.sequence,bestVGene,bestJGene,ind),]),quote=F,row.names = F,sep = "\t",file = "tmp2.tsv")
  fn <- "tmp_out.tsv"
  if (file.exists(fn)) file.remove(fn)
  fn <- "tmp2_out.tsv"
  if (file.exists(fn)) file.remove(fn)
  system("olga-compute_pgen --humanTRA --display_off --seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i tmp.tsv -o tmp_out.tsv")
  system("olga-compute_pgen --humanTRA --display_off --seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i tmp2.tsv -o tmp2_out.tsv")
  probstmp<-fread("tmp_out.tsv")
  tmp$Pgen<-probstmp$V2
  probstmp2<-fread("tmp2_out.tsv")
  tmp2$Pgen<-probstmp2$V2
  tmp$Pgen1<-tmp2[,sum(Pgen),ind]$V1 #add all nums... add p-values...
  tmp[,Pgen3:=Pgen1-Pgen*(nchar(CDR3.amino.acid.sequence)-2),]
  tmp[,p_value:=ppois(D,lambda = 3*Q*n_total*Pgen3/(OLGAVJA[cbind(bestVGene,bestJGene)]),lower.tail = F),]# 007 is VJ_comb coeff!!!
  tmp
}  

#run pipeline function.
ALICE_pipeline<-function(DTlist,folder="",cores=8,iter=50,nrec=5e5,P_thres=0.001,cor_method="BH")
{
  make_rda_folder(DTlist,folder) #generate .rda files for CDR3aa gen prob estimation for each VJ
  compute_pgen_rda_folder(folder,cores=cores,nrec=nrec,iter=iter) #estimate CDR3aa gen prob for each sequence and save to separate res_ files
  results<-parse_rda_folder(DTlist,folder,volume = cores*iter*nrec/3) #parse res_ files
  results<-convert_comblist_to_df(results) #convert to single dataset from VJ-combs
  select_sign(results,P_thres=P_thres,cor_method=cor_method) #filter for significant results
}
