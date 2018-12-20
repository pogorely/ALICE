#import functions into ALICE format.

#import adaptive:

aavec<-structure(c("(TTT|TTC)", "(TTA|TTG|CTT|CTC|CTA|CTG)", "(TCT|TCC|TCA|TCG|AGT|AGC)", 
                   "(TAT|TAC)", "(TAA|TAG|TGA)", "(TGT|TGC)", "(TGG)", "(CCT|CCC|CCA|CCG)", 
                   "(CAT|CAC)", "(CAA|CAG)", "(CGT|CGC|CGA|CGG|AGA|AGG)", "(ATT|ATC|ATA)", 
                   "(ATG)", "(ACT|ACC|ACA|ACG)", "(AAT|AAC)", "(AAA|AAG)", "(GTT|GTC|GTA|GTG)", 
                   "(GCT|GCC|GCA|GCG)", "(GAT|GAC)", "(GAA|GAG)", "(GGT|GGC|GGA|GGG)"
), .Names = c("F", "L", "S", "Y", "*", "C", "W", "P", "H", "Q", 
              "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))

Vlist<-c("TCRBV10-01","TCRBV10-02","TCRBV10-03","TCRBV11-01","TCRBV11-02","TCRBV11-03","TCRBV12-03","TCRBV12-04","TCRBV12-05","TCRBV13-01","TCRBV14-01","TCRBV15-01","TCRBV16-01","TCRBV18-01","TCRBV19-01","TCRBV02-01","TCRBV20-01","TCRBV21-01","TCRBV23-01","TCRBV24-01","TCRBV25-01","TCRBV27-01","TCRBV28-01","TCRBV29-01","TCRBV03-01","TCRBV30-01","TCRBV04-01","TCRBV04-02","TCRBV04-03","TCRBV05-01","TCRBV05-04","TCRBV05-05","TCRBV05-06","TCRBV05-08","TCRBV06-01","TCRBV06-02","TCRBV06-03","TCRBV06-04","TCRBV06-05","TCRBV06-06","TCRBV06-07","TCRBV07-01","TCRBV07-02","TCRBV07-03","TCRBV07-04","TCRBV07-06","TCRBV07-07","TCRBV07-08","TCRBV07-09","TCRBV09-01")
Jlist<-c("TCRBJ01-01","TCRBJ01-02","TCRBJ01-03","TCRBJ01-04","TCRBJ01-05","TCRBJ01-06","TCRBJ02-01","TCRBJ02-02","TCRBJ02-03","TCRBJ02-04","TCRBJ02-05","TCRBJ02-06","TCRBJ02-07")

#V-gene conversion tables

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


trim_aa<-function(x,trim=0){
  substr(x,1,nchar(x)-trim)
}

identify_cdr3aa_in_nuc<-function(aa,nuc){
  reg<-paste0(aavec[unlist(strsplit(aa,""))],collapse="")
  unlist(regmatches(nuc,gregexpr(reg,nuc)))
}

identify_cdr3aa_in_nuc_vector<-function(aa,nuc){
  resn<-character(length = length(nuc))
  for (i in 1:length(nuc))
  {
    resn[i]<-identify_cdr3aa_in_nuc(aa[i],nuc[i])   
  }
  resn
}

add_cdr3nt_to_immunoseq<-function(dt,trim=0){
  dt$CDR3nt<-identify_cdr3aa_in_nuc_vector(trim_aa(dt$CDR3.amino.acid.sequence,trim=trim),dt$CDR3.nucleotide.sequence)
  dt
}

cluster_immunoseq_by_CDR3<-function(x){
  x[,.(Read.count=sum(Read.count),Read.proportion=sum(Read.proportion),CDR3.amino.acid.sequence=unique(CDR3.amino.acid.sequence)),by=.(CDR3nt,bestVGene,bestJGene)]
}

read_immunoseq_folder<-function(folder,Read_thres=0){
  DTlist<-lapply(list.files(folder,full.names = T),fread)
  names(DTlist)<-list.files(folder,full.names = F)
  #print(names(DTlist))
  lapply(DTlist,setkey)
  DTlist<-lapply(DTlist,unique)
  #DTlist<-lapply(DTlist,function(x){if ("count"%in%colnames(x))setnames(x,"count","count (templates/reads)")})
  DTlist<-lapply(DTlist,function(x)x[aminoAcid!=""&`count (templates/reads)`>Read_thres,,])
  DTlist
}

read_immunoseq_folder_old<-function(folder,Read_thres=0){#parser for immunoseq legacy format, as in Emerson et al Plos ONE 2014
  DTlist<-lapply(list.files(folder,full.names = T),fread)
  names(DTlist)<-list.files(folder,full.names = F)
  lapply(DTlist,setnames,c("V1","V2","V3","V6","V8","V22"),c("nucleotide","aminoAcid","count","vMaxResolved","vGeneName","jGeneName"))
  #print(names(DTlist))
  lapply(DTlist,setkey)
  DTlist<-lapply(DTlist,unique)
  DTlist<-lapply(DTlist,function(x)x[aminoAcid!=""&`count`>Read_thres,,])
  DTlist
}

resolve_immunoseq_V<-function(dt){ #provides alleles for few V-families
  dt[vMaxResolved=="TCRBV20",vGeneName:="TCRBV20-01",]
  dt[vMaxResolved=="TCRBV12",vGeneName:="TCRBV12-03",]
  dt[vMaxResolved=="TCRBV06",vGeneName:="TCRBV06-02",]
  dt
}

convert_immunoseq<-function(dt){
  dt<-resolve_immunoseq_V(dt)
  dt[,Read.count:=`count (templates/reads)`,]
  dt[,Read.proportion:=prop.table(`count (templates/reads)`),]
  dt[,CDR3.nucleotide.sequence:=`nucleotide`,]
  dt[,CDR3.amino.acid.sequence:=`aminoAcid`,]
  dt[,bestVGene:=Vv[`vGeneName`],] #V gene conversion
  dt[,bestJGene:=Jv[`jGeneName`],] #J gene conversion
  dt[,.(Read.count, Read.proportion,CDR3.amino.acid.sequence, bestVGene, bestJGene,CDR3.nucleotide.sequence),]
}

convert_immunoseq_old<-function(dt){
  dt<-resolve_immunoseq_V(dt)
  dt[,Read.count:=`count`,]
  dt[,Read.proportion:=prop.table(`count`),]
  dt[,CDR3.nucleotide.sequence:=`nucleotide`,]
  dt[,CDR3.amino.acid.sequence:=`aminoAcid`,]
  dt[,bestVGene:=Vv[`vGeneName`],] #conversion!!!
  dt[,bestJGene:=Jv[`jGeneName`],] #conversion!!!
  dt[,.(Read.count,Read.proportion, CDR3.amino.acid.sequence, bestVGene, bestJGene,CDR3.nucleotide.sequence),]
}

import_immunoseq_pipeline<-function(folder,trim=0,Read_thres=1){
  DTlist<-read_immunoseq_folder(folder,Read_thres)
  DTlist<-lapply(DTlist,convert_immunoseq)
  DTlist<-lapply(DTlist,add_cdr3nt_to_immunoseq,trim=trim)
  DTlist<-lapply(DTlist,cluster_immunoseq_by_CDR3)
  DTlist
}

import_immunoseq_pipeline_old<-function(folder,trim=0,Read_thres=1){#special parser for legacy immunoseq format
  DTlist<-read_immunoseq_folder_old(folder,Read_thres)
  DTlist<-lapply(DTlist,convert_immunoseq_old)
  DTlist<-lapply(DTlist,add_cdr3nt_to_immunoseq,trim=trim)
  DTlist<-lapply(DTlist,cluster_immunoseq_by_CDR3)
  DTlist
}
#import vdjtools: 

