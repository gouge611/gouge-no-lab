#Args <- commandArgs()

###########function()###################

# ir_whole <- function(intron_counts,intron_length_list, exon_counts, exon_length_list){
# 	intron_length = sum(as.numeric(unlist(strsplit(as.character(intron_length_list),','))))
# 	exon_length = sum(as.numeric(unlist(strsplit(as.character(exon_length_list),','))))
# 	iri = (as.numeric(intron_counts) / intron_length) / (as.numeric(exon_counts) / exon_length)
# 	return(iri)
# }

#ir_neigh <- function(){}

cor_sense_and_antisense <- function(sense_list, antisense_list,sense_symbol,antisense_symbol,chrom){
	sense=as.numeric(unlist(strsplit(as.character(sense_list),',')))
	antisense=as.numeric(unlist(strsplit(as.character(antisense_list),',')))
	d_matrix=cbind(sense,antisense)
	d_matrix[which(d_matrix==0)] <-NA
	df=as.data.frame(d_matrix)
	if (dim(na.omit(df))[1] >= 3){
		ht=cor.test(df$sense,df$antisense)
		output=cbind(chrom,sense_symbol,antisense_symbol,ht$estimate,ht$parameter+2,ht$p.value)
		return(output)
	}
	else{
	  return(cbind(rep(NA,6)))
	}
}

exec_cor <- function(dat){
  len=length(dat)
  #   cor_sIR_asExpr=apply(dat,1,function(x) {
  #   cor_sense_and_antisense(x[5],x[9],x[1],x[6],x[2])})
  cor_sIR_asExpr=apply(dat,1,function(x) {
    cor_sense_and_antisense(x[11],x[13],x[1],x[6],x[2])})
  dat=cbind(dat,t(cor_sIR_asExpr)[,4:6])
  a=len+1
  b=len+3
  names(dat)[a:b]=c('corr_sIR_asExpr','sample_num','p_value')
  
  cor_sExpr_asExpr=apply(dat,1,function(x) {
    cor_sense_and_antisense(x[12],x[13],x[1],x[6],x[2])})
  dat=cbind(dat,t(cor_sExpr_asExpr)[,4:6])
  a=len+4
  b=len+6
  names(dat)[a:b]=c('corr_sExpr_asExpr','sample_num','p_value')
  
  return(dat)
}

suru <- function(filename){
  #update in 30/3/2017: field index change
  dat=read.table(filename,header = F)
  labels = c('sense_gene_symbol', 'chr', 'strand',
             'sense_exon_count_normalization', 
             'sense_intron_count_normalization',
             'antisense_gene_symbol', 'chr', 'strand',
             'antisense_exon_count_normalization', 
             'antisense_intron_count_normalization',
             'sense_expr_corrected', 'sense_ir_corrected',
             'antisense_expr_corrected')
  names(dat)=labels
  
  dat=exec_cor(dat)
  dat=na.omit(dat)
  write.table(dat,sub("_for_R",".txt",filename),sep = '\t',quote = F,row.names = F)
  # return(dat)
}

figure <- function(filename){
  df=read.table(filename,header = T)
  png(sub(".txt",".png",filename))
  par(mfrow=c(1,2))
  hist(df$corr_sIR_asExpr)
  hist(df$corr_sExpr_asExpr)
  dev.off()
}

group_fisher_test <- function(dat, threshold){
  q1=subset(dat,corr_sIR_asExpr>=threshold & corr_sExpr_asExpr>=threshold)
  q2=subset(dat,corr_sIR_asExpr<=-threshold & corr_sExpr_asExpr>=threshold)
  q3=subset(dat,corr_sIR_asExpr<=-threshold & corr_sExpr_asExpr<=-threshold)
  q4=subset(dat,corr_sIR_asExpr>=threshold & corr_sExpr_asExpr<=-threshold)
  
  m=matrix(c(dim(q1)[1], dim(q2)[1], dim(q3)[1], dim(q4)[1]), ncol = 2)
  return(fisher.test(m))
}

# r=apply(dat,1,function(x) {cor_sense_and_antisense(x[4],x[8],x[1],x[5],x[2])})
# write.table(t(r),"23333",quote=F,sep='\t',row.names=F,col.names=F)

############main()################

#dat = read.table(Args[6],sep="\t",header=FALSE)
#labels = c('gene_symbol', 'chr', 'strand', 
#            'shared_exon_start', 'shared_exon_end', 
#            'shared_intron_start', 'shared_intron_end', 
#            'shared_exon_length_list', 'shared_intron_lengh_list', 
#            'exon_count_total', 'exon_count', 'exon_count_normalization', 
#            'intron_count_total', 'intron_count', 'intron_count_normalization')
#names(dat)=labels
#dat$iri_whole=ir_whole(dat$intron_count_total,dat$shared_intron_lengh_list,dat$exon_count_total,dat$shared_exon_length_list)

#output=subset(dat,iri_whole!=0,select=c('gene_symbol', 'chr', 'strand','iri_whole','intron_count','exon_count'))
#write.table(output,paste(Args[6],"ir",sep="_"),quote=F,sep="\t",row.names=F,col.names=T)
