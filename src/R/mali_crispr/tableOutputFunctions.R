writeProbeRanksFile<-function(probes,p_rank){
  resp<-data.frame(probes,p_rank)
  # TODO: refactor file suffix, separator into easier to manage location
  write.table(resp,file=paste(project,"_p.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
}

writeSingleGeneFitnessesFile<-function(genes, f, f_sd){
  resy<-data.frame(genes,f,f_sd)
  # TODO: refactor col names, file suffix, separator into easier to manage location
  colnames(resy)<-c("gene","f","sd")
  write.table(resy[-1,],file=paste(project,"_f.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
}

writeConstructFitnessesFile<-function(pA_pB, fc, sdfc){
  resp<-data.frame(pA_pB,fc,sdfc)
  # TODO: refactor file suffix, separator into easier to manage location
  write.table(resp,file=paste(project,"_fc.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)  
}

writePosteriorProbabilityFile<-function(pi_mean, pi_sd, pi_iter_null, g1names, g2names, ggnames, uutri, f, fdr_left,fdr_right){
  z<-pi_mean/sd(pi_mean)
  
  #pi_iter_null<-pi_iter-pi_mean
  abspi<-abs(pi_mean)
  PP<-apply(abs(pi_iter_null)<abspi,1,mean)
  
  oPP<-order(z)
  
  #get names of these gene pairs
  # These g1//g2, etc hardcodes are ok because represent DUAL crispr genes
  names_of_g1<-g1names[uutri]
  fg1<-f[names_of_g1]
  names_of_g2<-g2names[uutri]
  fg2<-f[names_of_g2]
  fg12<-fg1+fg2
  names_of_gg<-ggnames[uutri]
  
  res<-data.frame(names_of_gg,names_of_g1,fg1,names_of_g2,fg2,fg12,pi_mean,pi_sd,PP,abspi,fdr_left,fdr_right,z)
  # TODO: Let's refactor these colnames, file suffix, and separator somewhere easier to manage
  colnames(res)<-c("gene_gene","geneA","fA","geneB","fB","fA+fB","pi","sd","PP","abs pi","FDR left","FDR right","z")
  #TODO: remove format call here; added to aid in testing
  write.table(format(res[oPP,], digits=14),file=paste(project,"_pi.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)  
}