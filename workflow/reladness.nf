workflow subsample_snp_rel{
 take :
  ch_plkfile
  ch_pheno
 main :
  if(params.listsnps_buildrelat!=""){
    getsnpbuilrelat(ch_plkfile)
    snpfilers=getsnpbuilrelat.out.posfile
  }else{
   getsnpexcluderelat(ch_plkfile)
   getsnpincluderelat(ch_plkfile)
   subsample_snps(ch_plkfile, getsnpexcluderelat.out.pos_chr,getsnpincluderelat.out.pos_chr,ch_pheno)
   checkposrsfile(subsample_snps.out.subsample_snps_list, ch_plkfile,'list_posrs_rel')
   snpfilers=checkposrsfile.out
  }
 plinkextractpos(ch_plkfile, snpfilers, ch_pheno)
 emit:
   bed_pos_rel=snpfilers
   plk_rel=plinkextractpos.out
}


