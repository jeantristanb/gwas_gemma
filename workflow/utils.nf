process getListeChro{
  input :
   tuple path(bed), path(bim), path(fam)
  output :
     stdout
  script:
     """
     cat $bim|awk '{print \$1}'|uniq|sort|uniq
     """
}




/*process insure that file is in chr pos pos rs*/
process checkposrsfile{
  input :
    path(pos)
    tuple path(bed), path(bim), path(fam)
    val(out)
  output :
    path("$out")
  script :
     """
     check_filpos.py $pos $bim $out
     """
}

/*check snp*/
workflow getsnpbuilrelat{
 take :
   ch_plkfile
 main :
   checkposrsfile(channel.fromPath(params.listsnps_buildrelat, checkIfExists:true), ch_plkfile,'snpbuildrelat.bed')
 emit :
   pos_chr=checkposrsfile.out
}
workflow getsnpexcluderelat{
  take:
   ch_plkfile
  main:
  if(params.snps_exclude_rel!="")ch_snps_exclude_rel=checkposrsfile(channel.fromPath(params.snps_exclude_rel,checkIfExists:true),  ch_plkfile,'snpexcluderelat.bed').out
   else ch_snps_exclude_rel=channel.fromPath("${dummy_dir}/00")
  emit :
   pos_chr=ch_snps_exclude_rel
}

workflow getsnpincluderelat{
  take :
   ch_plkfile
  main:

  if(params.snps_include_rel!=""){
     checkposrsfile(channel.fromPath(params.snps_include_rel,checkIfExists:true),  ch_plkfile,'snpincluderelat.bed')
     ch_snps_include_rel=checkposrsfile.out
   }else ch_snps_include_rel=channel.fromPath("${dummy_dir}/01")
  emit :
   pos_chr=ch_snps_include_rel
}

