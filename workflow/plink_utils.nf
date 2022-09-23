process plinkextractpos{
 cpus max_plink_cores
 input:
  tuple path(bed), path(bim), path(fam)
  path(listsnp)
  path(indkeep)
 output :
  tuple path("${bfilef}.bed"), path("${bfilef}.bim"), path("${bfilef}.fam")
 script :
   bfile=bed.baseName
   bfilef=bfile+'_relsub'
   """
   sed '1d' $indkeep | awk '{print \$1\"\\t\"\$2}' > indkeep
   plink -bfile $bfile  --extract range $listsnp --keep indkeep --keep-allele-order --make-bed -out $bfilef
   """
}


process plinkextractind{
 cpus max_plink_cores
 input :
  tuple path(bed), path (bim), path(fam)
  path(indkeep)
 output :
  tuple path("${outbed}.bed"), path("${outbed}.bim"), path("${outbed}.fam"), emit: filterind
 script :
  bfile=bed.baseName
  outbed=bfile+"_subind"
  """
  sed '1d' $indkeep | awk '{print \$1\"\\t\"\$2}' > indkeep
  plink -bfile $bfile --keep-allele-order --make-bed -out $outbed --threads $max_plink_cores --keep $indkeep
  """
}


process subsample_snps{
 cpus max_plink_cores
 input:
  tuple path(bed), path(bim), path(fam)
  path(snp_exclude_bed)
  path(snp_include_bed)
  path(indkeep)
 output:
  path("${outfile}.prune.in"), emit: subsample_snps_list
 script:
  bfile=bed.baseName
  outfile="sub_indep_pairwise"
  rangeexclude=(params.snps_exclude_rel=="") ? "" : " --exclude range $snp_exclude_bed "
  rangeinclude=(params.snps_include_rel=="") ? "" : " --extract range $snp_include_bed "
  maf=(params.cut_maf_rel=="") ? "" : " --maf ${params.cut_maf_rel}"
  balisethin=(params.thin_snp_rel=="") ? "0" : "1"
  """
  sed '1d' $indkeep | awk '{print \$1\"\\t\"\$2}' > indkeep
  plink --bfile $bfile --threads $max_plink_cores --autosome $rangeexclude --indep-pairwise ${params.plink_indep_pairwise} --out $outfile $maf --keep indkeep --keep-allele-order $rangeinclude
  if [ "$balisethin" == "1" ]
  then
   cp ${outfile}.prune.in ${outfile}.prune.tmp.in
   shuf ${outfile}.prune.tmp.in | head -${params.thin_snp_rel} > ${outfile}.prune.in
  fi
  """
}

process computeN_plink{
  label 'R'
  cpus params.max_plink_cores
  memory { strmem(other_mem_req) + 5.GB * (task.attempt -1) }
  maxRetries 10
  input :
    tuple path(data), path(bed), path(bim), path(fam), val(this_pheno),val(covar)
  output :
     tuple val(our_pheno2),file("${headout}.frq")
  script :
    our_pheno2          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
    our_pheno          = this_pheno.replaceAll(/[0-9]+@@@/,"")
    headout=our_pheno2+'_statn'
    plkf=bed.baseName
    covar = (covar=="") ? "" : " --covar $covar "

    """
    formatpheno_plink.r --data $data --pheno ${our_pheno} --out pheno_plink --binary 0  $covar
    plink -bfile $plkf --keep pheno_plink --freq -out $headout"_tmp" --keep-allele-order --threads ${params.max_plink_cores}
    merge_freqandbim.py  --freq  ${headout}_tmp.frq --bim $bim --out ${headout}.frq
    """
}


