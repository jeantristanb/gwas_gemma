plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores

process unzipdir{
  input :
     path(vcfzip) 
  publishDir "${params.output_dir}/format/vcfi",  mode:'copy'
  output :
   path("*.vcf.gz"), emit : vcf
   path("*.gz"), emit: gzip
  script :
    passd = (params.unzip_password!="") ? " -P ${params.unzip_password} " : ""
    """
    unzip $passd $vcfzip
    """
}

process computedstat{
 label 'py3utils'
 memory params.plink_mem_req
  time   params.big_time
  input :
     path(vcf) 
  output :
     path("${Ent}") 
  script :
    Ent=vcf.baseName+".stat"
    """
    bcftools query -f '%CHROM %REF %ALT %POS %INFO/${params.score_imp} ${params.statfreq_vcf}\n' $vcf > $Ent
    """
}

process dostat{
 memory params.plink_mem_req
 input :
    path(allstat) 
 publishDir "${params.output_dir}/format/statistics/",  mode:'copy'
 output :
   path("${fileout}*")
 script :
  fileout=params.output_pat+"_report"
  allfile=allstat.join(',')
  """
  stat_vcf.py  --out $fileout --min_score ${params.min_scoreinfo} --list_files $allfile
  """
}

process formatvcfscore{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     tuple path(ref),path(vcf)
  publishDir "${params.output_dir}/format/bcf_filter",  mode:'copy', pattern: "*.bcf"
  output :
     tuple path("${Ent}.bed"), path("${Ent}.bim"), path("${Ent}.fam"), emit: plk 
     path("${Ent}.bim"), emit: bim
     tuple val(Ent), path("${Ent}.bcf"), emit: bcf
  script :
    Ent=vcf.baseName.replaceAll(/.vcf$/, '')
    """
    bcftools view -Ou -i '${params.score_imp}>${params.min_scoreinfo}' $vcf | bcftools norm -Ou -m -any | bcftools norm -Ou -f $ref |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' > $Ent".bcf"
  cat $Ent".bcf" |plink --bcf /dev/stdin \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --make-bed \
    --out ${Ent}
    cp ${Ent}.bim ${Ent}.save.bim
    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
    """
}

process convertbcf_invcf{
  label 'py3utils'
  time   params.big_time
  publishDir "${params.output_dir}/format/vcf_filter",  mode:'copy', pattern: "*.vcf.gz"
  input :
      tuple val(Ent), file(bcf) 
   output :
      val(vcf)
   script :
     vcf=Ent+".vcf.gz"
     """
     bcftools convert $bcf -O z > $vcf
     """
}



process formatvcf{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     tuple path(ref),path(vcf) 
  output :
     tuple path("${Ent}.bed"), path("${Ent}.bim"), path("${Ent}.fam"), emit: plk
     path("${Ent}.bim"), emit: bim
  script :
    Ent=vcf.baseName
    """

    bcftools view -Ou $vcf | bcftools norm -Ou -m -any | bcftools norm -Ou -f $ref |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' |
  plink --bcf /dev/stdin \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --make-bed \
    --out ${Ent}


    cp ${Ent}.bim ${Ent}.save.bim
    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
    """
}

process GetRsDup{
    input :
      path(bim) 
    output :
      tuple path(outdel),path(out) 
    script :
      lbim=bim.join(",")
      out="snpfile_red.rs"
      outdel="snpfile_del.rs"
      """
      search_dup_bim.py $lbim $out $outdel
      """
}

process TransformRsDup{
    input :
     tuple path(delrange),path(rstochange), path(bed),path(bim),path(fam) 
    output :
       tuple path("${newheader}.bed"),path("${newheader}.bim"),path("${newheader}.fam"), emit : plk
       val("$newheader"), emit : plkHead
    script :
       header=bed.baseName
       newheader=header+"_rsqc"
       """
       cp $fam ${header}_save.fam
       cp $bed ${header}_save.bed
       replacers_forbim.py ${bim} $rstochange $header"_save.bim"
       plink -bfile ${header}_save --keep-allele-order --make-bed -out $newheader --exclude range $delrange
       rm ${header}_save*
       """
}

process AddedCM{
    cpus params.max_plink_cores
    memory params.plink_mem_req
    time   params.big_time
    input :
       tuple path(map),path(bedi),path(bimi),path(fami) 
    output :
       tuple path(bedf),path(bimf),path(famf), emit : plk
       val(header),  emit plkHead
    script :
       headeri=bedi.baseName
       header=headeri+"_map"
       bedf=header+".bed"
       bimf=header+".bim"
       famf=header+".fam"
       cm_shap=header+".shape"
       """
       chro=`head $bimi|awk '{print \$1}'|uniq`
       sed '1d' $map|awk -v chro=\$chro '{if(chro==\$1)print \$2"\\t"\$3"\\t"\$4}' >> $cm_shap
       awk '{print \$2}' $headeri".bim" | sort | uniq -d > duplicated_snps.snplist
       plink --bfile $headeri --exclude duplicated_snps.snplist --make-bed --keep-allele-order --out $headeri"_tmp"
       plink --bfile $headeri"_tmp" --list-duplicate-vars ids-only suppress-first
       plink --bfile $headeri"_tmp" --keep-allele-order --cm-map $cm_shap \$chro   --threads ${params.max_plink_cores} --make-bed --out $header  --exclude plink.dupvar
       """

}

process MergePlink{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
       path(lplk)
       //val(hplk)
  output :
     tuple path("${params.output_pat}_merge.bed"), path("${params.output_pat}_merge.bim"),path("${params.output_pat}_merge.fam")
  script :
       hplk2=lplk.join(',')
       """
       echo $hplk2 | awk -F',' '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' | sed 's/\\.[^.]*\$//'  | sort |uniq |sed '1d'> fileplk
       hplkFirst=`echo $hplk2 | awk -F',' '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' | sed 's/\\.[^.]*\$//'  | sort |uniq |head -1`
       plink --bfile \$hplkFirst --keep-allele-order --threads ${params.max_plink_cores} --merge-list fileplk --make-bed --out ${params.output_pat}"_merge"
       """
}

process check_names_plkconvert{
   label 'R'
   input :
       tuple path(bed), path(bim), path(fam)
       path(data)
   publishDir "${params.output_dir}/format/plink/",  mode:'copy'
   output : 
       tuple path("${newplk}.bed"), path("${newplk}.bim"), path("${newplk}.fam")
   script :
     plk=bed.baseName
     newplk=plk+'_idupdate'
     """
     cp $bed $newplk".bed"
     cp $bim $newplk".bim"
     change_names_plkconvert.r $fam $data $newplk".fam"
     """ 
}



workflow format_vcfinplk{
  /*case where file is unzip*/
   if(params.unzip_zip){
     unzipdir(channel.fromPath(params.listfile_vcf))
     list_vcf=unzipdir.out.vcf
   }else{
    if(params.listfile_vcf!=""){
       println("used a list of vcf")
       list_vcf=Channel.fromPath(file(params.listfile_vcf,  checkIfExists:true).readLines(), checkIfExists:true)
    }else if(params.file_vcf!=''){
       list_vcf = Channel.fromPath(params.file_vcf, checkIfExists:true)
    }
   }
   if(params.do_stat){
     computedstat(list_vcf)
     dostat(computedstat.out.collect())
   }
 if(params.reffasta=="" || params.reffasta==true){
  println("to format file need a fasta file : args --reffasta null")
      exit 1
 }
  ref_ch=Channel.fromPath(params.reffasta, checkIfExists:true)
  if(params.min_scoreinfo>0){
    formatvcfscore(ref_ch.combine(list_vcf))
    convertbcf_invcf(formatvcfscore.out.bcf)
    GetRsDup(formatvcfscore.out.bim.collect())
    TransformRsDup(GetRsDup.out.combine(formatvcfscore.out.plk))
  }else{
   formatvcf(ref_ch.combine(list_vcf))
   GetRsDup(formatvcf.out.bim.collect())
   TransformRsDup(GetRsDup.out.combine(formatvcf.out.plk))
 } 
  if(params.genetic_maps!=""){
    AddedCM(TransformRsDup.out.plk)
    listplink=AddedCM.out.plk.collect()
    headplink=AddedCM.out.plkHead.first()
  }else{
    listplink=TransformRsDup.out.plk.collect()
    //headplink=TransformRsDup.out.plkHead.collect()
  }
  if(params.listfile_vcf!=""){
    MergePlink(listplink)
    plk_merge=MergePlink.out
  }else{
    plk_merge=listplink
  }
  check_names_plkconvert(plk_merge, channel.fromPath(params.data, checkIfExists:true)) 
  emit :
    plk = check_names_plkconvert.out
 }  

process formatvcfinbimbam{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     tuple val(chro), path(vcf)
  publishDir "${params.output_dir}/format/bimbam", mode:'copy'
  output :
     tuple val(chro),path("${Ent}.bimbam"), path("${fileind}")
  script :
    headvcf=vcf.baseName
    Ent=(chro!=-1) ? "${headvcf}_${chro}" :  "$headvcf"
    chroparam=(chro!=-1) ?  " --regions $chro" : ""
    fileind=Ent+".ind"
    """
    zcat $vcf |head -10000|grep "#"|tail -1| awk '{for(Cmt=10;Cmt<=NF;Cmt++)print \$Cmt}' > $fileind
    bcftools index -f $vcf
    ${params.bcftools_bin} view -i '${params.score_imp}>${params.min_scoreinfo}' $chroparam $vcf |${params.qctoolsv2_bin} -g - -vcf-genotype-field ${params.genotype_field} -ofiletype bimbam_dosage -og ${Ent}.bimbam -filetype vcf
    """
}

process formatvcfinbimbam_ind{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     tuple val(chro), path(vcf), path(fileind)
  publishDir "${params.output_dir}/format/bimbam", mode:'copy'
  output :
     tuple val(chro),path("${Ent}.bimbam"), path("${fileind}")
  script :
    headvcf=vcf.baseName
    Ent=(chro!=-1) ? "${headvcf}_${chro}" :  "$headvcf"
    chroparam=(chro!=-1) ?  " --regions $chro" : ""
    fileind=Ent+".ind"
    """
    zcat $vcf |head -10000|grep "#"|tail -1| awk '{for(Cmt=10;Cmt<=NF;Cmt++)print \$Cmt}' > $fileind
    keep_vcfiid.r --data $fileind --vcf $vcf --out sample_vcf.keep
    bcftools index $vcf
    ${params.bcftools_bin} view --samples-file sample_vcf.keep -i '${params.score_imp}>${params.min_scoreinfo}' $chroparam $vcf |${params.qctoolsv2_bin} -g - -vcf-genotype-field ${params.genotype_field} -ofiletype bimbam_dosage -og ${Ent}.bimbam -filetype vcf
    """
}


