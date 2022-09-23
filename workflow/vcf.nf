 process get_chrovcf{
       errorStrategy { task.exitStatus in 142..144 ? 'retry' : 'terminate' }
       maxRetries 5

       input:
          path(vcf)
       output :
           tuple env(chro), path(vcf), emit : chro_vcf
       script :
          """
          chro=`zcat $vcf|head -1000|grep -v "#"|awk '{print \$1}'|uniq`
          """
}

process cleanvcf{
          label 'py3utils'
          errorStrategy { task.exitStatus in 142..144 ? 'retry' : 'terminate' }
          maxRetries 5
           input :
            tuple path(filevcf), path(IndTokeep)
          publishDir "${params.output_dir}/format/vcffilter", mode:'copy'
          output :
            path(newfilevcf)
          script :
            newfilevcf='filt_'+filevcf
            indkeep=(params.keep_vcf=="") ? "" : " --keep $IndTokeep "
            """
            ${params.vcfftools_bin} --gzvcf $filevcf --maf ${params.cut_maf} $indkeep  --recode --recode-INFO-all  --stdout | bgzip -c > $newfilevcf
            """
}

workflow cleanvcfwf{
         main :
         if(params.file_vcf!='' & params.listfile_vcf!=''){
          println "file_vcf != '' and listfile_vcf != ''"
          System.exit(-2);
         }
         if(params.file_vcf!=''){
          filevcf=channel.fromPath(params.file_vcf, checkIfExists:true)
          if(params.keep_vcf!=''){
            cleanvcf(filevcf.combine(channel.fromPath(params.keep_vcf, checkIfExists:true)))
            filevcf=cleanvcf.out
          }
         }else {
          filevcf=channel.fromPath("${dummy_dir}/02", checkIfExists:true)
         }
         if(params.listfile_vcf!=''){
            listfilevcf=channel.fromPath(file(params.listfile_vcf, checkIfExists:true).readLines(), checkIfExists:true)
            if(params.keep_vcf!=''){
              cleanvcf(listfilevcf.combine(channel.fromPath(params.keep_vcf, checkIfExists:true)))
              listfilevcf = cleanvcf.out
            }
          }else{
            listfilevcf=channel.fromPath("${dummy_dir}/03", checkIfExists:true)
          }
         emit :
          filevcf = filevcf
          listfilevcf= listfilevcf
}

