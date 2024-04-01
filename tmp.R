# Compute FGr

rule compute_FGr_chr:
  input:
  r="data/pga_paper/{rep}/{nsnp}/r/{contrasts}_chr{chr}.rvec",
  gp_genos="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr{chr}_v3.pgen",
  overlap_snps="data/pga_paper/{rep}/{nsnp}/variants/overlappingSNPs_chr{chr}.txt",
  IDs="data/pga_paper/{rep}/ids/gwas.ids"
  output:
    "output/calculate_FGr/pga_paper/{rep}/{nsnp}/FGr_{contrasts}_{chr}.txt"
  resources:
    mem_mb=100000,
  time="03:00:00"
  threads: 16
  params:
    gp_prefix = "/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr{chr}_v3",
  out_prefix = "output/calculate_FGr/pga_paper/{rep}/{nsnp}/{contrasts}_{chr}"
  shell:
    """
        Rscript code/calculate_FGr/calc_FGr_chr.R {params.gp_prefix} {params.out_prefix} {input.r} {input.overlap_snps} {output} {input.IDs}
        """

  rule concat_chr_FGr:
    input:
    expand("output/calculate_FGr/pga_paper/{{rep}}/{{nsnp}}/FGr_{{contrasts}}_{chr}.txt", chr = CHR)
  output:
    "output/calculate_FGr/pga_paper/{rep}/{nsnp}/FGr_{contrasts}.txt"
  shell:
    """
        Rscript code/calculate_FGr/concat_FGr.R {output} {input}
        """

  # Calculate Error in FGr

  rule compute_FGr_error:
    input:
    expand("output/calculate_FGr/pga_paper/{{rep}}/{{nsnp}}/FGr_{{contrasts}}_{chr}.txt", chr = CHR)
  output:
    "output/calculate_FGr/pga_paper/{rep}/{nsnp}/Var_{contrasts}.txt"
  params:
    prefix_in = "output/calculate_FGr/pga_paper/{rep}/{nsnp}/FGr_{contrasts}"
  shell:
    """
        Rscript code/calculate_FGr/compute_error_jacknife.R {params.prefix_in} {output}
        """

  # GWAS panel PCA

  rule make_tmp_plink2:
    input:
    gp="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr{chr}_v3.pgen",
  overlap_snps="data/ukbb/ukbb_pc_snps.txt",
  IDs="data/pga_paper/{rep}/ids/gwas.ids"
  output:
    "/scratch/jgblanc/ukbb/plink2-files/tmp/{rep}/tmp_{chr}.pgen",
  "/scratch/jgblanc/ukbb/plink2-files/tmp/{rep}/tmp_{chr}.pvar",
  "/scratch/jgblanc/ukbb/plink2-files/tmp/{rep}/tmp_{chr}.psam"
  params:
    prefix_out = "/scratch/jgblanc/ukbb/plink2-files/tmp/{rep}/tmp_{chr}",
  prefix_plink = "/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr{chr}_v3"
  threads: 16
  resources:
    mem_mb=38000,
  time="00:30:00"
  shell:
    """
        plink2 --pfile {params.prefix_plink} \
        --keep {input.IDs} \
        --extract {input.overlap_snps} \
	--make-pgen \
        --out {params.prefix_out}
        """

  rule merge_GWAS_PCA:
    input:
    gp_genos=expand("/scratch/jgblanc/ukbb/plink2-files/tmp/{rep}/tmp_{chr}.pgen", chr=CHR)
  output:
    vec="data/pga_paper/{rep}/pca/pca.eigenvec",
  val="data/pga_paper/{rep}/pca/{gwas}/pca.eigenval"
  params:
    prefix_out = "data/pga_paper/{rep}/pca/pca",
  merge_prefix = expand("/scratch/jgblanc/ukbb/plink2-files/tmp/{{rep}}/tmp_{chr}", chr=CHR, newline="\n"),
  remove_path = "/scratch/jgblanc/ukbb/plink2-files/tmp/{rep}/"
  threads: 16
  resources:
    mem_mb=38000,
  time="03:00:00"
  shell:
    """
	echo {params.merge_prefix} > {wildcards.rep}-tmp.txt
        tr ' ' '\n' < "{wildcards.rep}-tmp.txt" > "{wildcards.rep}-tmp_chrm_list.txt"
        plink2 --pmerge-list {wildcards.rep}-tmp_chrm_list.txt \
        --pca 40 approx \
	--memory 50000 \
	--threads 16 \
        --out {params.prefix_out}
        rm {wildcards.rep}-tmp*
	rm {params.remove_path}*
	rm {params.prefix_out}.p*
        """


  # Assemble covariates

  rule assemble_fixed_covars:
    input:
    sex="data/phenotypes/genetic_sex_22001.txt",
  batch="data/phenotypes/genotype_measurement_batch_22000.txt"
  output:
    "output/run_gwas/fixed_covars.txt"
  shell:
    """
        Rscript code/run_gwas/assemble_fixed_covar.R {input.sex} {input.batch} {output}
        """

  rule assemble_quant_covars_FGr:
    input:
    age="data/phenotypes/age_at_recruitment_21022.txt",
  covar="output/calculate_FGr/pga_paper/{rep}/{nsnp}/FGr_{contrasts}.txt",
  pcs="data/pga_paper/{rep}/pca/pca.eigenvec"
  output:
    "output/run_gwas/pga_paper/{rep}/{nsnp}/covariates/{covar}_{contrasts}_{chr}-PC{pc}.txt"
  params:
    pcNum = lambda wildcards: get_chr_num(wildcards.pc)
  shell:
    """
        Rscript code/run_gwas/assemble_quant_covar.R {input.age} {input.covar} {input.pcs} {params.pcNum} {output}
        """

  rule assemble_quant_covars_noFGr:
    input:
    age="data/phenotypes/age_at_recruitment_21022.txt",
  pcs="data/pga_paper/pca/{rep}/pca.eigenvec"
  output:
    "output/run_gwas/pga_paper/{rep}/{nsnp}/covariates/noFGr_{contrasts}_{chr}-PC{pc}.txt"
  params:
    pcNum = lambda wildcards: get_chr_num(wildcards.pc)
  shell:
    """
        Rscript code/run_gwas/assemble_quant_covar_control.R {input.age} {input.pcs} {params.pcNum} {output}
        """

  # Run GWAS

  rule run_fastGWA:
    input:
    pheno="data/phenotypes/{phenotype}.txt",
  qcovar="output/run_gwas/pga_paper/{rep}/{nsnp}/covariates/{covar}_{contrasts}_{chr}-PC{pc}.txt",
  gp_genos="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr{chr}_v3.pgen",
  overlap_snps="data/pga_paper/{rep}/{nsnp}/variants/overlappingSNPs_chr{chr}.txt",
  covar="output/run_gwas/fixed_covars.txt",
  IDs="data/pga_paper/{rep}/ids/gwas.ids
    output:
        "/scratch/jgblanc/stratification-data_analysis/output/run_gwas/pga_paper/{rep}/{nsnp}/{phenotype}/raw_LR/{covar}_{contrasts}_{chr}-PC{pc}.fastGWA",
    params:
        prefix_out = "/scratch/jgblanc/stratification-data_analysis/output/run_gwas/pga_paper/{rep}/{nsnp}/{phenotype}/raw_LR/{covar}_{contrasts}_{chr}-PC{pc}",
        prefix_plink = "/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr{chr}_v3"
    threads: 16
    resources:
        mem_mb=100,
        time="00:30:00"
    shell:
        """
  gcta  --pheno {input.pheno} \
  --keep {input.IDs} \
  --qcovar {input.qcovar} \
  --covar {input.covar} \
  --extract {input.overlap_snps} \
  --fastGWA-lr \
  --threads {threads} \
  --pfile {params.prefix_plink} \
  --out {params.prefix_out}
  """

# Pick minimum p-value per block

rule assign_snps_to_blocks:
    input:
        ldBlocks="data/LD_blocks/fourier_ls-all_parsed.bed",
        ss="/scratch/jgblanc/stratification-data_analysis/output/run_gwas/pga_paper/{rep}/{nsnp}/{phenotype}/raw_LR/{covar}_{contrasts}_{chr}-PC{pc}.fastGWA",
        r="data/pga_paper/{rep}/{nsnp}/r/{contrasts}_chr{chr}.rvec",
    output:
        betas="/scratch/jgblanc/stratification-data_analysis/output/pga_test/pga_paper/{rep}/{nsnp}/{phenotype}/formated_ss_LR/{covar}_{contrasts}_{chr}-PC{pc}.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        """
  Rscript code/pga_test/assign_snps_to_blocks_LR.R {input.ldBlocks} {input.ss} {input.r} {output.betas} {params.pt}
  rm {input.ss}
  """

rule concat_snps:
    input:
        expand("/scratch/jgblanc/stratification-data_analysis/output/pga_test/pga_paper/{{rep}}/{{nsnp}}/{{phenotype}}/formated_ss_LR/{{covar}}_{{contrasts}}_{chr}-PC{{pc}}.betas", chr=CHR)
    output:
        "output/pga_test/pga_test/{rep}/{nsnp}/{phenotype}/formated_ss_LR/{covar}_{contrasts}.all-PC{pc}.betas"
    shell:
        """
  cat {input} > {output}
  rm {input}
  """










#rule get_ID_lists_A:
#    input:
#        sex="data/phenotypes/genetic_sex_22001.txt",
#        batch="data/phenotypes/genotype_measurement_batch_22000.txt",
#        north="data/phenotypes/PlaceOfBirthNorthCord_129.txt",
#        east="data/phenotypes/PlaceOfBirthEastCord_130.txt",
#        age="data/phenotypes/age_at_recruitment_21022.txt",
#        genotyped="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_imp_chr22_v3.psam",
#        wbs="data/ukbb/WBS.ids"
#    output:
#        gwas="data/pga_paper/{rep}/ids/gwas.ids",
#        test="data/pga_paper/{rep}/ids/test.ids"
#    params:
#        gs=100000,
#        ts=100000
#    shell:
#        """
  #        Rscript code/PGA_paper/get_IDs_A.R {input.sex} {input.batch} {input.north} {input.east} {input.age} {input.genotyped} {output.gwas} {output.test} {params.gs} {params.ts} {input.wbs}
  #        """
