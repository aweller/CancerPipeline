config = {"test_Quasar_A1_VC1_VF01": 
              dict(root = "/home/andreas/bioinfo/core/test_ruffus/data/",
                    bed = "/home/andreas/bioinfo/core/test_ruffus/data/core/roi_148gene_panel_HP.bed",
                    ref = "/home/andreas/bioinfo/core/general/data/hg19.fasta",
                    raw_vcf_folder = "/home/andreas/bioinfo/core/test_ruffus/data/vcfs_A1_VC1_VF00/",
                    max_snv = 50,
                    min_cov = 400,
                    min_varfreq = 0.05),
            
            "test_Quasar_A1_VC1_VF02":   
            dict(root = "/home/andreas/bioinfo/core/test_ruffus/data/",
                    bam_folder = "/home/andreas/bioinfo/core/test_ruffus/data/bams/",
                    bed = "/home/andreas/bioinfo/core/test_ruffus/data/core/roi_148gene_panel_HP.bed",
                    ref = "/home/andreas/bioinfo/core/general/data/hg19.fasta",
                    raw_vcf_folder = "/home/andreas/bioinfo/core/test_ruffus/data/vcfs_A1_VC1_VF00/",
                    max_snv = 50,
                    min_cov = 100,
                    min_varfreq = 0.01),
            
            "Quasar_A1_VC1_VF01":   
            dict(root = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/",
                    raw_vcf_folder = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/vcfs_A1_VC1_VF00/",
                    bam_folder = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/data/quasar_bams/",
                    bed = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/core/roi_148gene_panel_HP.bed",
                    ref = "/home/andreas/bioinfo/core/general/data/hg19.fasta",
                    blacklist = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/vcfs_A1_VC1_VF00/deaminated_samples.txt",
                    max_snv = 50,
                    min_cov = 400,
                    min_varfreq = 0.05),
            
            "Quasar_A1_VC1_VF02":   
            dict(root = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/",
                    raw_vcf_folder = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/vcfs_A1_VC1_VF00/",
                    bam_folder = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/data/quasar_bams/",
                    bed = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/core/roi_148gene_panel_HP.bed",
                    ref = "/home/andreas/bioinfo/core/general/data/hg19.fasta",
                    whitelist = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/test_whitelist_134.txt",
                    version_numbers_not_in_blacklist = True,
                    max_snv = 50,
                    min_cov = 400,
                    min_varfreq = 0.05),

            "Quasar_A1_VC1_VF03":   
            dict(root = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/",
                    raw_vcf_folder = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/vcfs_A1_VC1_VF00/",
                    bam_folder = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/data/quasar_bams/",
                    bed = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/core/roi_148gene_panel_HP.bed",
                    ref = "/home/andreas/bioinfo/core/general/data/hg19.fasta",
                    whitelist = "/home/andreas/bioinfo/projects/wtc_quasar_pipeline/data/test_whitelist_134.txt",
                    version_numbers_not_in_blacklist = True,
                    cpus = 1, 
                    vcf_type = "iontorrent",
                    min_cov = 100,
                    min_varfreq = 0.05),
            
            "CLL_A1_VC1_VF01":   
            dict(root = "/home/andreas/bioinfo/projects/nhs_cll/data/",
                    raw_vcf_folder = "/home/andreas/bioinfo/projects/nhs_cll/data/vcfs_CLL_A1_VC1_VF00/",
                    #whitelist = "/home/andreas/bioinfo/projects/nhs_cll/data/vcfs_CLL_A1_VC1_VF00/test_whitelist_5.txt",
                    no_filtering = True,
                    vcf_type = "illumina_strelka",
                    bed = None,
                    ref = None),  

            "CLL_A1_VC1_VF02":   
            dict(root = "/home/andreas/bioinfo/projects/nhs_cll/data/",
                    raw_vcf_folder = "/home/andreas/bioinfo/projects/nhs_cll/data/vcfs_CLL_A1_VC1_VF00/",
                    #whitelist = "/home/andreas/bioinfo/projects/nhs_cll/data/vcfs_CLL_A1_VC1_VF00/test_whitelist_5.txt",
                    vcf_type = "illumina_strelka",
                    bed = None,
                    ref = None,
                    min_cov = 10,
                    min_qual = 30) 

}