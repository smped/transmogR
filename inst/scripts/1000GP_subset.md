# Creation of example vcf: 1000GP_subset.vcf

The source files were originally taken from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV.
The VCF included here is a very small subset of 100 variants with an AF > 0.5 amongst the unrelated samples.
The dummy header was first created using 

```bash
bcftools view -h /data/references/1000GP_NYGC_30x_GRCh38_Phased//1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz > original.vcf
```

This file was then filtered using a loop to step through all individual chroosome files

```bash
zcat ${CURFILE} |\
	  bcftools view -G -H -i '(INFO/AC_EUR_unrel + INFO/AC_EAS_unrel + INFO/AC_AMR_unrel + INFO/AC_SAS_unrel + AC_AFR_unrel) > 0.5*(INFO/AN_EUR_unrel + INFO/AN_EAS_unrel + INFO/AN_AMR_unrel + INFO/AN_SAS_unrel + AN_AFR_unrel)' |\
	  awk 'length($4) <= 50 && length($5) <= 50' |\
	  egrep -v 'HGSV' >> ${OUTFILE}
```

This final file was then subset to only the first 100 variants using

```bash
bcftools annotate -x INFO ~/create_variant_sets/vcf/1000GP_SNV_INDEL_panhuman_unrel_0.5.vcf.gz | head -n133 > inst/extdata/1000GP_subset.vcf
bgzip inst/extdata/1000GP_subset.vcf
tabix inst/extdata/1000GP_subset.vcf.gz
```
