tar -zxvf pgsc_HGDP+1kGP_v1.tar.zst

# GRCh38
# Subset to unrelated individuals, remove non-biallelic and ambiguous SNPs, and filter on genotype missingness, individual missingness, MAF, and HWE.
../plink2 \
    --threads 2 \
    --memory 16000 \
    --pfile GRCh38_HGDP+1kGP_ALL vzs \
    --remove GRCh38_HGDP+1kGP.king.cutoff.out.id \
    --max-alleles 2 \
    --snps-only just-acgt \
    --rm-dup exclude-all \
    --geno 0.1 \
    --mind 0.1 \
    --maf 0.01 \
    --hwe 0.0001 \
    --make-pgen vzs \
    --allow-extra-chr --autosome \
    --out ldref_38


# Further remove ambiguous variants
zstdcat ldref_38.pvar.zst |
awk 'BEGIN {OFS="\t"} !/^#/ {
    if (($4=="A" && $5=="T") || ($4=="T" && $5=="A") || \
        ($4=="C" && $5=="G") || ($4=="G" && $5=="C"))
        print $3
}' - > ambiguous_snps38.txt

../plink2 \
    --threads 2 \
    --memory 16000 \
    --pfile ldref_38 vzs \
    --exclude ambiguous_snps38.txt \
    --make-pgen vzs \
    --autosome \
    --out ldref_38

# Create LD matrices for each autosome (GRCh38)
for chr in {1..22}; do
    ../plink2 \
        --threads 4 \
        --memory 16000 \
        --chr $chr \
        --pfile ldref_38 vzs \
        --r2-unphased zs 'cols=-chrom,-pos,-ref,-maj' --ld-window-kb 10000 --ld-window-r2 0.099 \
        --freq zs 'cols=-nobs' \
        --out ldref_38_$chr
done


# GRCh37
../plink2 \
    --threads 2 \
    --memory 16000 \
    --pfile GRCh37_HGDP+1kGP_ALL vzs \
    --remove GRCh37_HGDP+1kGP.king.cutoff.out.id \
    --max-alleles 2 \
    --snps-only just-acgt \
    --rm-dup exclude-all \
    --geno 0.1 \
    --mind 0.1 \
    --maf 0.01 \
    --hwe 0.0001 \
    --make-pgen vzs \
    --allow-extra-chr --autosome \
    --out test_ref_37

zstdcat test_ref_37.pvar.zst |
awk 'BEGIN {OFS="\t"} !/^#/ {
    if (($4=="A" && $5=="T") || ($4=="T" && $5=="A") || \
        ($4=="C" && $5=="G") || ($4=="G" && $5=="C"))
        print $3
}' - > ambiguous_snps_37.txt

../plink2 \
    --threads 2 \
    --memory 16000 \
    --pfile test_ref_37 vzs \
    --exclude ambiguous_snps_37.txt \
    --make-pgen vzs \
    --autosome \
    --out test_ref_37

for chr in {1..22}; do
    ../plink2 \
        --threads 4 \
        --memory 16000 \
        --chr $chr \
        --pfile test_ref_37 vzs \
        --r2-unphased zs 'cols=-chrom,-pos,-ref,-maj' --ld-window-kb 10000 --ld-window-r2 0.099 \
        --freq zs 'cols=-nobs' \
        --out ldref_37_$chr
done