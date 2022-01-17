# Add "chr" to match UCSC style chr names
# compress using bcftools "-O z" -> compressed vcf
bcftools view /home/aquevedo/genomes/mouse/mm38/dbSNP_all.chr.vcf.gz | \
awk '
{
	# add chr to mathc UCSC style to lines not in header
	if ($0 !~ /^#/) {
		# Filter strange symbols like "-" or "(Total...[0-9]+) in REF ALT cols"	
		if ($4 !~ /[^ACTGN,]/ && $5 !~ /[^,ACTGN*]/) { 
			print "chr"$0;
		}
	# add "chr" to contig field header
	} else if ($0 ~ /^..contig=.ID=.+/) {
		r = gensub(/^(..contig=.ID=)(.+)/, "\\1chr\\2", "g", $0);
		print r;		
	} else {
	# rest of the lines unmodified
		print $0;
	}
}' | 
bcftools view --threads 6 -O z - > /home/aquevedo/genomes/mouse/mm38/dbSNP_UCSC-chr-names.vcf.gz

# Index using gatk, ensure that it will work later with BSQR
gatk IndexFeatureFile -I /home/aquevedo/genomes/mouse/mm38/dbSNP_UCSC-chr-names.vcf.gz

# index using tabix
# tabix -p vcf /home/aquevedo/genomes/mouse/mm38/dbSNP_UCSC-chr-names.vcf.gz
