## mm9, GRCm38

rsync -azvp aquevedo@172.22.238.48:/data/References/AgilentSureSelect/Mouse_Exome_V1/S0276129_{Regions.bed,Targets.txt} \
~/genomes/mouse/mm38/

# Check chr names mathc UCSC style
cut -f 1 ~/genomes/mouse/mm38/S0276129_Regions.bed | uniq -c

## mm10, GRCm39
# Copy file (mouse genome)
rsync -azvp aquevedo@172.22.238.48:/data/References/AgilentSureSelect/Mouse_Exome_V1/liftOver_mm10/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz* \
~/genomes/mouse/mm39/

# Check chr names
zcat ~/genomes/mouse/mm39/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz | head -n 20

# add chr to mathc UCSC style
zcat ~/genomes/mouse/mm39/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz | \
awk '\
{if ( $0 !~ /^#/ )
	if ( $0 ~ /^MT/ ){
		r = gensub(/^MT(.+)/, "chrM\\1", "g", $0 );
		print r;
	} else {
		print "chr"$0;
	}
}' | bgzip > /genomes/mouse/mm39/MouseExome_mm10_liftover_Covered_UCSC-chr-names_sorted.bed.bgz

