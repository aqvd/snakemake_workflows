rsync aquevedo@172.22.238.48:/data/References/AgilentSureSelect/Mouse_Exome_V1/liftOver_mm10/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz* \
~/genomes/mouse/mm38/

zcat ~/genomes/mouse/mm38/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz | head -n 20
zcat ~/genomes/mouse/mm38/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz | \
awk '\
{if ( $0 !~ /^#/ )
	if ( $0 ~ /^MT/ ){
		r = gensub(/^MT(.+)/, "chrM\\1", "g", $0 );
		print r;
	} else {
		print "chr"$0;
	}
}' | cut -f 1 | uniq -c