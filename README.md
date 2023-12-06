# ddRAD_Apodemus

narr_scaffold.rb

<- apo.final.assembly.fasta

<- ddRAD_EnzymeSelect.csv

-> tgtfragment.oj

 
parse_augustus.rb

<- Apo_augustus_ver2.gff

-> scaffcds.oj


make_overlap_refrag.rb 

<- tgtfragment.oj

<- scaffcds.oj

-> cds_frg_all.oj

overlap_refrag.rb

<- cds_frg_all.oj

-> rereg_hash_all.oj

-> overlap_set.oj

-> overlap_len_all.oj
