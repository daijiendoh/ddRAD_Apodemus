require 'csv'
require 'bio'
require 'oj'

tgtfragment=Oj.load_file("tgtfragment.oj",:mode=>:compat)
#tgtfragment[k][reset]=tfrag
p "end of import tgtfragment"
scaffcds=Oj.load_file("scaffcds.oj",:mode=>:compat)

# p scaffcds
tgt_sets=Hash.new
i=0
scaffcds.each{|scf,d1|

		tgt_sets[scf]=Hash.new
		if tgtfragment[scf] then
			tgt_sets[scf]["cds"]=d1
			tgt_sets[scf]["rfrg"]=tgtfragment[scf]
	
		end
}
Oj.to_file("cds_frg_all.oj",tgt_sets,:mode=>:compat)
# tgtfragment.each{|scf,d1|
# 	p scf
# 	d1.each{|reset,d2|
# 		# p reset
# 		# p d2
# 	}
# }