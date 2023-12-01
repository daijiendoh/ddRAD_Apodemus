require 'csv'
require 'bio'
require 'oj'

# tgtfragment=Oj.load_file("tgtfragment.oj",:mode=>:compat)
# #tgtfragment[k][reset]=tfrag
# p "end of import tgtfragment"
# scaffcds=Oj.load_file("scaffcds.oj",:mode=>:compat)

tgt_rgs=Oj.load_file("cds_frg_all.oj",:mode=>:compat)

cdreg_hash=Hash.new
cdidx_scf=Hash.new
rereg_hash=Hash.new
tgt_rgs.each{|scf,d1|
	cdreg_hash[scf]=Hash.new
	rereg_hash[scf]=Hash.new
	cdidx_scf_arr=Array.new
	d1["cds"].each_with_index{|crg,idx|
		cdidx_scf_arr << idx
		cdreg_hash[scf][["cs",crg[0]]]=[idx,crg[1]-crg[0]]
		cdreg_hash[scf][["ce",crg[1]]]=[idx,crg[1]-crg[0]]
	}
	cdidx_scf[scf]=cdidx_scf_arr
	d1["rfrg"].each{|reset,rrg|
		rereg_hash[scf][reset]=Hash.new
		rrg.each_with_index{|rrg,idx|
			rereg_hash[scf][reset][["rs",rrg[0]]]=[idx,rrg[1]-rrg[0]]
			rereg_hash[scf][reset][["re",rrg[1]]]=[idx,rrg[1]-rrg[0]]
		}
	}
}
Oj.to_file("cdreg_hash_all.oj",cdreg_hash,:mode=>:compat)
Oj.to_file("rereg_hash_all.oj",rereg_hash,:mode=>:compat)

flat_poss=Hash.new
tgt_rgs.each{|scf,d1|
	flat_poss[scf]=Hash.new
	d1["rfrg"].each{|reset,d2|
		
		posarr=Array.new
		d2.each{|refr|
			posarr << ["rs",refr[0]]
			posarr << ["re",refr[1]]
		}
		d1["cds"].each{|crg|
			posarr << ["cs",crg[0]]
			posarr << ["ce",crg[1]]
		}
		flat_poss[scf][reset] =posarr.sort_by{|x| x[1]}

	}
}

overlap_set=Hash.new
flat_poss.each{|scf,d1|
	overlap_set[scf]=Hash.new
	d1.each{|reset,d2|
		overlap_set[scf][reset]=Array.new
		d2.each_with_index{|pos,idx|
			if idx>0 then
				p1=d2[idx-1]
				p2=pos
				p3=d2[idx+1]
				p4=d2[idx+2]
				if ((p1[0]=="cs") && (p2[0]=="rs")) || ((p1[0]=="cs") && (p2[0]=="re")) || ((p1[0]=="rs") && (p2[0]=="ce")) || ((p1[0]=="rs") && (p2[0]=="cs")) then
					
					if ((p1[0]=="cs") && (p2[0]=="re")) ||  ((p1[0]=="rs") && (p2[0]=="ce")) then
						if p2[1]-p1[1] > 10 then
							overlap_set[scf][reset] << [p1,p2]
						end
					elsif (p1[0]=="cs" && p2[0]=="rs" && p3[0]=="re" && p4[0]=="ce") || (p1[0]=="rs" && p2[0]=="cs" && p3[0]=="ce" && p4[0]=="re") then
						# p "#{p1} #{p2} #{p3} #{p4}"
						overlap_set[scf][reset] << [p1,p2,p3,p4]
					end
				end
			end
		}
		
	}
}

Oj.to_file("overlap_set.oj",overlap_set,:mode=>:compat)

overlap_cds=Hash.new
overlap_len=Hash.new
overlap_set.each{|scf,d1|
	overlap_cds[scf]=Hash.new
	overlap_len[scf]=Hash.new
	d1.each{|reset,d2|
		overlap_len[scf][reset]=Array.new
		d2.each{|set|
			
			if set.map{|x| x[0]}==["rs", "ce"] then
				p set.map{|x| x[0]}
				p set[1][1]-set[0][1]
				p cdreg_hash[scf][set[1]][1]
				p rereg_hash[scf][reset][set[0]][1]
				overlap_len[scf][reset] << [set.map{|x| x[0]},set[1][1]-set[0][1],cdreg_hash[scf][set[1]][1],rereg_hash[scf][reset][set[0]][1]]
			elsif set.map{|x| x[0]}==["cs", "re"] then
				p set.map{|x| x[0]}
				p set[1][1]-set[0][1]
				p cdreg_hash[scf][set[0]][1]
				overlap_len[scf][reset] << [set.map{|x| x[0]},set[1][1]-set[0][1],cdreg_hash[scf][set[0]][1],rereg_hash[scf][reset][set[1]][1]]
			elsif set.map{|x| x[0]}==["cs","rs","re","ce"] then
				p set.map{|x| x[0]}
				p set[2][1]-set[1][1]
				p cdreg_hash[scf][set[0]][1]
				overlap_len[scf][reset] << [set.map{|x| x[0]},set[2][1]-set[1][1],cdreg_hash[scf][set[0]][1],rereg_hash[scf][reset][set[1]][1]]
			elsif set.map{|x| x[0]}==["rs","cs","ce","re"] then
				p set.map{|x| x[0]}
				p set[2][1]-set[1][1]
				p cdreg_hash[scf][set[1]][1]
				overlap_len[scf][reset] << [set.map{|x| x[0]},set[2][1]-set[1][1],cdreg_hash[scf][set[1]][1],rereg_hash[scf][reset][set[0]][1]]
			end

			# set.each{|pos|
			# 	if cdreg_hash[scf][pos] then
			# 		ovarr << cdreg_hash[scf][pos][0]
			# 	end
			# }
			
		}
		
	}
}

Oj.to_file("overlap_len_all.oj",overlap_len,:mode=>:compat)
# allscf_overlap=Hash.new
# overlap_cds.each{|scf,d1|
# 	d1.each{|reset,d2|
# 		unless allscf_overlap[reset] then
# 			allscf_overlap[reset]=Array.new
# 		end
# 		# p reset
# 		# p cdidx_scf[scf].length
# 		# p d2.length
# 		allscf_overlap[reset] <<  [cdidx_scf[scf].length,d2.length/cdidx_scf[scf].length.to_f]
# 	}
# }
# Oj.to_file("allscf_overlap.oj",allscf_overlap,:mode=>:compat)