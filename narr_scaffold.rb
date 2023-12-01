require 'csv'
require 'bio'
require 'oj'



b_no=Hash.new
CSV.foreach("base_number.csv") do |row|
	b_no[row[0]]=row[1].to_i
end
p b_no
ncomp=Hash.new
CSV.foreach("complement_base.csv") do |row|
	ncomp[row[0].to_i]=row[1].to_i
end

cpp=Hash.new
CSV.foreach("complement_pair_point.csv") do |row|
	cpp[[row[0],row[1]]]=row[2].to_i
end
ncpp=Hash.new
cpp.each{|set,pt|
	ncpp[[b_no[set[0]],b_no[set[1]]]]=pt
}

scffs=Bio::FlatFile.auto("apo.final.assembly.fasta")

tgt_scff=Hash.new
templ=Hash.new
scffs.each_entry do |f|
	scno= f.definition.sub("scaffold_","").to_i
	p scno
	templ[scno]=f.naseq.split("").map{|x| b_no[x]}
end
# p templ
Oj.to_file("ntemplate.oj",templ,:mode=>:compat)

=begin
i=0
renz=Hash.new
CSV.foreach("ddRAD_EnzymeSelect.csv") do |row|
	if i>0 then
		renz[row[0].to_i]=Hash.new
		renz[row[0].to_i]["name"]=row[1]
		renz[row[0].to_i]["recgsite"]=row[2]
		renz[row[0].to_i]["rsno"]=row[2].split("").map{|x| b_no[x]}
		renz[row[0].to_i]["cutpos"]=row[3].to_i
	end
	i+=1
end
# p renz

rsites=Hash.new
templ.each{|k,v|
	p k
	rsites[k]=Hash.new
	renz.each{|rno,d1|

		rname=d1["name"]
		p rname
		rsites[k][rname]=Hash.new
		rnseq= d1["rsno"]
		rrlen=d1["rsno"].length
		sarr=Array.new
		tgt_end=v.length - rrlen+1
		(0..tgt_end).to_a.each_with_index{|tpos,idx|
			if v[tpos,rrlen]==rnseq then
				# p "#{idx+1} #{v[tpos,rrlen]}"
				sarr << idx+1+d1["cutpos"]
			end
		}
		rsites[k][rname]=sarr
	}
}
# rsites.each{|k,d1|
# 	d1.each{|rname,d2|
# 		p "#{rname}  #{d2.length}  #{d2}"
# 	}
# }
rrsites=Hash.new
rsites.each{|k,d1|
	rrsites[k]=Hash.new
	d1.each{|rname1,d2_1|
		d1.each{|rname2,d2_2|
			unless rname1==rname2 then
				# p "#{rname1} #{rname2}"
				rrsites[k][[rname1,rname2]]= (d2_1 | d2_2).sort
			end
		}
	}
}

tgtfragment=Hash.new
rrsites.each{|k,rcsites|
	tgtfragment[k]=Hash.new
	
	# p k
	rcsites.each{|reset,csites|
		tfrag=Array.new
		csites.each_with_index{|pt,idx|
	 		# p "#{pt}  #{idx-1}"
		 	if idx>0 then
		 	# 	p pt
				# p csites[idx-1]
		 		flen=pt - csites[idx-1]
		 		if flen > 99 && flen < 201 then
		 			tfrag << [csites[idx-1],pt]
		 		end
		 	end
		}
		tgtfragment[k][reset]=tfrag
	}

	# d.each_with_index{|pt,idx|

	# }
}
Oj.to_file("tgtfragment.oj",tgtfragment,:mode=>:compat)
=end

# tgtfragment.each{|k,d1|
# 	d1.each{|reset,d2|
# 		p "#{reset} #{d2.length}"
# 	}
# }