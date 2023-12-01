require 'bio'
require 'csv'
require 'oj'

scaffcds=Hash.new
CSV.foreach("Apo_augustus_ver2.gff", col_sep: "\t") do |row|
	sfid = row[0].sub("scaffold_","").to_i
	unless scaffcds[sfid] then
		scaffcds[sfid]=Array.new
	end
	if row[2]=="CDS" then
		scaffcds[sfid] << [row[3].to_i,row[4].to_i]
	end
end
p scaffcds

Oj.to_file("scaffcds.oj",scaffcds,:mode=>:compat)

locus=Hash.new
locus_no=0
locus2=Hash.new

minifiles=Dir.glob("Apo_augustus_ver2.gff")
# p minifiles
minifiles.each{|mf|
	File.open(mf) do |file|
		p mf
	  file.each_line do |l|
	    if l =~ /\A\# ----- prediction on sequence number/ then
	      al= l.chomp.sub("# ----- prediction on sequence number","").split(" ")
	      locus_no=al[0].to_i
	      locus[locus_no]=Array.new
	      locus[locus_no] << al[6].sub(")","")
	      locus[locus_no] << al[3].to_i
	      locus2[locus_no]=al[3].sub(",","").to_i
	    else
	      if locus_no > 0 then
	 	  	  locus[locus_no] << l.chomp
	      end
	    end 
	  end
	end
}

features=Hash.new
f_tag=0
gene_no=""
locus3=Hash.new

locus.each{|locus_no,data|
	if data.length > 3 then
		data.each{|l|
			# p l
			if l =~ /\A\# start gene g/ then
				gene_no= l.split(" ")[3]
				unless features[gene_no] then
					features[gene_no]=Array.new
				end
				# p  features
				f_tag=1
			elsif l =~ /\A\# end gene g/ then
				f_tag=0
			end
			if f_tag==1 then
				features[gene_no] << l.sub("# start gene","").split(" ")
				locus3[gene_no]=[locus_no,locus2[locus_no]]
			end
    # p features
		}
	end
}

# p features
gdata2=Hash.new
features.each{|gno,data|
	gdata2[gno]=Array.new
	prot_seq=[]
	data.each{|l|
		if l[0]=~/\A#/ then
			if l[1]=="protein" then
				prot_seq << "protein"
				prot_seq << l[4].sub("[","").sub("]","")
			else
				prot_seq << l[1].sub("[","").sub("]","")
			end
		elsif l[1] == "AUGUSTUS" then
					gdata2[gno] << l[2,6]
		end
		# p gno
	}
	gdata2[gno]<< prot_seq
}

# p gdata2

gdata3=Hash.new

gdata2.each{|gno,data|
	gdata3[gno]=Hash.new
	data.each{|l|
		if l then
			unless gdata3[gno][l[0]] then
				gdata3[gno][l[0]]=Array.new
			end
			  gdata3[gno][l[0]] << l[1]
			  gdata3[gno][l[0]] << l[2]  
		    gdata3[gno][l[0]] << l[4]
		    # p l
		end    
		# p l
	}
}

# p gdata3

gdata3.each{|gno,data|
  data.each{|k,l|
		if k == "gene" || k == "intron" || k == "CDS" then
			unless data.has_key?("start_codon") then
				gdata3[gno][k] << "no_start"
			end
			unless data.has_key?("stop_codon") then
				gdata3[gno][k] << "no_stop"
			end
		end
		
		if k == "gene" then 
			if l[2] == "-" then
				if l.include?("no_start")
					gdata3[gno][k] = [l[0],">#{l[1]}",l[2..-1]].flatten				
				end
				if l.include?("no_stop")
					gdata3[gno][k] = ["<#{l[0]}",l[1],l[2..-1]].flatten
				end	
				if l.include?("no_start") && l.include?("no_stop") then
					gdata3[gno][k] = ["<#{l[0]}",">#{l[1]}",l[2..-1]].flatten
				end
			end
		end
		if k == "gene" then
			if l[2] == "+" then
				if l.include?("no_stop")
					gdata3[gno][k] = [l[0],">#{l[1]}",l[2..-1]].flatten
				end
				if l.include?("no_start")
					gdata3[gno][k] = ["<#{l[0]}",l[1],l[2..-1]].flatten
				end	
				if l.include?("no_start") && l.include?("no_stop") then
					gdata3[gno][k] = ["<#{l[0]}",">#{l[1]}",l[2..-1]].flatten
				end
			end
		end

		if k == "intron" || k == "CDS" then
			aa=Array.new
			# p l
			# p "#{l.length}  #{l.length/3}"
			(1..l.length/3).to_a.each{|i|
				aa << l[(i-1)*3,2]
			}
			aa << l[l.length/3*3-1]
			aa << l[l.length/3*3..-1]
			 # p aa
			aa2=Array.new
			if aa[-2] == "-" then
				if aa[-1].include?("no_start") then
					if aa.length > 3 then
						aa2 = [[aa[0..-4],[aa[-3][0],">#{aa[-3][1]}"]],aa[-2,2]]
					else
						aa2 = [[aa[-3][0],">#{aa[-3][1]}"],aa[-2,2]]
					end
				end	
				if aa[-1].include?("no_stop") then
					if aa.length > 3 then
						aa2 = [[["<#{aa[0][0]}",aa[0][1]],aa[1..-3]],aa[-2,2]]
					else
						aa2 = [["<#{aa[-3][0]}",aa[0][1]],aa[-2,2]]
					end
				end	
				if aa[-1].include?("no_start") && aa[-1].include?("no_stop") then	
					if aa.length > 3 then
						aa2 = [["<#{aa[0][0]}",">#{aa[-3][1]}"],aa[-2,2]]
					else
						aa2 = [["<#{aa[-3][0]}",">#{aa[-3][1]}"],aa[-2,2]]
					end
				end
			end

		 if aa[-2] == "+" then
				if aa[-1].include?("no_stop") then
					if aa.length > 3 then
						aa2=[[aa[0..-4],[aa[-4][0],">#{aa[-3][1]}"]],aa[-2,2]]
					else
						aa2=[[aa[-3][0],">#{aa[-3][1]}"],aa[-2,2]]
					end
				end	
				if aa[-1].include?("no_start") then
					if aa.length > 3 then
						aa2=[[["<#{aa[0][0]}",aa[0][1]],aa[1..-3]],aa[-2,2]]
					else
						aa2=[["<#{aa[-3][0]}",aa[0][1]],aa[-2,2]]
					end
				end	
				if aa[-1].include?("no_start") && aa[-1].include?("no_stop") then	
					if aa.length > 3 then
						aa2=[["<#{aa[0][0]}",">#{aa[-3][1]}"],aa[-2,2]]
					else
						aa2=[["<#{aa[-3][0]}",">#{aa[-3][1]}"],aa[-2,2]]
					end
				end
			end	

			unless aa[-1].include?("no_start") || aa[-1].include?("no_stop") then
				if aa.length > 3 then
						aa2=aa[0..-3],aa[-2,2]
			  else
						aa2=aa[0],aa[-2,2]
				end	
			end
    gdata3[gno][k]= aa2
    # p gno
    # p aa2
    end
    aa3=Array.new
		if k == "intron" || k == "CDS" then
			# if aa3.length >
			if gdata3[gno][k][-1].include?("-") then
        aa3=["complement",[gdata3[gno][k][0..-2]],gdata3[gno][k][-1]]
      else
      	aa3=["",[gdata3[gno][k][0..-2]],gdata3[gno][k][-1]]
      end 
    gdata3[gno][k]= aa3 
    end 

    if k == "gene"
    	if gdata3[gno][k].include?("-") then
    		aa3=["complement",[gdata3[gno][k][0..1]],gdata3[gno][k][2..-1]]
    	else
    	 	aa3=["",[gdata3[gno][k][0..1]],gdata3[gno][k][2..-1]]
    	end	
    gdata3[gno][k]= aa3
    end
    # p k
    # p gno
    # p aa3[1]
  }
}

gdata3.each{|gno,data|
  data.each{|k,l|
   if k == "exon" || k == "CDS" then
    p gno
    p k
    p l
   end
 }
} 	

