module myBioXAM
using Reexport
@reexport using XAM, GenomicFeatures

using my
mypth=joinpath(homedir(), "myjulia")
include(joinpath(mypth, "compatible.jl"))

#{{ bamReadCount
#[[ bamReadCount ]]
# using GenomicFeatures, XAM #or using Bio.Align #Required.
# readnum, site_len = bamReadCount(indexed_bam_file, (SiteChr, SitePos[, SiteStrand]); revstrand=false, UMI=nothing, UMI2=nothing, UMI_UMI2_linker=" ", rm_overlap_sites=true, read_span_site=false, readfilter= rd->rd["NH"]::UInt8==0x01)
# readnum, site_len, site_grp = bamReadCount(...; sitegrp=Group, ...)
# readnum, site_len, site_index::Int, smp_barcode|smpBCfun_output = bamReadCount(...; smpBC="BC_field", smpBCfun=..., ...)
# readnum, site_len, site_grp, smp_barcode = bamReadCount(...; smpBC="BC_field",  sitegrp=Group, ...)
# Counting overlapped unique-mapped read number in each site or site group. A read overlapped with more than one site in a site group will only be counted once.
# The read is counted by "touch" model, i.e., a read which only parly overlapped with a site will also be counted. Only read with tag "NH:i:1" can be counted, which aims to filter uniquely mapped reads.
# `UMI` is the field of UMI barcode. e.g., for 10x data, set UMI="UB". When UMI is setted, the outputted readnum is actually the unique UMI number for the reads touched a site or a site group. Otherwise, the outputted readnum is the number of unique read ID.
# `UMI2` It is allowed to set two UMI fields. The two UMI fields will be concatenated by `UMI_UMI2_linker`.
# `smpBC` is the sample barcode. e.g. for 10x data, set smpBC="CB" to get the gene read number of every cell x every gene.
# Tip: for 10x 5'-scRNAseq data, set UMI="UB", smpBC="CB", revstrand=true .
# `smpBCfun` custom function to modify sample barcode. Its output should either be a String or nothing (ignore this read). 
# When both `smpBCfun` and `UMI` are assigned, the default UMI2 value is smpBC in order to make sure that the result read number is restrictly equal to the sum of all cells in each type in case the sumBCfun merged barcodes (e.g.  converting BC to cell type in 10x data). In case smpBCfun is an one-to-one mapping function (e.g. UMI="UB", smpBCfun=x->replace(x, r"\-1$"=>"")), setting UMI2=nothing can save memery without accuracy compromise.
# When sitegrp is given, the read will be counted by group. Outputted readnum is the number of unique read ID or UMI for all the reads touched any sites of a group.
# When rm_overlap_sites=true, the site segments overlapped by different site groups will be removed in advance. These parts will also not be counted for `site_len`.
# when read_span_site=true, only the read fully span a site will be counted. e.g., the site could be splicing site +/- anchor length in order to count non-junction reads in this mode.
# readfilter=... assigns a function to filter reads. The input is read object in XAM.jl. The default one is used to filter out multiple mapped reads in STAR bam file. It need to be changed for other aligner. e.g., for Bowtie2 output, set:
# readfilter= rd -> !haskey(rd, "XS") || rd["AS"]::Union{UInt8, UInt16} > rd["XS"]::Union{UInt8, UInt16}
# for BWA output, set:
# readfilter=rd->(BAM.flags(rd) & 0x400)==0 && !haskey(rd, "XA") && !haskey(rd, "SA")
# BWA will assign an arbitrary location on unmapped reads by BWA. See https://www.biostars.org/p/141014/ . But it won't be a problem of bamReadCount() since function will check BAM.ismapped() anyhow. 10 Jun 2021
# See also: bamReadLoop, bamJuncCount, bamJuncCount_10x
# Xiong Jieyi, 27 Sep 2019 > 27 Feb 2020 >10 Jun 2021>12 Jan 2024

export bamReadCount
function bamReadCount(bamfn::AbstractString, Site::Tuple;
                      revstrand::Bool=false, rm_overlap_sites::Bool=true,
                      readfilter::Function=rd->rd["NH"]::UInt8==0x01,
                      UMI::Union{String, Nothing}=nothing,
                      smpBC::Union{String, Nothing}=nothing, smpBCfun=nothing,
                      UMI2::Union{String, Nothing}=!isnothing(smpBCfun) && !isempty(UMI) ? smpBC : nothing,
                      UMI_UMI2_linker::String=" ",
                      sitegrp::Union{Group, Nothing}=nothing, read_span_site::Bool=false,
                      is10x::Bool=false)
    if is10x #Abolished in 27 Mar 2020
        @error("is10x=true has abolished. Use UMI=\"UB\", smpBC=\"CB\" instead.")
        # UMI="UB"
        # smpBC="CB"
    end
    # BAM=importpkg(:XAM, preloaded=true).BAM::Module
    # eachoverlap=importpkg(:GenomicFeatures, preloaded=true).eachoverlap::Function
    isBC=!isnothing(smpBC)
    isSmpBCfun=!isnothing(smpBCfun)
    isstrand=length(Site)==3
    if isstrand
        Site::Tuple{AbstractVector{String}, AbstractMatrix{Int}, AbstractVector{Char}}
    else
        Site::Tuple{AbstractVector{String}, AbstractMatrix{Int}}
    end
    
    # Segment proprecess
    if isnothing(sitegrp)
        isSiteGrp=false
        sitegrp=1:rnum(Site)
        ositeidx=1:rnum(Site)
    else
        isSiteGrp=true
        Site, sitegrp=grpfunwith(sitegrp, Site) do csite
            pos=segunion(csite[2])[1]
            pos_len=size(pos, 1)
            if isstrand
                (fill(csite[1][1], pos_len), pos, fill(csite[3][1], pos_len))
            else
                (fill(csite[1][1], pos_len), pos)
            end
        end
    end
    if rm_overlap_sites
        segSt, t=segconc(Site)
        dupSt=getrow(segSt, t.>1)
        Site, t=segmentdiff(Site, dupSt, BnoOverlap=true)
        sitegrp=getrow(sitegrp, t)
    end

    if !isfile(bamfn*".bai")
        println("Bam file $bamfn haven't be indexed yet. Indexing....")
        sf"samtools index $1"(bamfn)
    end
    bamO=open(BAM.Reader, bamfn, index=bamfn*".bai")
    l=ismbr(Site[1], bamO.refseqnames)
    if !any(l)
        println("Site ref like: "*Site[1][1])
        println("bam file ref like: "*bamO.refseqnames[1])
        error("Site reference name does not match with bam file.")
    end
    Site=getrow(Site, l)
    sitegrp=getrow(sitegrp, l)
    siteLenTb=tb(c"sitelen, sitegrp"=>grpfun(x->sum(seglen(x)), sitegrp, Site[2]))
    println("Site proprecessing has done.")
    
    siteGrpNum=rownum(siteLenTb)

    isUMI=!isnothing(UMI)
    isUMI2=!isnothing(UMI2)
    if !isUMI && isUMI2
        error("UMI2 can be assigned only when UMI is assigned.")
    end
    
    if isBC
        bcpool=Dict{String, Set{String}}()
        grpfun1=grpfunwith
    else
        rdpool=Set{String}()
        grpfun1=grpfun
    end
    C, stgrp=grpfunwith(sitegrp, Site, no=true, id=true) do sitegrpno::Int, sitegrpid, csite::Tuple
        @timetip println(f"$1 out of $2 ($(3=.3g)%) site groups are done."(sitegrpno, siteGrpNum, 100*sitegrpno/siteGrpNum))
        empty!(isBC ? bcpool : rdpool)
        cstrand=isstrand ? csite[3][1] : '.'
        in(cstrand, Set(['+', '-', '.'])) || error("Invalid strand '$cstrand'.")
        chr=csite[1][1]::String
        for pos in eachrow(csite[2])
            for rd in eachoverlap(bamO, chr, max(pos[1], 1):pos[2]) # eachoverlap() will throw error for position<1.
                # if cstrand=='.' || xor(BAM.ispositivestrand(rd)==(cstrand=='+'), revstrand) && rd["NH"]::UInt8==0x01 && (!isUMI || haskey(rd, UMI)) #Version before 23 Sep 2020
                #Version after 23 Sep 2020 >>>
                # rd["NH"]::UInt8==0x01 || continue
                BAM.ismapped(rd) || continue #This filter is necessary for BWA mapping, because unmapped reads could stil be here as an arbitrary (wrong) location could be assigned on unmapped reads by BWA. See https://www.biostars.org/p/141014/ (10 Jun 2021)
                (!isBC || haskey(rd, smpBC)) || continue
                (!isUMI || haskey(rd, UMI)) || continue
                (!isUMI2 || haskey(rd, UMI2)) || continue
                readfilter(rd) || continue
                
                if cstrand!='.'
                    isPEread2=(BAM.flags(rd) & 0x81) == 0x81
                    readpositive=xor(BAM.ispositivestrand(rd), isPEread2)
                    # @assert !xor(rd["XS"]=='+', readpositive) #Only for debug. Should be removed. 23 Sep 2020
                    if xor(readpositive, xor((cstrand=='+'), revstrand))
                        continue
                    end
                end
                #<<<
                           
                if read_span_site
                    #Check if read fragment fully spanned site.
                    touched=false
                    tagpos=parseCIGARpos(BAM.cigar(rd), BAM.position(rd))
                    for i=1:size(tagpos, 1)
                        if tagpos[i, 1]<=pos[1] && pos[2]<=tagpos[i, 2]
                            touched=true
                            break
                        elseif tagpos[i, 1]>pos[2]
                            break
                        end
                    end
                else
                    #Check if read intron fully spanned site.
                    touched=true
                    if in('N', (cigar=BAM.cigar(rd);))
                        tagpos=parseCIGARpos(cigar, BAM.position(rd))
                        for i=1:size(tagpos, 1)-1
                            if tagpos[i, 2]<pos[1] && pos[2]<tagpos[i+1, 1]
                                touched=false
                                break
                            elseif tagpos[i+1, 2]>=pos[1]
                                break
                            end
                        end
                    end
                end
                if touched
                    rdid=if isUMI
                        if isUMI2
                            (rd[UMI]::String) * UMI_UMI2_linker * (rd[UMI2]::String)
                        else
                            rd[UMI]::String
                        end
                    else
                        BAM.tempname(rd)::String
                    end
                    if isBC
                        # if haskey(rd, smpBC)
                        cbc0=rd[smpBC]::String
                        cbc=isSmpBCfun ? smpBCfun(cbc0)::Union{String, Nothing} : cbc0
                        if !isnothing(cbc)
                            if haskey(bcpool, cbc)
                                push!(bcpool[cbc], rdid)
                            else
                                bcpool[cbc]=Set{String}([rdid])
                            end
                        end
                        # end
                    else
                        push!(rdpool, rdid)
                    end
                end
            end
        end
        if isBC
            (collect(keys(bcpool)), map(length, values(bcpool)))
        else
            length(rdpool)
        end
    end
    close(bamO)
    println("All finished.")
    if isBC
        smp_bc, rdnum=C
    else
        rdnum=C
    end
    if isBC
        rdnum, dtshift(stgrp, siteLenTb["sitegrp"], siteLenTb["sitelen"]), stgrp, smp_bc
    elseif isSiteGrp
        rdnum, dtshift(stgrp, siteLenTb["sitegrp"], siteLenTb["sitelen"]), stgrp
    else
        dtshift(ositeidx, stgrp, rdnum, 0), dtshift(ositeidx, siteLenTb["sitegrp"], siteLenTb["sitelen"], 0)
    end
end
#}}

#{{ bamReadLoop
#[[ bamReadLoop ]]
# using GenomicFeatures, XAM #or using Bio.Align #Required.
# bamReadLoop(bamfile, indexed_bam_file, (SiteChr, SitePos[, SiteStrand]); revstrand=false, ignore_jump_over=false, read_span_site=false, readfilter= rd->rd["NH"]::UInt8==0x01) do read, site_no, read_no_in_site, is_last_read_in_site, parsedCIGARpos|nothing
# # do-parameters may copy this: rd::BAM.Record, si::Int, rdi::Int, islast::Bool, _
#     ... do your own function ...
#     return nothing|:next_site|:exit
# end
# A framework to handle BAM file reads in specific region.
# read_span_site: When it is true, only the read fully span a site will be counted. e.g., the site could be splicing site +/- anchor length in order to count non-junction reads in this mode.
# ignore_jump_over: Whether to skip(T) the reads if its junctions ('N' in CIGAR) acrosses over the whole site without any touch.
# readfilter=... Assigns a function to filter reads. The input is read object in XAM.jl. The default one is used to filter out multiple mapped reads in STAR bam file. It need to be changed for other aligner. e.g., for Bowtie2 output, set:
#   readfilter= rd -> !haskey(rd, "XS") || rd["AS"]::Union{UInt8, UInt16} > rd["XS"]::Union{UInt8, UInt16}
# Custom function have five inputs:
#   read: The XAM.Record object;
#   site_no: Row number of matched site. The process order is restrict ascending.
#   read_no_in_site: The n-th passed reads in this site;
#   is_last_read_in_site: Whether it is the last read in the site. Can be used to summary the result of the whole site.
#   parsedCIGARpos|nothing: When ignore_jump_over=true, read CIGAR will be parsed. To avoid recalculation, this parameter is the parseCIGARpos(BAM.cigar(rd), BAM.position(rd)) result. When ignore_jump_over=false, this parameter is always nothing.
# Cunstom function can only return nothing, :next_site or :exit
#   :next_site -- ignore other reads in this site and jump to next site.
#   :exit      -- terminate bamReadLoop immediately.
# See also: bamReadCount, bamJuncCount, bamJuncCount_10x
# 13 Apr 2022

export bamReadLoop
function bamReadLoop(dofun::Function, bamfn::AbstractString, Site::Union{Tuple{AbstractVector{<:AbstractString}, AbstractMatrix{<:Integer}}, Tuple{AbstractVector{<:AbstractString}, AbstractMatrix{<:Integer}, AbstractVector{Char}}};
                     revstrand::Bool=false,
                     readfilter::Function=rd->rd["NH"]::UInt8==0x01,
                     ignore_jump_over::Bool=false,
                     read_span_site::Bool=false)
    # BAM=importpkg(:XAM, preloaded=true).BAM::Module
    # eachoverlap=importpkg(:GenomicFeatures, preloaded=true).eachoverlap::Function
    isstrand=length(Site)==3
    
    function isreadpassed(rd, pos::Matrix{Int}, cstrand::Char)
        BAM.ismapped(rd) || return (false, nothing)
        readfilter(rd) || return (false, nothing)
        if cstrand!='.'
            isPEread2=(BAM.flags(rd) & 0x81) == 0x81
            readpositive=xor(BAM.ispositivestrand(rd), isPEread2)
            if xor(readpositive, xor((cstrand=='+'), revstrand))
                return (false, nothing)
            end
        end
        if ignore_jump_over
            tagpos=parseCIGARpos(BAM.cigar(rd), BAM.position(rd))
            if read_span_site
                #Check if read fragment fully spanned site.
                touched=false
                for i=1:size(tagpos, 1)
                    if tagpos[i, 1]<=pos[1] && pos[2]<=tagpos[i, 2]
                        touched=true
                        break
                    elseif tagpos[i, 1]>pos[2]
                        break
                    end
                end
            else
                #Check if read intron fully spanned site.
                touched=true
                # if size(tagpos, 1)>1 #not needed
                for i=1:size(tagpos, 1)-1
                    if tagpos[i, 2]<pos[1] && pos[2]<tagpos[i+1, 1]
                        touched=false
                        break
                    elseif tagpos[i+1, 2]>=pos[1]
                        break
                    end
                end
                # end
            end
            return (touched, tagpos)
        else
            if read_span_site
                return (BAM.position(rd)<=pos[1] && BAM.rightposition(rd)>=pos[2], nothing)
            else
                return (true, nothing)
            end
        end
    end
    
    if !isfile(bamfn*".bai")
        println("Bam file $bamfn haven't be indexed yet. Indexing....")
        sf"samtools index $1"(bamfn)
    end
    open(BAM.Reader, bamfn, index=bamfn*".bai") do bamO
        for (ri, t) in enumerate(eachr(Site))
            chr, pos, strand = if isstrand
                t
            else
                (t..., '.')
            end
            oldrd=nothing
            oldtagpos=nothing
            dofunout=nothing
            rdno=0
            for rd in eachoverlap(bamO, chr, max(pos[1], 1):pos[2]) # eachoverlap() will throw error for position<1.
                ispass::Bool, tagpos::Union{Nothing, Matrix{Int}}=isreadpassed(rd, pos, strand)
                if ispass
                    if rdno>0
                        dofunout=dofun(oldrd, ri, rdno, false, oldtagpos)::Union{Symbol, Nothing}
                    end
                    oldrd=rd
                    oldtagpos=tagpos
                    rdno+=1
                    if dofunout==:next_site || dofunout==:exit
                        break
                    elseif !isnothing(dofunout)
                        error("Invalid custom function output. It can only be nothing, :next_site or :exit.")
                    end
                end
            end
            if dofunout==:next_site
                continue
            elseif dofunout==:exit
                break
            end
            if rdno>0 #For the last read in a site.
                dofunout=dofun(oldrd, ri, rdno, true, oldtagpos)::Union{Symbol, Nothing}
                if dofunout==:next_site
                    #do nothing
                elseif dofunout==:exit
                    break
                elseif !isnothing(dofunout)
                    error("Invalid custom function output. It can only be nothing, :next_site or :exit.")
                end
            end
        end
    end
    return nothing
end
#}}

#{{ bamJuncCount bamJuncCount_10x
#[[ bamJuncCount ]]
# using XAM
# juncT = bamJuncCount(bamfile; getstrand=:none(default)|:same|:rev, outpos=:jc(default)|:itn,
#                      readfilter::Function=rd->rd["NH"]::UInt8==0x01,
#                      region=( "chr1", [10000 20000] ) 
# Extract the junction information in a normal bam file, and count the read for each junction.
# Output is a tabel with 3-4 fields: chrno::Int8, itn_pos|jc_pos[, strand], readnum. strand field will miss when strand=:none.
# When outpos=:jc, function output the jc_pos field, where the positions are the boundary of exon. When outpos=:itn, function output the itn_pos field, where the positions are the boundary of intron, i.e., exon boundarys +/- 1.
# readfilter=... assigns a function to filter reads. The input is read object in XAM.jl. The default one is used to filter out multiple mapped reads in STAR bam file. It need to be changed for other aligner. e.g., for Bowtie2 output, set:
# readfilter= rd -> !haskey(rd, "XS") || rd["AS"]::Union{UInt8, UInt16} > rd["XS"]::Union{UInt8, UInt16}
# for BWA output, set:
# readfilter=rd->(BAM.flags(rd) & 0x400)==0 && !haskey(rd, "XA") && !haskey(rd, "SA")
# See also: bamJuncCount_10x, bamReadCount, bamReadLoop
# 19 Sep 2019 > 27 Sep 2019 > 20 Feb 2020 > 6 Apr 2020 > 26 Jul 2021

export bamJuncCount
function bamJuncCount(bamfile::AbstractString; getstrand::Symbol=:none, outpos::Symbol=:jc, readfilter::Function=rd->rd["NH"]::UInt8==0x01, region::Union{Tuple{<:AbstractString, <:AbstractMatrix{<:Integer}}, Nothing}=nothing)

    in(getstrand, Set([:same, :rev, :none])) || error("Invalid parameter strand=$strand. Should only be :same, :rev or :none.")
    in(outpos, Set([:jc, :itn])) || error("Invalid pos_style=$pos_style. Should only be :jc or :itn .")
    
    # BioAli=importpkg(:XAM, preloaded=true)
    # BAM=BioAli.BAM
    
    # bam=open(BAM.Reader, bamfile)
    bamf, rds=if isnothing(region)
        t=open(BAM.Reader, bamfile)
        (t, t)
    else
        if !isfile(bamfile*".bai")
            println("Bam file $bamfn haven't be indexed yet. Indexing....")
            sf"samtools index $1"(bamfile)
        end
        tbam=open(BAM.Reader, bamfile, index=bamfile*".bai")
        @assert(in(region[1], tbam.refseqnames), f"Site reference name (like \"$1\") does not match with bam file (like \"$2\")."(region[1], tbam.refseqnames[1]))
        @assert(size(region[2])==(1, 2), "Invalid region format. It should be like region=(\"chr1\", [100 200]).")
        trds=eachoverlap(tbam, region[1], region[2]|>x->max(x[1], 1):x[2])
        (tbam, trds)
    end
    
    jcNoD=Dict{Tuple{Int8, Matrix{Int}, Char}, Int}()
    maxJcNo=1
    oldchrno=Int8(0)
    for li in rds
        BAM.ismapped(li) || continue
        readfilter(li) || continue
        cg=BAM.cigar(li)
        in('N', cg) || continue
        pos1=BAM.position(li)
        chrno=Int8(chr2no(BAM.refname(li)))
        chrno>Int8(0) || continue
        if chrno!=oldchrno
            oldchrno=chrno
            isnothing(region) && println("Chr_no: $chrno")
        end
        strand=if getstrand==:same
            BAM.ispositivestrand(li) ? '+' : '-'
        elseif getstrand==:rev
            BAM.ispositivestrand(li) ? '-' : '+'
        else
            '.'
        end
        rdpos=parseCIGARpos(cg, pos1)
        jcpos=[rdpos[1:end-1, 2] rdpos[2:end, 1]]
        for pos in eachr(jcpos)
            jcNoD[(chrno, pos, strand)]=get(jcNoD, (chrno, pos, strand), 0)+1
        end
    end
    close(bamf)
    Tj=if getstrand==:none
        if outpos==:itn
            mapr(jcNoD) do ((chrno, jcpos, strand), readnum)
                ds(chrno=chrno, itn_pos=jcpos|>x->[x[:, 1].+1 x[:, 2].-1], readnum=readnum)
            end
        else
            mapr(jcNoD) do ((chrno, jcpos, strand), readnum)
                ds(chrno=chrno, jc_pos=jcpos, readnum=readnum)
            end
        end
    else
        if outpos==:itn
            mapr(jcNoD) do ((chrno, jcpos, strand), readnum)
                ds(chrno=chrno, itn_pos=jcpos|>x->[x[:, 1].+1 x[:, 2].-1], strand=strand, readnum=readnum)
            end
        else
            mapr(jcNoD) do ((chrno, jcpos, strand), readnum)
                ds(chrno=chrno, jc_pos=jcpos, strand=strand, readnum=readnum)
            end
        end
    end
    return Tj
end

#[[ bamJuncCount_10x ]]
# using XAM
# juncT, readT = bamJuncCount_10x(bamfile; outpos=:jc(default)|:itn, mixstrand=false,
#                      region=( "chr1", [10000 20000] ) )
# Extract the junction information from a 10x bam file (mapped by cellRanger), and count the read for each junction.
# juncT is a table with four fields: chrno::Int8, itn_pos|jc_pos, strand, i_jcno=jcno
# readT is a table with three fields: i_jcno, cellBC, uni_readnum
# When outpos=:jc, juncT has jc_pos field, where the positions are the boundary of exon. When outpos=:itn, juncT has the itn_pos field, where the positions are the boundary of intron, i.e., exon boundarys +/- 1.
# When mixstrand=true, juncT has no 'strand' field, and the 'uni_readnum' is the sum reads in both strands.
# See also: bamJuncCount, bamReadCount
# 19 Sep 2019 > 27 Sep 2019 > 6 Apr 2020 > 22 Mar 2023

export bamJuncCount_10x
function bamJuncCount_10x(bamfile::AbstractString; outpos::Symbol=:jc, region::Union{Tuple{<:AbstractString, <:AbstractMatrix{<:Integer}}, Nothing}=nothing, mixstrand::Bool=false)

    in(outpos, Set([:jc, :itn])) || error("Invalid pos_style=$pos_style. It can only be :jc or :itn .")
    
    # BioAli=importpkg(:XAM, preloaded=true)
    # BAM=BioAli.BAM
    
    bamf, rds=if isnothing(region)
        t=open(BAM.Reader, bamfile)
        (t, t)
    else
        if !isfile(bamfile*".bai")
            println("Bam file $bamfn haven't be indexed yet. Indexing....")
            sf"samtools index $1"(bamfile)
        end
        tbam=open(BAM.Reader, bamfile, index=bamfile*".bai")
        @assert(in(region[1], tbam.refseqnames), f"Site reference name (like \"$1\") does not match with bam file (like \"$2\")."(region[1], tbam.refseqnames[1]))
        @assert(size(region[2])==(1, 2), "Invalid region format. It should be like region=(\"chr1\", [100 200]).")
        trds=eachoverlap(tbam, region[1], region[2]|>x->max(x[1], 1):x[2])
        (tbam, trds)
    end
    
    JCuni=Dict{Tuple{Int, String}, Set{String}}()
    jcNoD=Dict{Tuple{Int8, Matrix{Int}, Char}, Int}()
    maxJcNo=1
    oldchrno=Int8(0)
    for li in rds
        BAM.ismapped(li) || continue
        cg=BAM.cigar(li)
        in('N', cg) || continue
        pos1=BAM.position(li)
        chrno=Int8(chr2no(BAM.refname(li)))
        chrno>Int8(0) || continue
        haskey(li, "CB") || continue
        haskey(li, "UB") || continue
        UMI=li["UB"]::String
        cellBC=split(li["CB"]::String, '-')[1]
        if chrno!=oldchrno
            oldchrno=chrno
            if isnothing(region)
                println(no2chr(chrno))
            end
        end
        strand=mixstrand ? '.' : BAM.ispositivestrand(li) ? '-' : '+' #The read direction is reversed.
        li["NH"]::UInt8==0x01 || continue #Filter unique mapped reads
        rdpos=parseCIGARpos(cg, pos1)
        jcpos=[rdpos[1:end-1, 2] rdpos[2:end, 1]]
        for pos in eachr(jcpos)
            jcno=get(jcNoD, (chrno, pos, strand), nothing)
            if isnothing(jcno)
                jcno=maxJcNo
                jcNoD[(chrno, pos, strand)]=maxJcNo
                maxJcNo+=1
            end
            hs=(jcno, cellBC)
            if isnothing((t=get(JCuni, hs, nothing);))
                JCuni[hs]=Set{String}([UMI])
            else
                push!(t, UMI)
            end
        end
    end
    close(bamf)

    Tc=mapr(JCuni) do ((jcno, cellBC), cset)
        ds(i_jcno=jcno, cellBC=cellBC, uni_readnum=length(cset))
    end
    Tj=if outpos==:itn
        mapr(jcNoD) do ((chrno, jcpos, strand), jcno)
            ds(chrno=chrno, itn_pos=jcpos|>x->[x[:, 1].+1 x[:, 2].-1], strand=strand, i_jcno=jcno)
        end
    else
        mapr(jcNoD) do ((chrno, jcpos, strand), jcno)
            ds(chrno=chrno, jc_pos=jcpos, strand=strand, i_jcno=jcno)
        end
    end
    Tj["uni_readnum"]=grpfun(sum, fastgrp(Tc["i_jcno"], Tj["i_jcno"]), Tc["uni_readnum"])
    if mixstrand
        delete!(Tj, "strand")
    end
    return (Tj, Tc)
end

#}}

addhelpfromfile(@__FILE__, inmodule=@__MODULE__)
end
