#{{ genomemap coordinatemap

#[[ genomemap coordinatemap ]]
#(SiteIdx, TagIdx)=genomemap(Site, Tag; fun=(si,ti,sp,tp), touch=false*, revstrand=false,
#                     minoverlap=0~1, no_match_output=(Int[],Int[]), mrows=false)
#Or ... = genomemap(fun, ...; ...)
#(SiteIdx, TagIdx)=coordinatemap(SitePos, TagPos; SiteId=index, TagId=index,
#                     sitesorted=false, tagsorted=false, sortcheck=true,
#                     minoverlap=0~1, no_match_output=(Int[],Int[]), mrows=false)
#Or ... = coordinatemap(fun, ... ; ... )
# Assignment of `minoverlap=` implys `touch=true`.
#Site and Tag are two tuples with (Chromosome, pos[, strand]). `pos` is a n x 2 matrix. Chromosoem could be other group IDs.
#When `revstrand=ture`, function will reverse Site's strand in the first.
#If the Tags is sorted by (chr[, strand], pos) order, set tagsorted as true to save time.
#In non-touch mode, any tag fall in site completely will be reported. In touch mode, the overlapped_region / min(site_length, tag_length) >= minoverlap will be reported.
#The input of `fun` is (SiteIdx::Int, TagIdx::Int, SitePos::Vector, TagPos::Vector). The function output will be pipeup by addrow!()(mrows=false) or addrows!()(mrows=true), except it is nothing or empty tuple ().
#In the case of no any match, function return `no_match_output`.
#`coordinatemap`: Mapping 1D coordinates without consider chromosome and strand.
#See also: coordinatemap, bamReadCount
#Xiong Jieyi, 8 Jul 2014 >December 8, 2014>February 8, 2015>Dec 22, 2015>Jan 4, 2016>8 Mar 2017>7 Apr 2022

# using my
export genomemap, coordinatemap
genomemap(Site::Tuple{Array,Matrix},Tag::Tuple{Array,Matrix}; kw...)=
    genomemap(tuple(Site...,Char[]),tuple(Tag...,Char[]); kw...)
function genomemap(Site::Tuple{Array{Tchr,N1},Matrix{Tpos},AbstractArray{Tstrand,N2}},
                   Tag::Tuple{Array{Tchr,N1},Matrix{Tpos},AbstractArray{Tstrand,N2}};
                   revstrand::Bool=false,no_match_output=(Int[],Int[]), kw...) where {Tchr,Tpos<:Real,Tstrand<:Char,N1,N2}
    
    @assert(length(Site)==length(Tag),"Input should both have or haven't strand information.")
    if isempty(Site[3]) && isempty(Tag[3])
        Site=Site[1:2]
        Tag=Tag[1:2]
    elseif revstrand
        Site=(Site[1],Site[2],dtshift(Site[3],['+','-'],['-','+']))
    end

    #Remove NaN
    SiteNo=collect(1:size(Site[2],1))
    TagNo=collect(1:size(Tag[2],1))

    l=.!isnan.(Site[2][:,1])
    Site=getrow(Site,l)
    SiteNo=SiteNo[l]

    l=.!isnan.(Tag[2][:,1])
    Tag=getrow(Tag,l)
    TagNo=TagNo[l]
    
    if length(Site)==2 #no strand
        tfg=fastgrp(Tag[1])
        sfg=fastgrp(Site[1],tfg.grpid)
    else
        tfg=fastgrp((Tag[1],Tag[3]))
        sfg=fastgrp((Site[1],Site[3]),tfg.grpid)
    end
    rlt=rowpile()
    for gi=1:sfg.grpnum
        ls=want(sfg,gi)
        if isempty(ls)#Add in 28 Aug 2014
            continue
        end
        lt=want(tfg,gi)
        crlt=coordinatemap(Site[2][ls,:], Tag[2][lt,:], SiteId=SiteNo[ls], TagId=TagNo[lt];
                           no_match_output=nothing, kw...)
        if crlt!=nothing
            addrows!(rlt, crlt)
        end
    end

    isvirgin(rlt) ? no_match_output : value(rlt)
end
genomemap(Site::Tuple{Array,Matrix,AbstractArray},Tag::Tuple{Array,Matrix,AbstractArray};args...)=
    genomemap(promote_tuple(Site,Tag)...;args...)
genomemap(fun::Function,args...; kw...)=genomemap(args...; kw..., fun=fun)

function coordinatemap(Site::Matrix{Ts},Tag::Matrix{Tt};
                       SiteId::Group=collect(1:size(Site,1)),TagId::Group=collect(1:size(Tag,1)),
                       sitesorted::Bool=false,tagsorted::Bool=false,sortcheck::Bool=true,
                       fun::Function=(si,ti,sp,tp)->(si,ti),
                       minoverlap::AbstractFloat=-2.0,touch::Bool=minoverlap>-1.0,
                       mrows::Bool=false, no_match_output=(Int[],Int[])) where {Ts<:Real,Tt<:Real}
    if sitesorted
        if sortcheck && !issortedr(Site)
            error("Site is not sorted.")
        end
    else
        (Site,t)=sortr(Site)
        SiteId=getrow(SiteId,t)
    end
    if tagsorted
        if sortcheck && !issortedr(Tag)
            error("Tag is not sorted.")
        end
    else
        (Tag,t)=sortr(Tag)
        TagId=getrow(TagId,t)
    end

    rlt=rowpile()
    ps=1
    pt=1
   if touch
        #Touch model
        while ps<=size(Site,1) && pt<=size(Tag,1)
            if Tag[pt,2]<Site[ps,1]
                pt+=1
            elseif Site[ps,2]<Tag[pt,1]
                ps+=1
            else
                oldpt=pt
                while true
                    if Site[ps,1]<=Tag[pt,2] #Hitted
                        sitepos=Site[ps,:]
                        tagpos=Tag[pt,:]
                        if (minoverlap<=0.0) || ( (min(sitepos[2],tagpos[2])-max(sitepos[1],tagpos[1])+1) / (min(sitepos[2]-sitepos[1],tagpos[2]-tagpos[1])+1) )>=minoverlap
                            crlt=fun(getrow(SiteId,ps),getrow(TagId,pt),sitepos,tagpos)
                            if !isnothing(crlt) && crlt!=()
                                if mrows
                                    addrows!(rlt,crlt)
                                else
                                    addrow!(rlt,crlt)
                                end
                            end
                        end
                    end
                    pt+=1
                    if pt>size(Tag,1) || Site[ps,2]<Tag[pt,1]
                        break
                    end
                end
                pt=oldpt
                ps+=1
            end
        end
    else
        #Match model
        while ps<=size(Site,1) && pt<=size(Tag,1)
            if Tag[pt,1]<Site[ps,1]
                pt+=1
            elseif Site[ps,2]<Tag[pt,1]
                ps+=1
            else
                oldpt=pt
                while true
                    if Tag[pt,2]<=Site[ps,2]#Hitted
                        crlt=fun(getrow(SiteId,ps),getrow(TagId,pt),Site[ps,:],Tag[pt,:])
                        if !isnothing(crlt) && crlt!=()
                            if mrows
                                addrows!(rlt,crlt)
                            else
                                addrow!(rlt,crlt)
                            end
                        end
                    end
                    pt+=1
                    if pt>size(Tag,1) || Site[ps,2]<Tag[pt,1]
                        break
                    end
                end
                pt=oldpt
                ps+=1
            end
        end
    end
    # out=value(rlt)
    # return out==nothing?no_match_output:out
    isvirgin(rlt) ? no_match_output : value(rlt)
 end
coordinatemap(fun::Function, args...;kw...)=coordinatemap(args...; kw..., fun=fun)

#}}

#{{ chr2no, no2chr, no2chr_celeg
#[[ chr2no no2chr no2chr_celeg ]]
# chrNO|chrNO = chr2no(chrID|chrID_vector)
# chrId = no2chr(chrNO; chr=TF|prefix="chr") | no2chr_celeg(chrNo; chr=TF|prefix="CHROMOSOME_") #C_elegans style: CHROMOSOME_III
#Convert chrID to chrNO or NO. NN=chrNN, 31=chrX, 32=chrY, 41=chrM( or 41=MT when chr=false or prefix is empty), 51=chr2a, 52=chr2b, 0=all others.
#See also: genomemap
#Xiong Jieyi, 13 Jul 2014>February 25, 2015>31 Aug 2019

export chr2no, no2chr, no2chr_celeg
function chr2no(chr::AbstractString)
    chr=replace(chr, r"^CHROMOSOME_|^Chr|^chr|^CHR" => "")
    Tb=Dict{AbstractString,Int}("I"=>1,"II"=>2,"III"=>3,"IV"=>4,"V"=>5,
        "X"=>31,"Y"=>32,
        "M"=>41,"MT"=>41,"Mt"=>41,"mt"=>41,"mtDNA"=>41,"MTDNA"=>41,"MtDNA"=>41,
        "2a"=>51,"2A"=>51,"2b"=>52,"2B"=>52)
    chrno=get(Tb,chr,0)
    if chrno==0
        try
            chrno=int(chr)
        catch
            chrno=0
        end    
    end
    return chrno
end
chr2no(chr::Vector{T}) where {T<:Union{AbstractString,Symbol}} = map(chr2no,chr)
chr2no(chr::Symbol)=chr2no(string(chr))
function no2chr(no::Real; chr::Bool=true, prefix::AbstractString=chr ? "chr" : "")
    Tb=Dict{Int,AbstractString}(31=>"X", 32=>"Y", 41=>isempty(prefix) ? "MT" : "M" ,
                                51=>"2a",52=>"2b")
    chr=get(Tb,no,"")
    if chr==""
        chr=string(int(no))
    end
    return prefix*chr
end
no2chr(chr::Vector{T};args...) where {T<:Real} = map(x->no2chr(x;args...),chr)
function no2chr_celeg(no::Real; chr::Bool=true, prefix::AbstractString=chr ? "CHROMOSOME_" : "")
    Tb=Dict{Int,AbstractString}(31=>"X",32=>"Y",41=>isempty(prefix) ? "MT" : "M",
                                1=>"I",2=>"II",3=>"III",4=>"IV",5=>"V")
    chr=get(Tb,no,"")
    if chr==""
        error("Unknown chrNO $no.")
    end
    return prefix*chr
end
no2chr_celeg(chr::Vector{T};args...) where {T<:Real} = map(x->no2chr_celeg(x;args...),chr)
#}}

#{{ seglen
#[[ seglen ]]
# lengths=seglen(pos)
#Equal to pos[:,2]-pos[:,1]+1, but seglen([0 0])=0.
#A error will be reported if any segment is less than 1nt.
#See also: chr2no, no2chr
#Xiong Jieyi, 11 Sep 2014>5 Oct 2014

#<v0.6# function seglen{T<:Real}(pos::Array{T})
function seglen(pos::Array{T}) where {T<:Real}
    l=.!isnan.(pos[:,1])
    all(pos[.!l,1].<=pos[.!l,2]) || error("Invalid position number.")
    len=pos[:,2]-pos[:,1]+1
    len[(pos[:,1].==0) .& (pos[:,2].==0)].=0
    return len
end
export seglen
#}}

#{{ readGeneralSamFile, readGeneralSamLines
#[[ readGeneralSamFile readGeneralSamLines ]]
# Table = readGeneralSamLines(iterator; pairend=false, strand=true, readidfun=ascii, chrfun=ascii, recordno_ref::Vector{Int}=[0], fields::Dict=ds(XX=default, ...), fieldfun::Dict=ds(XX=parsefun(string), ...), readfilter=fun(Dict))
# ... = readGeneralSamFile("alignment.sam" | "alignment.bam" | `pipe`; samtoolsparam="...",
#                          region=(chr|chrno|vec, [posstart posend])) or (chr|chrno|vec,), ...)
# Input:
#   iterator: Any iteratiable data, each item is a line of sam file (no newline in the end).
#   pairend: Paired reads?
#   strand:  Read has strand information?
#   readidfun & chrfun: function to handle readid and chr. Input is a SubString and output must be a row or nothing. if function output nothing, this record will be ignored.
#      e.g. chrfun=x->(chr2no(x)|>x->x==0 ? nothing : Int16(x)) will save reads as chrno.
#   recordno_ref: A 1-element Int vector to assign and access the start and end recordno.
#   fields: a one row table, keys are the needed optional fields in sam file, and values are the default values for the missing fields. If any default is undef, read missing this tag will be discarded.
#   fields data type: :i:-Int32, :f:-Float32, :Z:-String, :A:-Char
#   fieldfun: a dict, keys are the fields and values are parsing functions. Note that the input of this funciton are always a string. If any function output is undef, read missing this tag will be discarded.
#   readfilter: a callable object with a Bool output, to assign if the current read should be kept (true) or discarded (false). The input is a Dict with the keys of readid, chr, flag::Int16, seq::String, cigar::String, and all the required tags, where the readid, chr and tag values are processed by given functions.
#  region: only work in indexed bam file. region can be single one or mulitple regions. In the case of multiple regions, function will segunion them first to avoid duplicated output records (For the reasion, see samtools manual.). If input is chrno rather than chr-name, function will check the first record of bam file first, to determind which prefix ("chr", "CHR" or "") should be used for no2chr. If only chromosome is given, function will read all the reads in given chromosome(s).
#  samtoolsparam: a string for additional parameters of samtools view. Only work then input is a bam file.
#
# Output field:
#  record_no, readid, flag, chr, pos, juncno, read_len, ...
#
#See also: readStarSamFile
#Xiong Jieyi, 12 Apr 2018 > 16 May 2018 >30 May 2018 > 1 Apr 2020 > 14 Oct 2021

export readGeneralSamFile, readGeneralSamLines
function readGeneralSamLines(samlines;
                             pairend::Bool=false,
                             strand::Bool=true,
                             readidfun=ascii,
                             chrfun=ascii,
                             recordno_ref::Vector{Int}=[0],
                             fields::Dict=tb(),
                             fieldfun::Dict=tb(),
                             readfilter=nothing)
    #Bam file only save 32-bit numbers.
    tagcheckundef=String[]
    for (k, v) in fields
        if isa(v, Int64)
            fields[k]=Int32(v)
        elseif isa(v, Float64)
            fields[k]=Float32(v)
        elseif isa(v, UndefInitializer)
            push!(tagcheckundef, k)
        end
    end
    
    RECORD_NO=onlyone(recordno_ref)
    rlt=rowpile()
    rlt_fields=rowpile()
    
    for ctxt in samlines
        tc=split(chomp(ctxt),'\t')
        if length(tc)<3
            @warn("Ignore invalid line (record NO: $RECORD_NO): \"$ctxt\"")
            continue
        end
        FLAG=parse(Int16,tc[2])
        # if tc[3]=="*"
        if FLAG & 0x4!=0 # Unmapped reads
            continue
        end
        CHR=chrfun(tc[3])
        if isnothing(CHR)
            continue
        end
        
        READID=readidfun(tc[1])
        if isnothing(READID)
            continue
        end
        
        RECORD_NO+=1
        PAIRNO=pairend ? (FLAG&128>0)+1 : 0
        READLEN=Int32(length(tc[10])) #Add in 23 Sep 2014

        #Parse Appendences
        ADDFDS=deepcopy(fields)
        
        for ctc=tc[12:end]
            fdnm, fdty, fdval = split(ctc,':', limit=3)
            fdfun=get(fieldfun, fdnm, nothing)
            if isnothing(fdfun)
                fdfun=if fdty=="i"
                    x->parse(Int32, x)
                elseif fdty=="f"
                    x->parse(Float32, x)
                elseif fdty=="Z"
                    ascii
                elseif fdty=="A" #Char
                    x->x[1]
                else
                    error("Unknown field type: $fdty")
                end
            end
            if haskey(ADDFDS, fdnm)
                ADDFDS[fdnm]=fdfun(fdval)
            end
        end

        t=false
        for k in tagcheckundef
            if isa(ADDFDS[k], UndefInitializer)
                t=true
                break
            end
        end
        if t
            continue
        end
        
        if !isnothing(readfilter) && !readfilter(ds(ADDFDS, readid=READID, cigar=tc[6], chr=CHR, flag=FLAG, seq=tc[10]))
            continue
        end
        
        #Parse CIGAR
        STR=int(tc[4])
        POS=zeros(Int,1,2)
        LOC=zeros(Int,1,2)
        cigar=tc[6];p=1;rdp=1;rfp=STR;
        aL=findall(['A'<=x<='Z' for x in cigar])
        LOC_read_start=[];
        cjuncno=Int16(0);
        for cL=aL
            NO=int(cigar[p:cL-1]);p=cL+1;
            ccL=cigar[cL]
            if ccL=='S' || ccL=='H'
                rdp+=NO
                #S is not be calculated in the mapping start position.
                ## Remove this: ## rfp+=NO
                #Changed in Dec 14, 2015
            elseif ccL=='M'
                if !(POS[1]>0)
                    LOC[1]=rdp
		    if isempty(LOC_read_start)
		        LOC_read_start=rdp
		    end
                    POS[1]=rfp
                end
                rfp+=NO
                rdp+=NO
                LOC[2]=rdp-1
                POS[2]=rfp-1
            elseif ccL=='D'
                rfp+=NO;
            elseif ccL=='I'
                rdp+=NO
            elseif ccL=='N'
	        cjuncno+=1
                addrow!(rlt,RECORD_NO,READID,FLAG,CHR,POS,cjuncno,READLEN,PAIRNO)
                addrow!(rlt_fields, ADDFDS)
                # SEQ="" #For junction reads, only the first one assigned with sequence.
                # CIGAR4save=""
                rfp+=NO
                POS=zeros(1,2)
            else
                error("Unknown cigar char in read: $ccL")
            end
        end
        if cjuncno>0
            cjuncno+=1
        end
        addrow!(rlt,RECORD_NO,READID,FLAG,CHR,POS,cjuncno,READLEN,PAIRNO)
        addrow!(rlt_fields, ADDFDS)
        @timetip println(RECORD_NO)
    end
    rltvalue=value(rlt)
    if rltvalue==nothing
        return nothing
    end
    read=tb(value(rlt_fields), c"record_no,readid,flag,chr,pos,juncno,read_len,pairno"=>rltvalue)
    if !pairend
        delete!(read, "pairno")
    end
    if strand
        read["strand"]=map(x->x&16==0 ? '+' : '-',read["flag"])
        if pairend
            #Reverse second half-read in pair-end mode
            read["strand"][read["pairno"].==2]=map(x->x=='-' ? '+' : '-',read["strand"][read["pairno"].==2])
        end
    end
    recordno_ref[1]=RECORD_NO
    return read
end
function readGeneralSamFile(filename::Union{AbstractString, Cmd}, args...; region::Union{Tuple, Nothing}=nothing, recordno_ref::Vector{Int}=[0], samtoolsparam::AbstractString="", wargs...)
    if isa(filename, AbstractString) && splitext(filename)[2]==".bam"
        if isnothing(region)
            readGeneralSamLines(eachline(isempty(samtoolsparam) ? `samtools view $filename` : `samtools view $samtoolsparam $filename`), args...; wargs...)
        else
            #Check if bam file is indexed
            if !isfile(filename[1:end-1]*"i") && !isfile(filename*".bai")
                println("Indexing bam file $filename")
                run(`samtools index $filename`)
                println("Indexing done.")
            end
            if length(region)==1
                #Read all reads by chromosomes
                chrs=isa(region[1],Group) ? unival(region[1]) : [region[1]]
                if eltype(chrs)<:Real
                    #Check if chromosome in bam file is started with chr
                    S=readchomp(pipeline(`samtools view -H $filename`,`grep '^@SQ'`))
                    M=match(r"\bSN\:(chr|CHR|Chr|chromosome|CHROMOSOME|Chromosome)",S)
                    if M!=nothing
                        chrs=no2chr(chrs,prefix=M.captures[1])
                    else
                        chrs=no2chr(chrs,prefix="")
                    end
                end
                T,=rowfunwith(chrs) do chr
                    file=isempty(samtoolsparam) ? `samtools view $filename $chr` : `samtools view $samtoolsparam $filename $chr`
                    # println("Reading BAM file chromosome $chr ....")
                    readGeneralSamLines(eachline(file), args...; recordno_ref=recordno_ref, wargs...)
                end
                return T
            else
                #Read reads in some regions
                if isa(region[1],Group)
                    region,=segunion(region)
                else
                    region=([region[1]],region[2])
                end
                if eltype(region[1])<:Real
                    #Check if chromosome in bam file is started with chr
                    S=readchomp(pipeline(`samtools view -H $filename`,`grep '^@SQ'`))
                    M=match(r"\bSN\:(chr|CHR|Chr|chromosome|CHROMOSOME|Chromosome)",S)
                    if M!=nothing
                        region=(no2chr(region[1],prefix=M.captures[1]), Array{Int}(region[2]))
                    else
                        region=(no2chr(region[1],prefix=""), Array{Int}(region[2]))
                    end
                end
                T,=rowfunwith(region) do crg
                    rgstr="$(crg[1]):$(crg[2][1])-$(crg[2][2])"
                    file=isempty(samtoolsparam) ? `samtools view $filename $rgstr` : `samtools view $samtoolsparam $filename $rgstr`
                    # println("Reading BAM file region $rgstr ....")
                    readGeneralSamLines(eachline(file), args...; recordno_ref=recordno_ref, wargs...)
                end
                return T
            end
        end
    else # When input is not a bam file (sam file or pipeline):
        if isnothing(region) && isempty(samtoolsparam)
            readGeneralSamLines(eachline(filename), args...; wargs...)
        else
            error("`region=...' and `samtoolsparam=...' are only supported for bam file.")
        end
    end
end
#}}

#{{ readStarSamFile readStarSamLine
#[[ readStarSamFile readStarSamLine ]]
# read=readStarSamFile(filename|IOStream; savefile="", onlyuniq=false, pairend=false, strand=true, chrfun=int16(chr2no), keepHI=false, onlyjunc=false,
#                      [readprefixlen=auto, readiddm=auto] | readidfun=fun(Field1_Str),
#                      region=(chr|chrno|vec, [posstart posend])) or (chr|chrno|vec,),
#                      addfields=(c"name1,name2,...", [N1, N2, ...])
#Read Star/Tophat output sam/bam file (with fields AS, NH, nM(or NM) and HI(optional)). Note that this function will not check if the sam/bam file was aligned by Star/Tophat. For general purpose, readGeneralSamFile() is more recommended.
#
# read|nothing=readStarSamLine(String/String_Vector; recordno_ref=[0], pairend=false, strand=true, readprefixlen=auto | readidfun=fun(Field1_Str), chrfun=int16(chr2no), addfields=...)
#Read lines of sam file.
#
#Output field:
#  record_no, readid, flag, chrno, pos, juncno, read_len,
#  align_score(AS), mapped_loci_num(NH), alignment_index(HI, Only for STAR), mismatch_num(nM in STAR or NM in Tophat).
#
#In readStarSamFile, you can use addfields to add other fields in SAM file. for example, readStarSamFile(..., addfields=(c"sequence,CIGAR", [9,6])).
#
#When onlyuniq=true, function will remove low-score multiple mapping reads, and the multiple mapped reads with equal highest score will still be keeped. There is a "is_uniq_map" field showing the real unique mapped reads. This function only ultilized the align_score (AS field) and need to analyze the whole mapped reads. For this reason, in region mode or when onlyjunc=true, onlyuniq=true is not supported.
#In pair-end mode, paired reads just be regarded as two indipendent reads, with field pairno=1 or 2. If strand is needed, the strand of the second half pair will be reversed to make sure both half pairs have consistent strand information.
#If onlyjunc is ture, function will only keep junction reads (e.g. for alternative splicing analysis). Note that as the non-junction reads will be ignored before calculate uniquity, the duplicated reads selection (onlyuniq=true) will only considered among junction reads, and the is_uniq_map field will just be replaced by mapped_loci_num.==1.
#readprefixlen assigned the id of reads which will be removed. It will be detected automaticly if omitted. readiddm is the delimiter. If ignored, function will check if readid has ':', '.' or '_' (and error if none is found.)
#Read will be sorted as (chrno,strand,pos).
#chrfun: The function convert chromosome field (#3) to a number. If the number is zero, this record will be discarded.
#
#region only work in indexed bam file. region can be single one or mulitple regions. In the case of multiple regions, function will segunion them first to avoid duplicated output records (For the reasion, see samtools manual.). If input is chrno rather than chr-name, function will check the first record of bam file first, to determind which prefix ("chr", "CHR" or "") should be used for no2chr. If only chromosome is given, function will read all the reads in given chromosome(s).
#
#For readStarSamFile, recordno_ref should be a reference (one element vector) for record number counter. Use [0] for the start value. If it is missed, the read in all output will start with 1.
#
#When no any valid read was found, readStarSamFile throw an error, while readStarSamLine return nothing.
#
#See also: readGeneralSamFile, readGeneralSamLines, seglen, segconc, genomemap, parseGTF, parseDemulReport, parseBamMutation
#Xiong Jieyi, 13 Jul 2014>29 Jul 2014>25 Aug 2014>23 Sep 2014>January 12, 2015>May 17, 2015>Jun 12, 2015>Aug 4, 2015>Aug 21, 2015>Aug 23, 2015>Nov 23, 2015>Dec 10, 2015>Dec 14, 2015>Apr 19, 2016>Apr 25, 2016>May 17, 2016>Sep 9, 2016

function readStarSamFile(fp::IO;
                         pairend::Bool=false,
                         savefile::AbstractString="",
                         onlyuniq::Bool=false, #Changed the default value to false in 30 Jun 2021.
                         strand::Bool=true,
                         readprefixlen::Int=-1,readiddm='\0',
                         readidfun::Union{Function,Void}=nothing,
                         chrfun::Function=x->Int16(chr2no(x)),
                         keepHI::Bool=false,
                         addfields::Tuple{Vector,Vector}=([],[]),
                         onlyjunc::Bool=false )#Add in Sep 9, 2016

    if onlyjunc && onlyuniq
        error("When onlyjunc=true, onlyuniq=true is not supported.")
    end
    
    isaddfields=!isempty(addfields[1])
    rlt=rowpile()
    ctxt=chomp(readline(fp))
    if isempty(ctxt)
        @warn("File is empty.")
        return nothing
    end
    while ctxt[1]=='@'
        ctxt=chomp(readline(fp))
    end
    RECORD_NO=0
    while true
        tc=split(ctxt,'\t')
        if length(tc)<3
            @warn("Ignore invalid line (record NO: $RECORD_NO): \"$ctxt\"")
            if eof(fp)
                break
            end
            ctxt=chomp(readline(fp))
            continue
        end
        FLAG=parse(Int16,tc[2])
        if (FLAG & 0x4 != 0) || (CHRNO=chrfun(tc[3]))==0
            if eof(fp)
                break
            end
            ctxt=chomp(readline(fp))
            continue
        end

        #Parse READID
        READID=nothing
        try
            if isa(readidfun,Void)
                if readiddm=='\0'
                    for t in [':','.','_']
                        if in(t,tc[1])
                            readiddm=t
                            break
                        end
                    end
                end
                in(readiddm,tc[1])||error("readiddm is not exist in readid.")
                if readprefixlen<0
                    t=match(Regex("\\$readiddm[\\$readiddm\\d]+\$"),tc[1])
                    # t=match(r"\:[\:\d]+$",tc[1])
                    # if t==nothing
                    #     readprefixlen=findonly(collect(tc[1]).=='.')
                    #     readidfun=X->parse(Int32,X[readprefixlen+1:end])
                    # else
                    readprefixlen=length(tc[1])-length(t.match)+1
                    readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],readiddm))
                    # end
                else
                    @assert(tc[1]==0 || tc[1][readprefixlen]==readiddm || tc[1][readprefixlen]=='.', "Inproperate readprefixlen. Try to set readidfun and readprefixlen for the specific SAM file.")
                    if tc[1][readprefixlen]==readiddm
                        readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],':'))
                    else
                        # readidfun=X->parse(Int32,X[readprefixlen+1:end])
                        error("readprefixlen($readprefixlen) is not point to a readiddm($readiddm)")
                    end
                end
            end
            READID=readidfun(tc[1])
        catch err
            println("Error in parsing read-ID with prefix length $readprefixlen: "*tc[1])
            rethrow(err)
        end
        
        RECORD_NO+=1
        # FLAG=parse(Int16,tc[2])
        if pairend
            PAIRNO=Int8((FLAG&128>0)+1)
        else
            PAIRNO=Int8(0)
            # READID=[PAIRNO;READID]
        end
        READLEN=Int32(length(tc[10])) #Add in 23 Sep 2014

        if isaddfields
            ADDFDS=map(x->ascii(tc[x]),addfields[2])
        else
            ADDFDS=ASCIIString[]
        end
        
        #Parse Appendences
        AS=Int32(-1)
        NH=Int16(-1);HI=Int16(-1);nM=Int16(-1)
        for ctc=tc[12:end]
            if ctc[1:5]=="AS:i:"
                AS=parse(Int32,ctc[6:end])
            elseif ctc[1:5]=="NH:i:"
                NH=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="HI:i:" #Tophat have no HI.
                HI=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="nM:i:" || ctc[1:5]=="NM:i:" #STAR using nM for mismatch number, but Tophat use NM.
                nM=parse(Int16,ctc[6:end])
            end
        end
        
        #Parse CIGAR
        STR=int(tc[4])
        POS=zeros(Int,1,2)
        LOC=zeros(Int,1,2)
        # t=zeros(1,MAXMIS)
        # INS=copy(t);DEL=copy(t);INSLEN=copy(t);DELLEN=copy(t)
        # Ip=1;Dp=1;
        cigar=tc[6];p=1;rdp=1;rfp=STR;
        # aL=find(['A'.<=x.<='Z' for x in cigar])
        aL=findall(['A'<=x<='Z' for x in cigar])
        LOC_read_start=[];
        cjuncno=Int16(0);
        for cL=aL
            NO=int(cigar[p:cL-1]);p=cL+1;
            ccL=cigar[cL]
            if ccL=='S' || ccL=='H' #Added H in Dec 14, 2015. H used in BWA.
                rdp+=NO
                #S is not be calculated in the mapping start position.
                ## Remove this: ## rfp+=NO
                #Changed in Dec 14, 2015
            elseif ccL=='M'
                if !(POS[1]>0)
                    LOC[1]=rdp
		    if isempty(LOC_read_start)
			LOC_read_start=rdp
		    end
                    POS[1]=rfp
                end
                rfp+=NO
                rdp+=NO
                LOC[2]=rdp-1
                POS[2]=rfp-1
            elseif ccL=='D'
                # if Dp<=MAXMIS
                #     DEL[Dp]=rdp
                #     DELLEN[Dp]=NO
                # end
                # Dp+=1
                rfp+=NO;
            elseif ccL=='I'
                # if Ip<=MAXMIS
                #     INS[Ip]=rdp
                #     INSLEN[Ip]=NO
                # end
                # Ip+=1
                rdp+=NO
            elseif ccL=='N'
		cjuncno+=1
                addrow!(rlt,RECORD_NO,READID,PAIRNO,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,ADDFDS...)
                rfp+=NO
                POS=zeros(1,2)
            else
                error("Unknown cigar char in read: $ccL")
            end
        end
        if cjuncno>0
            cjuncno+=1
        end
        if (!onlyjunc) || cjuncno>0
            addrow!(rlt,RECORD_NO,READID,PAIRNO,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,ADDFDS...)
        end
        if eof(fp)
            break
        else
            ctxt=chomp(readline(fp))
        end
        
        if mod(RECORD_NO,1000000)==0
            println("Reading record: ", RECORD_NO);
        end
    end

    rltvalue=value(rlt)
    @assert(rltvalue!=nothing,"No any valid records. Please check, e.g., if the chromosome fields (#3) is valid.")
    read=tb([c"record_no,readid,pairno,flag,chrno,pos,juncno,align_score,mapped_loci_num,alignment_index,mismatch_num,read_len";addfields[1]],rltvalue)

    # Tophat and Bowtie are not support HI field.
    if !keepHI && all(read["alignment_index"]==-1)
        delete!(read,"alignment_index")
    end
    
    if onlyuniq
        println("Finding unique-mapped reads....")
        read=rec(read,sortri((read["readid"],read["pairno"],read["record_no"],-read["align_score"])))
        read["is_uniq_map"]=fill(true,rownum(read))
        lkeep=fill(true,rownum(read))

        preadid=pairend ? d"readid,pairno"read : read["readid"]
        cscore=read["align_score"][1]
        flag=true
        recordstrp=1
        for i=2:rownum(read)
            # if read["readid"][i,:]==read["readid"][i-1,:]
            if getrow(preadid,i)==getrow(preadid,i-1)
                if read["record_no"][i,:]==read["record_no"][i-1,:]
                    read["is_uniq_map"][i]=read["is_uniq_map"][i-1]
                    lkeep[i]=lkeep[i-1]
                else
                    if flag && read["align_score"][i]==cscore
                        read["is_uniq_map"][recordstrp:i].=false
                    else
                        flag=false
                        lkeep[i]=false
                    end
                    recordstrp=i
                end
            else
                flag=true
                recordstrp=i
                cscore=read["align_score"][i]
            end
        end
        read=rec(read,lkeep)
        if onlyjunc
            error("WTF")
            read["is_uniq_map"]=read["mapped_loci_num"].==1
        end
    end


    if !pairend
        delete!(read,"pairno")
    end

    if strand
        read["strand"]=map(x->x&16==0 ? '+' : '-',read["flag"])
        if pairend
            #Reverse second half-read in pair-end mode
            read["strand"][read["pairno"].==2]=map(x->x=='-' ? '+' : '-',read["strand"][read["pairno"].==2])
        end
        read=rec(read,sortri(d"chrno,strand,pos"read))
    else
        read=rec(read,sortri(d"chrno,pos"read))        
    end
    if !isempty(savefile)
        println("Saving file to $savefile ....")
        jldsave(savefile,read)
    end
    return read
end
function readStarSamFile(filename::Union{AbstractString,Cmd};region::Tuple=(), onlyuniq::Bool=false, param...)
    if isa(filename,AbstractString) && ismatch(r"\.bam$",filename)
        if isempty(region)
            file=`samtools view $filename`
            println("Reading BAM file ....")
        else
            if onlyuniq
                error("In region mode, onlyuniq=true is not supported. Considering set rd[\"is_uniq_map\"]=rd[\"mapped_loci_num\"].==1")
            end
            #Check if bam file is indexed
            if !isfile(filename[1:end-1]*"i") && !isfile(filename*".bai")
                println("Indexing bam file $filename")
                run(`samtools index $filename`)
                println("Indexing done.")
            end
            if length(region)==1
                #Read all reads by chromosomes
                chrs=isa(region[1],Group) ? unival(region[1]) : [region[1]]
                if eltype(chrs)<:Real
                    #Check if chromosome in bam file is started with chr
                    S=readchomp(pipeline(`samtools view -H $filename`,`grep '^@SQ'`))
                    M=match(r"\bSN\:(chr|CHR|Chr|chromosome|CHROMOSOME|Chromosome)",S)
                    if M!=nothing
                        chrs=no2chr(chrs,prefix=M.captures[1])
                    else
                        chrs=no2chr(chrs,prefix="")
                    end
                end
                T,=rowfunwith(chrs) do chr
                    file=`samtools view $filename $chr`
                    println("Reading BAM file chromosome $chr ....")
                    open((x)->readStarSamFile(x;keepHI=true,param...),file)
                end
                return T
            else
                #Read reads in some regions
                if isa(region[1],Group)
                    region,=segunion(region)
                else
                    region=([region[1]],region[2])
                end
                if eltype(region[1])<:Real
                    #Check if chromosome in bam file is started with chr
                    # fstbamli=open(readline,pipeline(`samtools view $filename`,`head -n1`))
                    # S=split(fstbamli,'\t')[3]|>x->x[1:min(3,length(x))]
                    # if lowercase(S)=="chr"
                    #     region=(no2chr(region[1],prefix=S), Array{Int}(region[2]))
                    # else
                    #     region=(no2chr(region[1],prefix=""), Array{Int}(region[2]))
                    # end
                    #Check if chromosome in bam file is started with chr
                    S=readchomp(pipeline(`samtools view -H $filename`,`grep '^@SQ'`))
                    M=match(r"\bSN\:(chr|CHR|Chr|chromosome|CHROMOSOME|Chromosome)",S)
                    if M!=nothing
                        region=(no2chr(region[1],prefix=M.captures[1]), Array{Int}(region[2]))
                    else
                        region=(no2chr(region[1],prefix=""), Array{Int}(region[2]))
                    end
                end
                T,=rowfunwith(region) do crg
                    rgstr="$(crg[1]):$(crg[2][1])-$(crg[2][2])"
                    file=`samtools view $filename $rgstr`
                    println("Reading BAM file region $rgstr ....")
                    open((x)->readStarSamFile(x; keepHI=true, onlyuniq=false, param...), file)
                end
                return T
            end
        end
    else
        @assert(isempty(region) || isempty(region[1]),"region parameter only support bam file so far.")
        file=filename
        println("Reading SAM file ....")
    end
    open((x)->readStarSamFile(x;param...),file)
end
export readStarSamFile

#readStarSamLine
export readStarSamLine
readStarSamLine(samline::AbstractString,args...;wargs...)=readStarSamLine([samline],args...;wargs...)
function readStarSamLine(samlines::Vector{T};
                         pairend::Bool=false,
                         strand::Bool=true,
                         readprefixlen::Int=-2,
                         readidfun::Union{Function,Void}=nothing,
                         chrfun::Function=x->Int16(chr2no(x)),
                         recordno_ref::Vector{Int}=[0],
                         addfields::Tuple{Vector,Vector}=([],[])) where {T<:AbstractString}
    RECORD_NO=recordno_ref[1]
    isaddfields=!isempty(addfields[1])
    rlt=rowpile()
    # SEQ=""
    # CIGAR4save=""
    for ctxt in samlines
        tc=split(chomp(ctxt),'\t')
        if length(tc)<3
            @warn("Ignore invalid line (record NO: $RECORD_NO): \"$ctxt\"")
            continue
        end
        FLAG=parse(Int16,tc[2])
        if FLAG & 0x4 != 0 || (CHRNO=chrfun(tc[3]))==0
            continue
        end
        
        #Parse READID
        if isa(readidfun,Void)
            if readprefixlen<0
                t=match(r"\:[\:\d]+$",tc[1])
                if t==nothing
                    readprefixlen=findonly(collect(tc[1]).=='.')
                    readidfun=X->parse(Int32,X[readprefixlen+1:end])
                else
                    readprefixlen=length(tc[1])-length(t.match)+1
                    readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],':'))
                end
            else
                @assert(tc[1]==0 || tc[1][readprefixlen]==':' || tc[1][readprefixlen]=='.', "Inproperate readprefixlen. Try to set readidfun and readprefixlen for the specific SAM file.")
                if tc[1][readprefixlen]==':'
                    readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],':'))
                else
                    readidfun=X->parse(Int32,X[readprefixlen+1:end])
                end
            end
        end
        READID=readidfun(tc[1])
        
        # if isa(readidfun,Void)
        #     if readprefixlen<0
        #         if readprefixlen<-1
        #             error("readprefixlen must be set in readStarSamLine. If you want it be signed automatically, set it as -1.")
        #         end
        #         readprefixlen=length(tc[1])-length(match(r"\:[\:\d]+$",tc[1]).match)+1
        #     end
        #     # READID=readidfun(tc[1],readprefixlen,readsplitchar)
        #     @assert(tc[1]==0 || tc[1][readprefixlen]==':', "Inproperate readprefixlen. Try to set readidfun and readprefixlen for the specific SAM file.")
        #     READID=map(x->parse(Int32,x),split(tc[1][readprefixlen+1:end],':'))
        # else
        #     if readprefixlen>0
        #         tc1=tc[1][readprefixlen+1:end]
        #     else
        #         tc1=tc[1]
        #     end
        #     READID=readidfun(tc1)
        # end
        
        RECORD_NO+=1
        # if getseq
        #     SEQ=ascii(tc[10])
        # end
        # if getCIGAR
        #     CIGAR4save=ascii(tc[6])
        # end
        if pairend
            PAIRNO=(FLAG&128>0)+1
            READID=[PAIRNO;READID]
        end
        READLEN=Int32(length(tc[10])) #Add in 23 Sep 2014

        if isaddfields
            ADDFDS=map(x->ascii(tc[x]),addfields[2])
        else
            ADDFDS=ASCIIString[]
        end
        
        #Parse Appendences
        AS=Int32(-1)
        NH=Int16(-1);HI=Int16(-1);nM=Int16(-1)
        for ctc=tc[12:end]
            if ctc[1:5]=="AS:i:"
                AS=parse(Int32,ctc[6:end])
            elseif ctc[1:5]=="NH:i:"
                NH=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="HI:i:" #Tophat have no HI.
                HI=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="nM:i:" || ctc[1:5]=="NM:i:" #STAR using nM for mismatch number, but Tophat use NM.
                nM=parse(Int16,ctc[6:end])
            end
        end
        
        #Parse CIGAR
        STR=int(tc[4])
        POS=zeros(Int,1,2)
        LOC=zeros(Int,1,2)
        cigar=tc[6];p=1;rdp=1;rfp=STR;
        aL=findall(['A'<=x<='Z' for x in cigar])
        LOC_read_start=[];
        cjuncno=Int16(0);
        for cL=aL
            NO=int(cigar[p:cL-1]);p=cL+1;
            ccL=cigar[cL]
            if ccL=='S' || ccL=='H'
                rdp+=NO
                #S is not be calculated in the mapping start position.
                ## Remove this: ## rfp+=NO
                #Changed in Dec 14, 2015
            elseif ccL=='M'
                if !(POS[1]>0)
                    LOC[1]=rdp
		    if isempty(LOC_read_start)
		        LOC_read_start=rdp
		    end
                    POS[1]=rfp
                end
                rfp+=NO
                rdp+=NO
                LOC[2]=rdp-1
                POS[2]=rfp-1
            elseif ccL=='D'
                rfp+=NO;
            elseif ccL=='I'
                rdp+=NO
            elseif ccL=='N'
	        cjuncno+=1
                addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,SEQ,CIGAR4save)
                # SEQ="" #For junction reads, only the first one assigned with sequence.
                # CIGAR4save=""
                rfp+=NO
                POS=zeros(1,2)
            else
                error("Unknown cigar char in read: $ccL")
            end
        end
        if cjuncno>0
            cjuncno+=1
        end
        addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,ADDFDS...)
        # addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,SEQ,CIGAR4save)
        # SEQ=""
        # CIGAR4save=""
    end
    rltvalue=value(rlt)
    if rltvalue==nothing
        return nothing
    end
    read=tb([c"record_no,readid,flag,chrno,pos,juncno,align_score,mapped_loci_num,alignment_index,mismatch_num,read_len";addfields[1]],rltvalue)
    # read=tb(c"record_no,readid,flag,chrno,pos,juncno,align_score,mapped_loci_num,alignment_index,mismatch_num,read_len,sequence,CIGAR",rltvalue)
    # if !getseq
    #     delete!(read,"sequence")
    # end
    # if !getCIGAR
    #     delete!(read,"CIGAR")
    # end
    if pairend
        read["pairno"]=Array{Int8}(read["readid"][:,1])
        read["readid"]=read["readid"][:,2:end]
    end
    if strand
        read["strand"]=map(x->x&16==0 ? '+' : '-',read["flag"])
        if pairend
            #Reverse second half-read in pair-end mode
            read["strand"][read["pairno"].==2]=map(x->x=='-' ? '+' : '-',read["strand"][read["pairno"].==2])
        end
    end
    recordno_ref[1]=RECORD_NO
    return read
end

#}}

#{{ parseBamMutation
#[[ parseBamMutation ]]
# T = parseBamMutation(Cigar, MD, read_seq; vcf_indel=false)
# Parse the mutation information according to the Cigar and MD fields.
# Output a table with fields: mut_typ, ref_loc, read_loc, ref_allele, alt_allele
# Where mut_typ is a charactor vector indicating the SNP ('S'), Insert ('I') or Delete ('D') mutations. ref_loc and read_loc are N x 2 matrix showing the mutation happened at the position on reference sequences (based on the left alignment point N0) and on read respectively. The absolute position of mutation can be calculated by N0+ref_loc-1.
# If vcf_indel is true, the left one nt adjacent to indel region will be included in the record.
#See also: readStarSamFile, readStarSamLine
#Xiong Jieyi, 7 Aug 2017

export parseBamMutation
function parseBamMutation(cigar::String, MD0::String, read::String; vcf_indel::Bool=false)
    rpl=rowpile()

    cgN,cgC=mapr(eachmatch(r"(\d+)([A-Z])",cigar)) do M
        (int(M.captures[1]), M.captures[2])
    end
    cg_readloc=zeros(Int,length(cgN),2)
    cg_refloc=zeros(Int,length(cgN),2)
    readIdxNoIns=collect(1:length(read)) #The index in MD is not include inserts and the "S" part of CIGAR
    rdp=1
    rfp=1
    for i=1:length(cgN)
        if cgC[i]=="S"
            cg_readloc[i,:]=[rdp,rdp+cgN[i]-1]
            cg_refloc[i,:]=[rfp,rfp-1]
            rdp+=cgN[i]
            deleteat!(readIdxNoIns, cg_readloc[i,1].<=readIdxNoIns.<=cg_readloc[i,2])
        elseif cgC[i]=="M"
            cg_readloc[i,:]=[rdp,rdp+cgN[i]-1]
            cg_refloc[i,:]=[rfp,rfp+cgN[i]-1]
            rdp+=cgN[i]
            rfp+=cgN[i]
        elseif cgC[i]=="I"
            cg_readloc[i,:]=[rdp,rdp+cgN[i]-1]
            cg_refloc[i,:]=[rfp,rfp-1]
            rdp+=cgN[i]
            addrow!(rpl, 'I', cg_refloc[i,:], cg_readloc[i,:], "", read[cg_readloc[i,:]|>x->x[1]:x[2]])
            deleteat!(readIdxNoIns, cg_readloc[i,1].<=readIdxNoIns.<=cg_readloc[i,2])
        elseif cgC[i]=="D"
            cg_refloc[i,:]=[rfp,rfp+cgN[i]-1]
            cg_readloc[i,:]=[rdp,rdp-1]
            rfp+=cgN[i]
        else
            error("Unknown GIGAR $cgC")
        end
    end
    @assert sum(seglen(cg_readloc))==length(read)
    
    MD=replace(MD0, r"(\d+[A-Z])0([A-Z])" => s"\1\2")
    mts=eachmatch(r"(\d+)(\^?[A-Z]+)",MD)
    if !isempty(mts)
        mdN,mdC=mapr(mts) do M
            (int(M.captures[1]), ascii(M.captures[2]))
        end
        cgi=1
        rdp_noins=1
        for i=1:length(mdN)
            rdp_noins+=mdN[i]
            rdp=readIdxNoIns[rdp_noins]
            
            if mdC[i][1]=='^'
                while cg_readloc[cgi,1]<rdp
                    cgi+=1
                end
                @assert cg_readloc[cgi,2]==rdp-1
                addrow!(rpl, 'D', cg_refloc[cgi,:], cg_readloc[cgi,:], mdC[i][2:end], "")
            else
                while cg_readloc[cgi,2]<rdp
                    cgi+=1
                end
                @assert cgC[cgi]=="M" && cg_readloc[cgi,1]<=rdp
                len=length(mdC[i])
                offset=cg_refloc[cgi,1]-cg_readloc[cgi,1]
                addrow!(rpl, 'S', [rdp rdp+len-1]+offset, [rdp rdp+len-1], mdC[i], read[rdp:rdp+len-1])
                rdp_noins+=len
            end
        end
        @assert int(match(r"\d+$",MD).match)+rdp_noins-1==length(readIdxNoIns)
    end
    
    T=tb(c"mut_typ, ref_loc, read_loc, ref_allele, alt_allele", value(rpl))

    if vcf_indel
        #Add the very left NT for INDEL
        for i=1:rownum(T)
            if T["mut_typ"][i]=='I' || T["mut_typ"][i]=='D'
                C=read[[T["read_loc"][i,1]-1]]
                T["ref_allele"][i]=C*T["ref_allele"][i]
                T["alt_allele"][i]=C*T["alt_allele"][i]
                T["ref_loc"][i,1]-=1
                T["read_loc"][i,1]-=1
            end
        end
    end
    
    T
end

#}}

#{{ parseCIGARpos
#[[ parseCIGARpos ]]
# pos_mx = parseCIGARpos(cigar::String, start_position::Int)
# Parse CIGAR code in bam/sam file, and output the position segments. If N exists in the CIGAR, the output will have multi-rows for each segment of junction read. This funciton is useful along with the XAM.jl.
# If cigar don't contain any mapping information, an error 'Ummapped CIGAR' will reported.
# See also:
# Xiong Jieyi, 17 Apr 2019 > 10 Jun 2021

export parseCIGARpos
function parseCIGARpos(cigar::AbstractString, STR::Int)
    POS=[0 0]
    LOC=[0 0]
    p=1;rdp=1;rfp=STR;
    aL=findall(['A'<=x<='Z' for x in cigar])
    LOC_read_start=[];
    cjuncno=Int(0);
    outpos=rowmx(Int, 2)
    for cL=aL
        NO=int(cigar[p:cL-1]);p=cL+1;
        ccL=cigar[cL]
        if ccL=='S' || ccL=='H'
            rdp+=NO
            #S is not be calculated in the mapping start position.
            ## Remove this: ## rfp+=NO
            #Changed in Dec 14, 2015
        elseif ccL=='M'
            if !(POS[1]>0)
                LOC[1]=rdp
	        if isempty(LOC_read_start)
		    LOC_read_start=rdp
	        end
                POS[1]=rfp
            end
            rfp+=NO
            rdp+=NO
            LOC[2]=rdp-1
            POS[2]=rfp-1
        elseif ccL=='D'
            rfp+=NO;
        elseif ccL=='I'
            rdp+=NO
        elseif ccL=='N'
	    cjuncno+=1
            if POS==[0 0]
                error("Invalid CIGAR \"$cigar\".")
            end
            addrow!(outpos, POS)
            rfp+=NO
            POS=zeros(1,2)
        else
            error("Unknown cigar char in read: $ccL")
        end
    end
    if POS==[0 0]
        error("Unmapped CIGAR \"$cigar\".")
    end
    addrow!(outpos, POS)
    value(outpos)
end
#}}

#{{ parseGTF
#[[ parseGTF ]]
# exon_table = parseGTF( GTF_filename[, additional_field1, additional_field2, ... ])
# Parse GTF file. Output will contain below fields:
#       chrno, exn_pos, cds_pos, strand,
#       geneid="",
#       transid="",
#       exn_rank=0,
#       gene_name="",
#       gene_typ="",
#       trans_name="",
#       protid="",
#       exnid=""
# Only "exon" rows in GTF files are parsed. Any other rows are ignored.
# `additional_field` is either a string for the name of field, or a NamedTuple in the format as (<field="...">, rename="...", fun=ascii, default=""). All keywords except `field=` are optional. For the rows missed this field, the default value will be filled.
#See also: readStarSamFile, genome_fasta2jld, parseDemulReport
#Xiong Jieyi, May 21, 2015 >Jul 30, 2015 >Aug 8, 2015 >31 Jan 2024

export parseGTF
function parseGTF(fn, field1::Union{AbstractString, NamedTuple}, fields::Union{AbstractString, NamedTuple}...)
    function fielddata0(;field::AbstractString="", rename::AbstractString=field, fun=ascii, default="")
        isempty(field) && error("If 2-N arguments is given as NamedTuple, key `field=' is mandatory.")
        (rename, field, fun, default)
    end
    fielddata(field::AbstractString)=(field, field, ascii, "")
    fielddata(field::NamedTuple)=fielddata0(; field...)
    parseGTF(fn, fielddata(field1), fielddata.(fields)...)
end
function parseGTF(fn, fields::Tuple{AbstractString,AbstractString,Function,Any}...)
    function phasefd(X)
        m=match(r"^\s*(\S+)\s+\"([^\"]*)\"\s*$",X)
        if m==nothing
            m=match(r"^\s*(\S+)\s+(\S*)\s*$",X)
        end
        if m==nothing
            @warn("Invalid item: $X")
            ("","")
        else
            (m.captures[1],m.captures[2])
        end
    end
    fielddict=Dict{AbstractString, Tuple{AbstractString,Function}}("gene_id"=>("geneid",ascii),
               "transcript_id"=>("transid",ascii),
               "exon_number"=>("exn_rank",int),
               "gene_name"=>("gene_name",ascii),
               "gene_biotype"=>("gene_typ",ascii),
               "gene_type"=>("gene_typ",ascii),
               "transcript_name"=>("trans_name",ascii),
               "protein_id"=>("protid",ascii),
               "exon_id"=>("exnid",ascii))
    empexn=tb(geneid="",
              transid="",
              exn_rank=0,
              gene_name="",
              gene_typ="",
              trans_name="",
              protid="",
              exnid="")
    for c in fields
        # @assert(isa(c,(AbstractString,AbstractString,Function,Any)),"Invalid additional fields. The input should be (OutFieldName::AbstractString, GTFFieldName::AbstractString, Function, DefaultValue::Any).")
        fielddict[c[2]]=(c[1],c[3])
        empexn[c[1]]=c[4]
    end
    exn=rowpile()
    cds=rowpile()
    fileloop(fn, ignore='#')do li
        C=split(rstrip(li),'\t')
        if C[3]!="exon" && C[3]!="CDS"
            return nothing
        end
        chrno=chr2no(C[1])
        if chrno==0
            return nothing
        end
        if length(C)>9
            Ced=join(C[9:end]," ")
        else
            Ced=C[9]
        end
        D=split(rstrip(Ced,';'),';')
        cexn=copy(empexn)
        d"chrno,exn_pos,strand"cexn=(chrno,int(C[4:5])',C[7][1])
        for cfd in D
            k,v=phasefd(cfd)
            if isempty(k)
                return nothing #When parsed an invalid item, give up this record.
            end
            if haskey(fielddict,k)
                kk,fun=fielddict[k]
                cexn[kk]=fun(v)
            end
        end
        if C[3]=="exon"
            addrow!(exn,cexn)
        else
            # C[3]=="CDS"
            addrow!(cds,cexn)
        end
     end
    exn=value(exn)
    cds=value(cds)
    if cds!=nothing
        (ie,ic)=genomemap(d"transid,exn_pos"exn,d"transid,exn_pos"cds)
        exn["cds_pos"]=zeros(Int,rownum(exn),2)
        exn["cds_pos"][ie,:]=cds["exn_pos"][ic,:]
        exn["protid"][ie,:]=cds["protid"][ic,:]
    end
    return exn
end
#}}

#{{ writeGTF
#[[ writeGTF ]]
#writeGTF(filename, chr|chrno, pos, geneid; strand='.', prefix="chr", trans="transcript_id")
#Write a simplest GTF file. This GTF file can be used by HTSeq-count.
#prefix is used in no2chr.
#When trans="..." is given, a tag like ' transcript_id "...";' will be added to each line. It is not necessary for HTSeq-count.
#See also:no2chr
#Xiong Jieyi, May 15, 2016 > 28 Mar 2024

export writeGTF
function writeGTF(filename::AbstractString, chr::AbstractVector{<:Union{Integer, AbstractString}}, pos::AbstractMatrix{<:Integer}, geneid::AbstractVector{<:AbstractString}; strand::Union{Nothing, AbstractVector{Char}}=nothing, prefix="chr", trans::Union{AbstractVector{<:AbstractString}, Nothing}=nothing)
    if eltype(chr)<:Integer
        chr=no2chr(chr, prefix=prefix)
    end
    if isnothing(strand)
        strand=fill('.', length(chr))
    end
    addtxts=if isnothing(trans)
        fill("", length(chr))
    else
        f" transcript_id \"$1\";".(trans)
    end
    open(filename, "w") do fid
        rowloop(chr,pos,geneid,strand,addtxts) do cchr,cpos,cgeid,cstrand,addtxt
            println(fid,"$cchr\t.\texon\t$(cpos[1])\t$(cpos[2])\t.\t$cstrand\t.\tgene_id \"$cgeid\";$addtxt")
        end
    end
    nothing
end
#}}

#{{ segunion

#[[ segunion ]]
# (seg_pos, score)=segunion(pos; score=1s, weld=true, mergegap=0)
# ((seg_chr, seg_pos), score)=segunion(chr, pos;...)
# ((seg_chr, seg_pos, seg_strand), score)=segunion(chr, pos, strand;...)|segunion((chr,pos,strand);...)
#Union the segments. The output score is the sum of inputed scores for each union segment.
#If weld is true, [1 3] will be merged with [4 6].
#If mergegap>0, the map no longer than this cutoff will be merged.
# e.g. segunion([1 3; 6 8], mergegap=2) => ([1 8], [2])
# but  segunion([1 3; 6 8], mergegap=2, weld=false) => ([1 3; 6 8], [1, 1])
#See also: segconc, seg2cvg
#Xiong Jieyi, 1 Sep 2014 >October 21, 2014 >2 Mar 2017>14 May 2018>22 Jan 2021

function segunion(pos::Matrix{T};score::AbstractVector{T2}=ones(Int,size(pos,1)), weld::Bool=true, mergegap::Int=0) where {T<:Real,T2<:Real}
    mergegap>=0 || error("mergegap should >=0.")
    if isempty(pos)
        return zeros(T,0,2)
    end
    len=size(pos,1)
    W=mergegap + (weld ? 1 : 0)
    t,idx=sortr([pos[:,1]-W ones(T,len); pos[:,2] 2*ones(T,len)])
    st=t[:,1]
    yy=t[:,2]
    yy[yy.==2].=-1
    dscore=[score;zeros(typein(score),length(score))][idx]
    
    cur=0
    count=0
    bd=st[1]+W
    rlt=rowpile()
    for i=1:length(st)
        cur+=yy[i]
        if cur==0
            addrow!(rlt,[bd, st[i]],count)
            count=0
            maxcover=0
            if i<length(st)
                bd=st[i+1]+W
            end
        elseif yy[i]==1
            count+=dscore[i]                
        end
    end
    return value(rlt)
end
#<v0.6# function segunion{T<:Real}(chr::Group,pos::Array{T};args...)
function segunion(chr::Group,pos::Array{T};args...) where {T<:Real}
    (seg_pos,seg_score),seg_chr=grpfunwith(x->segunion(x;args...),chr,pos)
    return ((seg_chr,seg_pos),seg_score)
    # O,t=grpfun(x->segunion(x;args...),chr,pos)
    # seg_chr,t=vcatr_with(t,O...)
    # return ((seg_chr,t[1]),t[2])
end
#<v0.6# function segunion{T<:Real}(chr::Array,pos::Array{T},strand::Array;args...)
function segunion(chr::Array,pos::Array{T},strand::Array;args...) where {T<:Real}
    (seg_pos,seg_score),(seg_chr,seg_strand)=grpfunwith(x->segunion(x;args...),(chr,strand),pos)
    return ((seg_chr,seg_pos,seg_strand),seg_score)
    # O,t=grpfun(x->segunion(x;args...),(chr,strand),pos)
    # tt,t=vcatr_with(t,O...)
    # return ((tt[1],t[1],tt[2]),t[2])
end
segunion(X::Tuple;args...)=segunion(X...;args...)
export segunion

#}}

#{{ findintron
#[[ findintron ]]
#intron_pos = findintron(exon_pos; out_jc_pos=false)
#(transid, intron_pos)=findintron(transid, exon_pos; out_jc_pos=false)
#Looking for intron using given exon position.
#When out_jc_pos=true, the function output junction position, i.e., the boundaries of exons. Which equals to [intron_pos1.-1 intron_pos2.+1].
#See also: segunion, segconc, findorf
#Xiong Jieyi, 16 Sep 2014 > 6 Apr 2020

export findintron
function findintron(pos::Matrix{T}; out_jc_pos::Bool=false) where {T<:Real}
    allpos=[minimum(pos[:,1]) maximum(pos[:,2])]
    outpos=segconc(allpos,1,pos,-1)[1]
    if out_jc_pos
        outpos=[outpos[:, 1].-1 outpos[:, 2].+1]
    end
    outpos
end
function findintron(transid::Group, pos::Matrix{T}; out_jc_pos::Bool=false) where {T<:Real}
    grpfunwith(x->findintron(x; out_jc_pos=out_jc_pos), transid, pos)[[2, 1]]
end
#}}

#{{ getseq
#[[ getseq ]]
#Seqs = getseq((chr, pos[, strand]); path="path_of_sequence_lib", fun=x->typeof(x)[x], funarg2=...)
#Or Seqs = getseq( fun(sequence[, getrow(funarg2, i)]), ...)
#Seqs = getseq(pos[, strand];path="path_of_sequence_lib",fun=x->x)
#Fetch the sequences from reference genomes. If the strand is given, the '-' sequence will be reverse-complement transferred. You can build the path_of_sequence_lib by genome_fasta2jld.
#pos is a two-column integer matrix, but can also be a vector of coordinate vector.
#When funarg2 is assigned and not nothing, the given function will receive two inputs: the sequence, and cooresponding row slice of funarg2.
#See also: catseq, segunion, genome_fasta2jld
#Xiong Jieyi, 5 Sep 2014>February 26, 2015>Jan 1, 2016>29 Apr 2020>17 Oct 2020

function getseq(chrpos::Tuple{Group, AbstractArray};
                path::ASCIIString="", fun::Function=x->typeof(x)[x], funarg2::Union{Nothing, Group}=nothing)
    chr,pos=chrpos
    if isnothing(funarg2)
        grpfunexp(chr, pos; input_grpid=true) do cchr, cpos
            if isa(cchr,Real)
                cchr=no2chr(cchr)
            end
            file=joinpath(path,cchr)
            getseq(cpos, file; fun=fun)
        end
    else
        grpfunexp(chr, pos, funarg2; input_grpid=true) do cchr, cpos, cfunarg
            if isa(cchr,Real)
                cchr=no2chr(cchr)
            end
            file=joinpath(path,cchr)
            getseq(cpos, file; fun=fun, funarg2=cfunarg)
        end
    end
end
function getseq(chrposstrand::Tuple{Group, AbstractArray, AbstractVector{Char}};
                path::ASCIIString="", fun::Function=x->typeof(x)[x], funarg2::Union{Nothing, Group}=nothing)
    chr,pos,strand=chrposstrand
    if isnothing(funarg2)
        grpfunexp((chr, strand), pos; input_grpid=true) do cchr_strand, cpos
            cchr,cstrand=cchr_strand
            if isa(cchr,Real)
                cchr=no2chr(cchr)
            end
            file=joinpath(path,cchr)
            getseq(cpos, file; strand=cstrand, fun=fun)
        end
    else
        grpfunexp((chr, strand), pos, funarg2; input_grpid=true) do cchr_strand, cpos, cfunarg
            cchr,cstrand=cchr_strand
            if isa(cchr,Real)
                cchr=no2chr(cchr)
            end
            file=joinpath(path,cchr)
            getseq(cpos, file; strand=cstrand, fun=fun, funarg2=cfunarg)
        end
    end
end
function getseq(pos::AbstractMatrix{T},file::AbstractString; strand::Char='+',fun::Function=x->typeof(x)[x], funarg2::Union{Nothing, Group}=nothing) where {T<:Integer}
    BioSeq=importpy("Bio.Seq")
    seq=jldload(file, "seq") #Data still saved in JLD format, since JLD2 do not support long string. 15 Jun 2019
    rlt=rowpile()
    for i=1:size(pos,1)
        csq=seq[pos[i,1]:pos[i,2]]
        if strand=='-'
            csq=@il BioSeq.reverse_complement(csq)
        end
        funout=if isnothing(funarg2)
            fun(csq)
        else
            fun(csq, getrow(funarg2, i))
        end
        addrow!(rlt, funout)
    end
    return value(rlt)
end
function getseq(pos::AbstractVector{T}, file::AbstractString; strand::Char='+',fun::Function=x->typeof(x)[x], funarg2::Union{Nothing, Group}=nothing) where {T<:AbstractVector}
    BioSeq=importpy("Bio.Seq")
    seq=jldload(file, "seq")
    rlt=rowpile()
    for i=1:length(pos)
        csq=seq[pos[i]]
        if strand=='-'
            csq=@il BioSeq.reverse_complement(csq)
        end
        funout=if isnothing(funarg2)
            fun(csq)
        else
            fun(csq, getrow(funarg2, i))
        end
        addrow!(rlt, funout)
    end
    return value(rlt)
end
getseq(fun::Function,args...;wargs...)=getseq(args...;fun=fun,wargs...)
export getseq
#}}

#{{ genome_fasta2jld
#[[ genome_fasta2jld ]]
# genome_fasta2jld(fasta_filename, output_path)
#Convert genome fasta file to binrary sequence library. This path can be used by getseq.
#The reason this function still saves data using JLD.jl rather than JLD2.jl is JLD2.jl do not support storage long string so far (15 Jun 2019).
#See also: getseq, parseGTF, readStarSamFile
#Xiong Jieyi,February 26, 2015 > 14 May 2018

export genome_fasta2jld
function genome_fasta2jld(fn::AbstractString,pth::AbstractString)
    if !isdir(pth)
        println("Create path $pth")
        mkdir(pth)
    end
    open(fn)do fp
        chr=""
        iob=IOBuffer()
        while !eof(fp)
            ln=strip(chomp(readline(fp)))
            if ln[1]=='>'
                if !isempty(chr)
                    # seq=takebuf_string(iob)
                    seq=String(take!(iob))
                    jldsave(joinpath(pth,chr), "seq"=>seq, "chr"=>chr, "length"=>length(seq))
                end
                chr=ascii(split(strip(ln[2:end]))[1])
                iob=IOBuffer()
            else
                print(iob,ln)
            end
        end
        # seq=takebuf_string(iob)
        seq=String(take!(iob))
        jldsave(joinpath(pth,chr), "seq"=>seq, "chr"=>chr, "length"=>length(seq))
    end
end
#}}

#{{ catseq

#[[ catseq ]]
#(trans_seq, exn_loc_in_trans)=catseq(exn_seqs,exn_pos,strand='+')
#((trans_seq, transid), exn_loc_in_trans)=catseq(exn_seq,(transid,pos,strand))
#Caternate exon sequence together, and also return the location of each exon in the transcrition.
#exn_loc_in_trans is a N_exn x 2 matrix with the corresponding line as input.
#Supported sequence could be AbstractString or any vector.
#The coverlaps of exn_pos will be checked automatically, and a error will occur if happened.
#See also: getseq, segmask
#Xiong Jieyi, 4 Oct 2014>May 28, 2015>Aug 19, 2015

function catseq(seq::Vector{T},pos::Matrix,strand::Char='+') where {T}
    if size(pos,2)>1
        if strand=='+'
            idx=sortri(pos)
            spos=pos[idx,:]
            @assert(all(spos[1:end-1,2].<spos[2:end,1]),"Segments have overlap.")
        elseif strand=='-'
            idx=sortri(pos,rev=true)
            spos=pos[idx,:]
            @assert(all(spos[2:end,2].<spos[1:end-1,1]),"Segments have overlap.")
        else
            error("Strand could neither be + or -.")
        end
    else
        spos=pos
    end
    sseq=seq[idx]
    # spos=pos[idx,:]
    if T<:AbstractString
        trseq=*(sseq... )
    else
        trseq=[sseq...]
    end
    t=cumsum(spos[:,2]-spos[:,1]+1)
    trloc=[[1;t[1:end-1]+1] t][rvorder(idx),:]
    return (trseq,trloc)
end
function catseq(seq::Vector{T},trid_pos_strand::Tuple{Group,Matrix,Vector{Char}}) where {T}
    (transid,pos,strand)=trid_pos_strand
    trT=rowpile()
    trans_loc=grpfunexp((transid,strand),seq,pos;input_grpid=true) do trid_strand, cseq, cpos
        (trseq,trloc)=catseq(cseq,cpos,trid_strand[2])
        if T<:AbstractString
            addrow!(trT,trseq,trid_strand[1])
        else
            addrow!(trT,Any[trseq],trid_strand[1])
        end
        trloc
    end
    return (value(trT),trans_loc)
end
export catseq

#}}

#{{ findorf
#[[ findorf ]]
#pos_mx = findorf(Seq)
#Greedy search all open reading frame (ATG-TAA|TGA|TAG) in the sequence.
#See also: getseq, segmask
#Xiong Jieyi, 7 Sep 2014>7 Oct 2014

function findorf(Seq::ASCIIString)
    stri=1
    pos=rowmx()
    while (m=match(r"ATG(?:[ATGC]{3})*(?:TAA|TGA|TAG)"i,Seq,stri))!=nothing
        addrow!(pos,[m.offset,m.offset+length(m.match)-1])
        stri=m.offset+1
    end
    return value(pos)
end
export findorf
#}}

#{{ segmask
#[[ segmask ]]
#(mask_pos, mask_loc)=segmask(seg_pos, mask_start, mask_len; strand='+')
#(mask_pos, mask_loc)=segmask(seg_pos, [mask_start mask_end]; strand='+')
#Mask the given segments, and output the masked absolute-position and relative-location. 0 in output means this segment is not in the mask region. When strand='-', segment will be masked from right. The inputed seg_pos isn't needed to be sorted, and the output is always consistent with inputed seg_pos rowly-by-rowly.
#e.g., cds_pos=segmask(exn_pos,orf_loc_in_trans; strand=strand)
#See also: findorf, getseq, catseq
#Xiong Jieyi, 5 Oct 2014

#<v0.6# function segmask{T}(pos::Matrix{T},mask_str::Real,mask_len::Real;strand::Char='+')
function segmask(pos::Matrix{T},mask_str::Real,mask_len::Real;strand::Char='+') where {T}
    if strand=='-'
        (mask_pos, mask_loc)=segmask(-fliplr(pos),mask_str,mask_len;strand='+')
        return (-fliplr(mask_pos), mask_loc)
    end

    idx=sortri(pos)
    pos=pos[idx,:]
    
    mask_skp=mask_str-1
    mask_pos=zeros(T,size(pos))
    mask_loc=zeros(T,size(pos))
    seg_len=seglen(pos)
    p=1
    while p<=size(pos,1) && seg_len[p]<=mask_skp
        mask_skp-=seg_len[p]
        p+=1
    end
    next_mask_loc=1
    while p<=size(pos,1)
        mask_loc[p,1]=next_mask_loc
        mask_pos[p,1]=pos[p,1]+mask_skp
        if mask_len<=seg_len[p]-mask_skp
            mask_pos[p,2]=mask_pos[p,1]+mask_len-1
            mask_loc[p,2]=mask_loc[p,1]+mask_len-1
            break
        else
            mask_pos[p,2]=pos[p,2]
            mask_loc[p,2]=mask_loc[p,1]+seg_len[p]-mask_skp-1
            mask_len-=seg_len[p]-mask_skp
            next_mask_loc=mask_loc[p,2]+1
        end
        p+=1
        mask_skp=0
    end
    vidx=rvorder(idx)
    (mask_pos[vidx,:], mask_loc[vidx,:])
end
segmask(pos::Matrix,mask_pos::Array;wargs...)=segmask(pos,mask_pos[1],mask_pos[2]-mask_pos[1]+1;wargs...)
export segmask
#}}

#{{ segconc

#[[ segconc ]]
# (pos,score)=segconc(pos,score=1;cutoff=1)
# ((chr,pos),score)=segconc((chr,pos),score=1;cutoff=1)
# ((chr,pos,strand),score)=segconc((chr,pos,strand),score=1;cutoff=1)
# ...=segconc(A,score_A,B,score_B;cutoff)
#Sort, merge segments and tide scores. The result bed are sorted and unoverlapped.
#See also: segunion, seg2cvg, segintersect segmentdiff, segmask
#Xiong Jieyi, 8 Sep 2014>22 Sep 2014

function segconc(pos::Matrix{T1},
               score::Union{AbstractVector,Real}=1;cutoff::Real=1) where {T1<:Real}
    if isa(score,Real)
        score=fill(score,size(pos,1))
    end
    P0,t=sortr([pos[:,1];pos[:,2]+1])
    tt=[score;-score]
    S0=tt[t]
    ldel=fill(false,length(S0))
    if S0[1]==0 && cutoff>0
        ldel[1]=true
    end
    for i=2:length(S0)
        if P0[i]==P0[i-1]
            ldel[i-1]=true
            S0[i]+=S0[i-1]
        end
        if S0[i]==0 && (cutoff>0 || i<length(S0))
            ldel[i]=true
        end
    end
    P=P0[.!ldel]
    S=S0[.!ldel]

    rlt=rowpile()
    if isempty(S)
        return (zeros(typein(pos),0,2),ones(typein(score),0))
    else
        V=S[1]
        for i=2:length(S)
            if V>=cutoff
                addrow!(rlt,[P[i-1] P[i]-1],V)
            end
            V+=S[i]
        end
        ret=value(rlt)
        if ret==nothing
            return (zeros(typein(pos),0,2),ones(typein(score),0))
        else
            return ret
        end
    end
end
function segconc(chrpos::Tuple{Group,Matrix,AbstractVector},
                 score::Union{AbstractVector,Real}=1;cutoff::Real=1,args...)
    (chr,pos,strand)=chrpos
    if isa(score,Real)
        score=fill(score,size(pos,1))
    end
    # (pos_score,chr_strand)=grpfun((chr,strand),pos,score) do cpos,cscore
    #     segconc(cpos,cscore;args...)
    # end
    # ((ochr,ostrand),(opos,oscore))=vcatr_with(chr_strand,pos_score...)
    (opos,oscore),(ochr,ostrand)=grpfunwith((chr,strand),pos,score) do cpos,cscore
        segconc(cpos,cscore;args...)
    end
    return ((ochr,opos,ostrand),oscore)
end
function segconc(chrpos::Tuple{Group,Matrix},score::Union{AbstractVector,Real}=1;args...)
    (chr,pos)=chrpos
    if isa(score,Real)
        score=fill(score,size(pos,1))
    end
    # (pos_score,gchr)=grpfun(chr,pos,score) do cpos,cscore
    #     segconc(cpos,cscore;args...)
    # end
    # (gchr,(opos,oscore))=vcatr_with(gchr,pos_score...)
    (opos,oscore),ochr=grpfunwith(chr,pos,score) do cpos,cscore
        segconc(cpos,cscore;args...)
    end
    return ((ochr,opos),oscore)
end
function segconc(Ap::Union{Matrix,Tuple}, As::Union{AbstractVector,Real},
                 Bp::Union{Matrix,Tuple}, Bs::Union{AbstractVector,Real};args...)
    if !isa(As,AbstractVector)
        As=fill(As,rownum(Ap))
    end
    if !isa(Bs,AbstractVector)
        Bs=fill(Bs,rownum(Bp))
    end
    segconc(vcatr(Ap,Bp),vcatr(As,Bs);args...)
end
export segconc

#}}

#{{ segintersect segmentdiff
#[[ segintersect segmentdiff ]]
#(pos_mx|(chrno,pos[,strand]), Idx_segA)=segintersect|segmentdiff(segA,segB; BnoOverlap=false)
# pos_mx, Idx_segA = segintersect(A_pos_mx, B_pos_mx; BnoOverlap=false, Asorted=false, Bsorted=false, sortcheck=true)
#Get the seg parts of segA in/not in region of segB. Input could be pos_matrix or (chr,pos[, strand]) tuple.
#If BnoOverlap is ture, segintersect will not union B first. Note that when B have overlap, the results are wrong. At the moment, when BnoOverlap=false, the segunion() step will not be optimized even if Bsorted=true (5 Sep 2019).
#See also: segconc, segunion, seg2cvg, segmask
#Xiong Jieyi, February 8, 2015 > 5 Sep 2019

export segintersect, segmentdiff
function segintersect(segA::Tuple,segB::Tuple; BnoOverlap::Bool=false)
    if !BnoOverlap && rownum(segB)>1
        segB,=segunion(segB)
    end
    difseg,Ai=genomemap(segA,segB,touch=true,
              fun=(ai,bi,ap,bp)->([max(ap[1],bp[1]) min(ap[2],bp[2])], ai))
    if length(segA)==2
        difseg=tuple(getrow(segA[1],Ai),difseg)
    else
        @assert(length(segA)==3)
        difseg=tuple(getrow(segA[1],Ai),difseg,getrow(segA[3],Ai))
    end
    return (difseg, Ai)
end
function segintersect(segA::Matrix{T1},segB::Matrix{T2}; BnoOverlap::Bool=false, Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true) where {T1<:Real,T2<:Real}
    if !BnoOverlap && size(segB,1)>1
        segB,=segunion(segB)
    end
    coordinatemap(segA,segB,touch=true, no_match_output=(zeros(T1, 0, 2), Int[]),
                  fun=(ai,bi,ap,bp)->([max(ap[1],bp[1]) min(ap[2],bp[2])], ai),
                  sitesorted=Asorted, tagsorted=Bsorted, sortcheck=sortcheck)
end
function segmentdiff(segA,segB; BnoOverlap::Bool=false)
    if !BnoOverlap && size(segB,1)>1
        segB,=segunion(segB)
    end
    segD,=segconc(segunion(segA)[1],1,segB,-1)
    segintersect(segA,segD;BnoOverlap=true)
end
#}}

#{{ seg2cvg
#[[ seg2cvg ]]
#(coverage::Vector{Int}, X::Range) = seg2cvg(pos[, start, end]; max_pixel=0)
#coverage::Vector{Int} = seg2cvg(pos, x_range::Range)
#Calculate coverage using segments.
#If max_pixel >0 and end-start+1>max_pixel, coverage will be conpressed, and the length of output will be between 1/2 max_pixel and max_pixel.
#See also: segconc, segunion
#Xiong Jieyi, 10 Sep 2014>December 12, 2014>Mar 3, 2016

function seg2cvg(pos::Array{T},
                 limsta::Real=minimum(pos[:,1]),limend::Real=maximum(pos[:,2]); max_pixel::Integer=0) where {T<:Real}
    len=limend-limsta+1
    if 0<max_pixel<len
        N=ceil(Int,len/max_pixel)
    else
        N=1
    end
    X=limsta:N:limend
    cvg=seg2cvg(pos,X)
    return (cvg,X)
end
function seg2cvg(pos::Array{T}, xrang::Range) where {T<:Real}
    len=length(xrang)
    cvg=zeros(Int,len)
    rpos=round.(Int,(pos-first(xrang))/step(xrang)+1)
    for i=1:size(rpos,1)
        cvg[max(rpos[i,1],1):min(rpos[i,2],len)]+=1
    end
    return cvg
end
export seg2cvg
#}}

#{{ writefasta
#[[ writefasta ]]
#writefasta("filename"|FileIO, seqs, heads; append=false, gapline=false, width=0)
#Write sequence to fasta file.
#If gapline=true, a blank line will be added to the end of each sequence.
#If width>0, sequence longer than width will be broke into multiple lines.
#See also: getseq, readfasta
#Xiong Jieyi, 16 Sep 2014 > 28 Mar 2023

function writefasta(filename, seq::AbstractVector{<:AbstractString}, head::AbstractVector{<:AbstractString}; append::Bool=false, gapline::Bool=false, width::Integer=0)
    gap=gapline ? "\n\n" : "\n"
    open(filename, append ? "a" : "w") do fp
        for (cseq, chead) in zip(seq, head)
            write(fp, ">"*chead*"\n")
            if width>0
                p=1
                while p<=length(cseq)
                    q=p+width-1
                    if q<length(cseq)
                        write(fp, cseq[p:q]*"\n")
                    else
                        write(fp, cseq[p:end]*gap)
                    end
                    p+=width
                end
            else
                write(fp, cseq*gap)
            end
        end
    end
end
export writefasta
#}}

#{{ readfasta
#[[ readfasta ]]
#(seqeuence, head)=readfasta(filename|pipe_cmd; skipline::Int=0, ignore::Char='\0')
#Read fasta file. The line breaks and side blanks in the sequence will be cleaned.
#See also: writefasta, getseq
#Xiong Jieyi, 16 Sep 2014

# function readfasta(filename)
#     # needBioSeq()
#     importpkg(:PyCall, preloaded=true)
#     BioSeqIO=importpy("Bio.SeqIO")
#     oo=BioSeqIO.parse(filename,"fasta")
#     str(x)=Main.PyCall.pycall(Main.PyCall.pybuiltin("str"),Main.PyCall.PyAny,x)
#     t=map(oo) do x
#         (str(x.seq), x.id)
#     end
#     return vcatr(t...)
# end
export readfasta
function readfasta(fn; skipline::Int=0, ignore::Char='\0')
    head=String[]
    seq=String[]
    buf=IOBuffer()
    fileloop(fn, no=true, skipline=skipline, ignore=ignore) do no, li
        isempty(strip(li)) && return nothing
        if li[1]=='>'
            no>1 && push!(seq, String(take!(buf)))
            push!(head, strip(li[2:end]))
        else
            write(buf, strip(li))
        end
    end
    push!(seq, String(take!(buf)))
    (seq, head)
end
#}}

#{{ balancemask balancemask2
#[[ balancemask balancemask2 ]]
#lA, lB = balancemask(A, B; N=1~10)
#lB, times_B/A, B_weight = balancemask2(A_small, B_big; N=1~10, toofewbg=:error|:warn)
#Divide [A,B] into N group, and mask in order to keep the number of A and B balance in the two groups. A[lA] and B[lB] is balanced. (lA and lB are BitVector).
#balancemask2: function will try to keep as many as B while keep A and B in balance. When toofewbg=:error, #B could no more than #A in any group. When toofewbg=:warn or :ignore, lB in this bin is all trues, but the balance only by subsetting #B is not completed.+#B_weight is a probability(0~1) list with the same length as B, indicating the probability of each B element be choosen as background. In the enrichment analysis, instead filtering B by lB, you can use round(Int, sum(B_weitht[is_B_hit])) to estimate background hits, which has the advantge to eliminate the random effect of lB. In case toofewbg=:warn, if in any bins #A even more than #B, the output A/B will be forcely set to 1.
#see also: drawer, GOenrichFB, GOenrichGrp
#Xiong Jieyi, October 19, 2014>March 9, 2015>Nov 12, 2015>8 Sep 2018>14 Apr 2021>14 Mar 2022

function balancemask(A::Vector,B::Vector;N::Int=0)
    Idx=sortri([A;B])
    if N==0
        N=min(max(floor(Int,length(Idx)./10),1),10)
        println("Group by $N bins.")
    end
    sL=[trues(length(A)); falses(length(B))][Idx]
    oL=trues(length(sL))
    grpelnum=drawer(length(sL),N)
    p=1
    for gi=1:N
        l=p:p+grpelnum[gi]-1
        cL=sL[l]
        nA=sum(cL)
        nB=length(cL)-nA
        if nA>nB
            oL[l[find(cL)[randperm(nA)[1:(nA-nB)]]]].=false
        elseif nA<nB
            oL[l[find(.!cL)[randperm(nB)[1:(nB-nA)]]]].=false
        end
        p=p+grpelnum[gi]
    end
    ooL=oL[rvorder(Idx)]
    return ooL[1:length(A)],ooL[length(A)+1:end]
end
function balancemask2(A::Vector,B::Vector;N::Int=0, toofewbg::Symbol=:error)
    Idx=sortri([A;B])
    if N==0
        N=min(max(floor(Int,length(Idx)./10),1),10)
        println("Group by $N bins.")
    end
    sL=[trues(length(A)); falses(length(B))][Idx]
    oL=trues(length(sL))
    grpelnum=drawer(length(sL),N)
    p=1
    minratio=Inf
    for gi=1:N
        l=p:p+grpelnum[gi]-1
        cL=sL[l]
        nA=sum(cL)
        nB=length(cL)-nA
        cratio=nB/nA
        if cratio<minratio
            minratio=cratio
        end
        p=p+grpelnum[gi]
    end
    if minratio<1
        if toofewbg==:warn
            @warn("B in some group is even less than A (min B/A=$minratio). Under toofewbg=:warn, output B/A has been forcely set to 1.")
            minratio=1.0
        else
            @error("B in some group is even less than A (min B/A=$minratio). Try to decrease N, set toofewbg=:warn, or using balancemask instead.")
        end
    end
    sWeight=fill(NaN, length(Idx))
    p=1
    for gi=1:N
        l=p:p+grpelnum[gi]-1
        cL=sL[l]
        nA=sum(cL)
        nB=length(cL)-nA
        sWeight[l].=nA*minratio/nB
        oL[l[find(.!cL)[randperm(nB)[1:(nB-round(Int,nA*minratio))]]]].=false
        p=p+grpelnum[gi]
    end
    rvIdx=invperm(Idx)
    ooL=oL[rvIdx]
    ooWeight=sWeight[rvIdx]
    return ooL[length(A)+1:end], minratio, ooWeight[length(A)+1:end]
end
export balancemask, balancemask2
#}}

#{{ enrich_test
#[[ enrich_test ]]
# T = enrich_test(fg, bg, item, cata; fdr_cutoff=0.1, fgN_cutoff=2,
#                 explain=(cata, explain_text), nowarn=false,
#                 balance=(fg, bg) | bg_weight=weight_vector, 
#                 balance_method=:bg_weight|:bg_subsample|:fg_and_bg_subsample,
#                 list_fg_hit=true, list_fg_max=100,
#                 list_fg_idconv=(item2, new_id[, default_id]),
#                 list_fg_fun=x->join(sort(x), ',')::AbstractString,
#                 savefile=" tsv file to save ")
#Enrichment calculations using one-tail Fisher exact test and FDR correction.
#Output: ds(rawP=Float64, fgN=Int, bgN=Int, exN=Float64, fdr=Float64[, explain=String])
#fg, bg and balance have the same meaning as GOenrichFB.
#fg, bg, item and item2 should be vectors with the same type.
#In default, enrich_test normalize by a background weight rather than subsample. This weight can be pre-calculated by balancemask2 and assigned in bg_weight. In case the background is slightly insufficient for minor bins, try balance_kw=(toofewbg=:warn,) to get an approximate result under imperfect balance.
#list_fg_hit: whether output `fg_hit' field for fg hitted item.
#list_fg_max: if fgN is larger than this number, `fg_hit' will just be an empty string.
#list_fg_idconv: `fg_hit' shows the converted id by dtshift(item, list_fg_idconv...; safe=true). e.g., list_fg_idconv=(exn["geneid"], exn["gene_name"], "??"). Function actually do dtshift(item, list_fg_idconv...).
#list_fg_fun: function used to convert fg hit vector to a string.
#Indeed list_fg_idconv and list_fg_fun follows the "row-rule". e.g.:
#          list_fg_idconv=(geneid, (gene_name, cite_num), ("", 0)),
#          list_fg_fun=((t, n),)->join(t[sortri(n, rev=true)], ',')
#will sort the gene names accoriding to the cite_num in T["fig_hit"].
#If no enrichment detected, function will throw a warning (unless nowarn=true) and output an zero-row table.
#See also: GOenrichFB, balancemask, balancemask2
#2 Apr 2021 > 1 Sep 2021

export enrich_test
function enrich_test(fg::AbstractVector{T1}, bg::AbstractVector{T1},
                     item0::AbstractVector{T1}, cata0::AbstractVector{T2};
                     fdr_cutoff::Float64=0.1,
                     fgN_cutoff::Int=2,
                     balance::Union{Tuple{AbstractVector, AbstractVector}, Nothing}=nothing,
                     bg_weight::Union{AbstractVector{Float64}, Nothing}=nothing,
                     balance_method::Symbol=:bg_weight, balance_kw=NamedTuple(),
                     explain::Union{Tuple{AbstractVector{T2}, Group}, Nothing}=nothing,
                     list_fg_hit::Bool=true, list_fg_idconv::Union{Tuple, Nothing}=nothing,
                     list_fg_max::Int=100, list_fg_fun::Function=x->join(sort(x), ','),
                     savefile::AbstractString="",
                     nowarn::Bool=false
                ) where {T1, T2}
    fg=sort(fg)
    bg=sort(bg)
    if !isnothing(balance)
        isnothing(bg_weight) || error("balance and bg_weight cannot be assigned simultaneously.")
        fgv,bgv=balance
        length(fg)==length(fgv) || error("Foreground items and their balance values have inconsistent lengths.")
        length(bg)==length(bgv) || error("Background items and their balance values have inconsistent lengths.")
        #First remove shared items in background.
        l=.!ismbr(bg, fg, Asorted=true, Bsorted=true)
        bg=bg[l]
        bgv=bgv[l]

        #Remove duplicated items.
        l=uniqr(fg, sorted=true)[1]
        fg=fg[l]
        fgv=fgv[l]
        
        l=uniqr(bg, sorted=true)[1]
        bg=bg[l]
        bgv=bgv[l]

        if balance_method==:bg_weight
            _, t, bg_weight=balancemask2(fgv, bgv; balance_kw...)
            println(f"Balance FG($1):BG($2) = 1:$3 (by weight)"(length(fg),length(bg),t))
            tt=bg_weight.>1
            if any(tt)
                @warn(f"$1 ($(2=.2g)%) of $3 background balance weight >1. They are been forcely set to 1."(count(tt), length(tt), 100*mean(tt)))
                bg_weight[tt].=1.0
            end
        elseif balance_method==:bg_subsample
            l, t, _=balancemask2(fgv, bgv; balance_kw...)
            bg=bg[l]
            println(f"Balance FG($1):BG($2) = 1:$3"(length(fg),length(bg),t))
        elseif balance_method==:fg_and_bg_subsample
            l1, l2=balancemask(fgv, bgv; balance_kw...)
            fg=fg[l1]
            bg=bg[l2]
            println(f"Balance FG($1=$(3=.3g)%) vs BG($2=$(4=.3g)%)"(length(fg), length(bg), 100*mean(l1), 100*mean(l2)))
        else
            error("Unknown balance_method $balance_method. It can only be :bg_by_weight, :bg_subsample or :fg_and_bg_subsample.")
        end
    elseif !isnothing(bg_weight)
        balance_method==:bg_weight || error("When bg_weight is given, balance_method must be :bg_weight.")
        length(bg)==length(bg_weight) || error("Background items and bg_weight have inconsistent lengths.")
        #First remove shared items in background.
        l=.!ismbr(bg, fg, Asorted=true, Bsorted=true)
        bg=bg[l]
        bg_weight=bg_weight[l]
        
        #Remove duplicated items.
        l=uniqr(fg, sorted=true)[1]
        fg=fg[l]
        
        l=uniqr(bg, sorted=true)[1]
        bg=bg[l]
        bg_weight=bg_weight[l]
    else
        #First remove shared items in background.
        l=.!ismbr(bg, fg, Asorted=true, Bsorted=true)
        bg=bg[l]

        #Remove duplicated items.
        l=uniqr(fg, sorted=true)[1]
        fg=fg[l]
        
        l=uniqr(bg, sorted=true)[1]
        bg=bg[l]
    end
    ug=unival(vcatr(fg, bg))

    t=sortperm(item0)
    item=item0[t]
    cata=cata0[t]
    
    l=ismbr(item, fg, Asorted=true)
    ll=ismbr(cata, cata[l])
    item=item[ll]
    cata=cata[ll]

    EMPTYOUT=isnothing(explain) ? zerorow(ds(cata=T2, rawP=Float64, fgN=Int, bgN=Int, exN=Float64, fdr=Float64)) : zerorow(ds(cata=T2, rawP=Float64, fgN=Int, bgN=Int, exN=Float64, fdr=Float64, explain=String))
    EMPTYOUT=zerorow(ds(cata=T2, rawP=Float64, fgN=Int, bgN=Int, exN=Float64, fdr=Float64))
    if !isnothing(explain)
        EMPTYOUT["explain"]=zerorow(String)
    end
    if !isnothing(list_fg_hit)
        EMPTYOUT["fg_hit"]=zerorow(String)
    end
    
    if isempty(cata)
        nowarn || @warn("No matched item.")
        return EMPTYOUT
    end
    flen=length(fg)
    if isnothing(bg_weight)
        blen=length(bg)
    else
        bg_weight_sum=sum(bg_weight)
    end
    T=grpfun(cata, item) do citem
        inlist=in(Set(citem))
        fg_hit=count(inlist, fg)
        isbgin=inlist.(bg)
        bg_hit=count(isbgin)
        fisherM, ex_hit_n=if isnothing(bg_weight)
            ([fg_hit flen-fg_hit; bg_hit blen-bg_hit], bg_hit*flen/blen)
        else
            t=sum(bg_weight[isbgin])
            ([fg_hit flen-fg_hit; round(Int, t) round(Int, sum(bg_weight[.!isbgin]))],
             t*flen/bg_weight_sum)
        end
        p=fishertest(fisherM, alternative="greater")
        (p, fg_hit, bg_hit, ex_hit_n)
    end |> tb((c"rawP, fgN, bgN, exN", "cata"))

    T=rec(T, T["fgN"].>=fgN_cutoff)
    if rnum(T)>0
        T["fdr"]=padjust(T["rawP"], "fdr")
        T=rec(T, T["fdr"].<=fdr_cutoff)
    end
    if rnum(T)==0
        nowarn || @warn("No enrichment passed cutoff.")
        return EMPTYOUT
    end
    T=rec(T, sortri((T["rawP"], T["fgN"])))
    if !isnothing(explain)
        T["explain"]=dtshift(T["cata"], explain[1], explain[2], nothing, safe=true)
    end
    
    if list_fg_hit
        idconv=if isnothing(list_fg_idconv)
            item
        else
            dtshift(item, list_fg_idconv..., safe=true)
        end
        T["fg_hit"]=grpfun(fastgrp(cata, T["cata"]), item, idconv, no=true) do Ti, citem, cidconv
            if T["fgN"][Ti]>list_fg_max
                ""
            else
                list_fg_fun(getrow(cidconv, ismbr(citem, fg, Asorted=true, Bsorted=true)))::AbstractString
            end
        end
    end
    
    if !isempty(savefile)
        hd=c"cata, fgN, bgN, exN, rawP, fdr"
        if !isnothing(explain)
            hd=[hd; "explain"]
        end
        if list_fg_hit
            hd=[hd; "fg_hit"]
        end
        writetb(savefile, T, hd)
    end
    T
end
#}}

#{{ GOenrichFB GOenrichGrp
#[[ GOenrichFB GOenrichGrp ]]
#Table|nothing = GOenrichFB(foreground, background; balance=(fg,bg), onlybalancebg=true, ann="org.Hs.eg.db", ID="ensembl", outfile="filename.tsv", sigP=0.05|0.01, topGO_algorithm=:classic, padj_method=:fdr|:none)
#Table|nothing, GOdata, GOtest = GOenrichFB(...; outRobj=true) #Output topR objects also.
#Table|nothing = GOenrichGrp(gene, group; balance=(fg,[ bg=fg]), bg=gene, ann="org.Hs.eg.db", ID="ensembl", outfile="filename.tsv", sigFDR=0.05)
#Check the GO enrichment using Fisher test in topGO. If outfile is given, a text file will be saved. If any GO terms are enriched, a table will be return, otherwise will return nothing instead.
#If the background is need to be balanced, set balance value( e.g. as gene expression). When onlybalancebg is true, function will use balancemask2(), i.e., only subsample background but not subsample forground. In the case that background is not enough, an error is occur. Set onlybalancebg=false will let function use balancemask() instead and will not trigger such error. The shared item in background, and the duplicated item will be remove, so does their balance values, before apply to balancemask function.
#The supported ID could be "entrez", "genbank", "alias", "ensembl", "symbol", "genename".
#padj_method: the method of p.adjust in R. Can be :holm, :hochberg, :hommel, :bonferroni, :BH, :BY, :fdr, :none. The default is :fdr if method is :classic, or :none otherwise.
#sigP: The default is 0.05 if p-value adjustment is used, or 0.01 otherwise.
#topGO_algorithm: can be :classic|:weight01|...
#* (i) 'elim' : this method processes the GO terms by traversing the GO hierarchy from bottom to top, ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later to more general (higher) GO terms. When it assesses a higher (more general) GO term, it discards any genes that are annotated with significantly enriched descendant GO terms (considered significant using a pre-defined P-value threshold). This method does tend to miss some true positives at higher (more general) levels of the GO hierarchy.
#* (ii) 'weight': this is designed to miss less true positives than the 'elim' method. The significance scores of connected nodes (a parent and its child in the GO hierarchy) are compared in order to detect the locally most significant terms in the GO hierarchy. 
#* (iii) 'classic': each GO term is tested independently, not taking the GO hierarchy into account.
#* (iv) 'weight01': this is the default method used by TopGO, and is a mixture of the 'elim' and 'weight' methods.
#* (v) 'parentchild': described by Grossman et al (2007) : when assessing a GO term, it takes into accoount the annotation of terms to the current term's parents, and so reduces false positives due to the inheritance problem. You can also use this algorithm in a very nice tool called Ontologizer , which has a lovely GUI (I used it in a collaboration on gene duplicates with Shane Woods, described in his paper Woods et al (2013)). 
# The 'elim' and 'weight' methods are designed to be more conservative (have less false positives) than the 'classic' method, so should have larger p-values in general. Simulation results reported by Alexa et al (2006) show that the 'weight' algorithm has less false positives than 'classic' and misses few true positives, while 'elim' has even less false positives than 'weight' but misses more true positives.  Grossman et al (2007) compare the 'elim', 'weight' and 'parentchild' algorithms, and say that each has advantages in certain scenarios. 
#See also: balancemask, balancemask2, GSVA
#Xiong Jieyi, October 27, 2014>Aug 28, 2015>Oct 2, 2015>23 Apr 2017>16 Oct 2019

export GOenrichFB,GOenrichGrp
function _GOenrichCore(val::BitVector,genes::AbstractVector{T},ann::AbstractString,ID::AbstractString; topGO_algorithm::Symbol=:classic, padj_method::Symbol=topGO_algorithm==:classic ? (:fdr) : (:none), sigP::Real=padj_method==:none ? 0.01 : 0.05) where {T<:AbstractString}
    myr=importpkg(:myR)
    @pkgfun(myR, eR, callr, relem, j2r, r2j, callrj, callr1)
    rlt=rowpile()
    GOdatas=Dict{String,myr.RObject}()
    GOtests=Dict{String,myr.RObject}()
    ontology_ls=c"BP,MF,CC"
    for ontology in ontology_ls
        println("*** TEST $ontology ONTOLOGY ***")
        GOdata=callr("new","topGOdata",
                     ontology=ontology,
                     allGenes=j2r(int(val),names=genes),
                     geneSel=eR("function(x){x==1}"),
                     annot=eR("annFUN.org"), mapping=ann, ID=ID,
                     nodeSize=5)
        GOdatas[ontology]=GOdata
        GOtest=callr("runTest",GOdata,algorithm=topGO_algorithm,statistic="fisher")
        GOtests[ontology]=GOtest
        allGOlen=callr("length", callr("usedGO",object = GOdata))
        GenTb=callr("GenTable",GOdata,classicFisher = GOtest, topNodes = allGOlen, numChar=9999)
        FDR=callr("p.adjust",relem(GenTb,"classicFisher"),method=padj_method)

        #For unknown reason, classicFisher and Term occationally become invalid StringVector after transferred to Julia DataFrame.
        GenTb=callr("function(y,z,w){y['classicFisher']<-as.vector(data.matrix(y['classicFisher']));subset(cbind(y, Padj_$padj_method=z), Padj_$padj_method<=w)}",GenTb,FDR,sigP)
        # end

        N=callrj("nrow", GenTb)
        if N>0
            println("$N GO item(s) is/are enriched.")
            addrows!(rlt, tb(r2j(GenTb), Ontology=fill(ontology, N)))
         else
            minFDR=callr1("min",FDR)
            println("No enriched item was found. The lowest p_adj($padj_method) is (p_adj=$minFDR): ")
        end
    end
    return (value(rlt), GOdatas, GOtests)
end
function GOenrichFB(fg::AbstractVector{T},bg::AbstractVector{T};
                    ann::AbstractString="org.Hs.eg.db",
                    balance::Union{Tuple{Vector,Vector}, Nothing}=nothing,
                    onlybalancebg::Bool=true,
                    ID::AbstractString="ensembl",
                    outfile::AbstractString="",
                    outRobj::Bool=false,
                    wargs...) where {T<:AbstractString}
    if isnothing(balance)
        #First remove shared items in background.
        l=.!ismbr(bg, fg, Asorted=true, Bsorted=true)
        bg=bg[l]

        #Remove duplicated items.
        l=uniqr(fg, sorted=true)[1]
        fg=fg[l]
        
        l=uniqr(bg, sorted=true)[1]
        bg=bg[l]
    else
        fgv,bgv=balance
        @assert(length(fg)==length(fgv),"Foreground items and their balance values have inconsistent lengths.")
        @assert(length(bg)==length(bgv),"Background items and their balance values have inconsistent lengths.")
        #First remove shared items in background.
        l=.!ismbr(bg,fg)
        bg=bg[l]
        bgv=bgv[l]

        #Remove duplicated items.
        l=uniqr(fg)[1]
        fg=fg[l]
        fgv=fgv[l]
        
        l=uniqr(bg)[1]
        bg=bg[l]
        bgv=bgv[l]

        if onlybalancebg
            l,t=balancemask2(fgv,bgv)
            bg=bg[l]
            println(f"Balance FG($1):BG($2) = 1:$3"(length(fg),length(bg),t))
        else
            l1, l2=balancemask(fgv, bgv)
            fg=fg[l1]
            bg=bg[l2]
            println(f"Balance FG($1=$(3=.3g)%) vs BG($2=$(4=.3g)%)"(length(fg), length(bg), 100*mean(l1), 100*mean(l2)))
        end
    end
    @pkgfun(myR,eR,callr)
    eR("library(topGO)")
    callr("library",ann)
    ug=unival(vcatr(fg,bg))
    rlt, GOdata, GOtest=_GOenrichCore(ismemberr(ug,fg),ug,ann,ID;wargs...)
    
    if !isempty(outfile)
        if isnothing(rlt)
            writeall(outfile, "No any enriched GO item was found.")
        else
            # 4 Dec 2019: In weight01 method, no Padj_fdr outputed.
            lbs=c"Ontology, GO.ID, Annotated, Significant, Expected, classicFisher, Padj_fdr, Term"
            rltkeys=collect(keys(rlt))
            lbs=lbs[ismbr(lbs, rltkeys)]
            lbs=[lbs[1:end-1]; setdiff(rltkeys, lbs); lbs[end]]
            writetb(outfile, rlt, lbs)
        end
    end
    if outRobj
        return (rlt, GOdata, GOtest)
    else
        return rlt
    end
end
function GOenrichGrp(gene::Vector{T},grp::Group;
                     bg::Vector{T}=gene,outfile::AbstractString="",
                     balance::Union{Tuple{Vector,Vector}, Nothing}=nothing,
                     wargs...) where {T<:AbstractString}
    if !isnothing(balance)
        if length(balance)==1
            allfgv=balance[1]
            bgv=balance[1]
        else
            allfgv,bgv=balance
        end
        rlt,_=grpfunwith(grp,gene,allfgv;input_grpid=true)do grp,fg,fgv
            println("--- Group $grp ---")
            cgo=GOenrichFB(fg,bg;balance=(fgv,bgv),wargs...)
            cgo["Group"]=fill(grp,size(cgo,1))
            cgo
        end
    else
        rlt,_=grpfunwith(grp,gene;input_grpid=true)do grp,fg
            println("--- Group $grp ---")
            cgo=GOenrichFB(fg,bg;wargs...)
            cgo["Group"]=fill(grp,size(cgo,1))
            cgo
        end
    end
    if !isempty(outfile)
        if isnothing(rlt)
            writeall(outfile, "No any enriched GO item was found.")
        else
            writetb(outfile, rlt, c"Group, Ontology, GO.ID, Annotated, Significant, Expected, classicFisher, Padj_fdr, Term")
        end
    end
    return rlt
end

#}}

#{{ parseDemulReport
#[[ parseDemulReport ]]
#T = parseDemulReport(demultiplex_path)
#Parse demultiplex_path/Reports and get information of %>=Q30, MeanQual and PASSchastity. If a sample were sequenced in multiple lanes, this funciton will output the mean value.
#Output fields: lane,PF_clusters,%ofTheLane,%PrefectBarcode,%OneMismatchBarcode,Yield_Mbases,%PF_clusters,%>=Q30bases,meanQualityScore
#See also: parseGTF
#Xiong Jieyi, Feb 5, 2016

export parseDemulReport
function parseDemulReport(reportpth)
    pth=getfiles(joinpath(reportpth,"Reports/html/"),typ=:dir)[1][1]
    t1,t2=getfiles(pth,typ=:dir)
    pth=t1[findonly(!ismemberr(t2,c"all,default"))]
    t1,t2=getfiles(pth,typ=:dir)
    l=t2.!="all"
    smp=tb(ID=t2[l],
           lane_report_html = map(f"$1/all/lane.html",t1[l]))

    t=rowfunwith(smp) do csmp
        txt=readchomp(csmp["lane_report_html"])
        rpl=rowpile()
        pso=ParseStr(txt)
        next(pso,"Lane Summary","Score</th>","</tr>")
        while next(pso,(:check,"<tr>"))[1]
            C=tb(c"lane,PF_clusters,%ofTheLane,%PrefectBarcode,%OneMismatchBarcode,Yield_Mbases,%PF_clusters,%>=Q30bases,meanQualityScore",
                 next(pso,
                     "<td>",("</td>",int), #Lane
                     # "<td>",("</td>",ascii), #Barcode
                     "<td>",("</td>",x->int(replace(x, "," => ""))), #PF Clusters
                     "<td>",("</td>",x->float(x)/100), # % of the lane
                     "<td>",("</td>",x->float(x)/100), # % Perfect barcode
                     "<td>",("</td>",x->float(x)/100), # % One mismatch barcode
                     "<td>",("</td>",x->int(replace(x, "," => ""))), #Yield (Mbases)
                     "<td>",("</td>",x->float(x)/100), #% PF Clusters
                     "<td>",("</td>",x->float(x)/100), #% >= Q30 bases
                     "<td>",("</td>",float), #Mean Quality Score
                     "</tr>"))
            C["sample"]=csmp["ID"]
            addrow!(rpl,C)
        end
        value(rpl)
    end
    t[1]
end
#}}

#{{ DESeq DESeqMF
#[[ DESeq DESeqMF ]]
# Table = DESeq(A_counts::Matrix|Vector, B_counts::Matrix|Vector; norep=false, sizeFactors=(vecA, vecB))
#     --Run DESeq to test different gene expression.
# (FDRs, raw_p-values, log2_coef::Dict) = DESeqMF( dataMatrix, keyFactor, otherFactor1, otherFactor2, ...; sizeFactors=vec)
# (FDRs, rawP, log2_coef) = DESeqMF( dataMatrix, keyFactor::BoolVec, ...; ...) #log2coef: true - false
#     --Run DESeq to test if the key factor model is better than non-key-factor model, with other factors under consideration.
# If no any repeats, set norep=true to use "blind" method in estimateDispersions.
# Each row of counts is a gene, while each column of counts is a sample. nbinomTest will do between A and B.
# DESeq output table contains fields:
# "id"            
# "baseMean"      
# "baseMeanA"     
# "baseMeanB"     
# "foldChange"    
# "log2FoldChange"
# "pval"          
# "padj"
# In the outpus, the positive direction of fold change is B-A in DESeq().
# For the log2_coef output of DESeqMF(), the alphabetical smallest key factor will not be given.
#
#See also: edgeR_anova, GOenrichFB, GOenrichGrp, myR, GSVA, DESeq2, DESeq2MF
#Xiong Jieyi, Mar 10, 2016 > 23 Jul 2017 >24 Jul 2017 >31 May 2018 >16 Aug 2018

export DESeq, DESeqMF
function DESeq(A::AbstractArray{T1}, B::AbstractArray{T2}; norep::Bool=false, sizeFactors::Tuple=()) where {T1<:Integer, T2<:Integer}
    @pkgfun(myR,callr,r2j,eR,rlibrary)
    X=hcat(A,B)
    L=[fill(:A,size(A,2));fill(:B,size(B,2))]
    rlibrary("DESeq")
    cds=callr("newCountDataSet",X,L)
    if isempty(sizeFactors)
        cds=callr("estimateSizeFactors", cds)
    else
        # callr("`sizeFactors<-`", cds, vcat(sizeFactors...))
        cds=callr("function(x, y){sizeFactors(x)<-y;x}", cds, vcat(sizeFactors...)|>x->x./mean(x))
    end
    if norep
        cds=callr("estimateDispersions",cds,method="blind",sharingMode="fit-only")
    else
        cds=callr("estimateDispersions",cds)
    end
    r2j(callr("nbinomTest",cds,"A","B"), NA=NaN, keepvec=true)
end
function DESeqMF(M::AbstractMatrix{T1}, keyF::AbstractVector{T2}, Fs::AbstractVector...; sizeFactors::AbstractVector=[]) where {T1<:Integer, T2<:Union{String, Symbol}}
    @pkgfun(myR,callr,callrj,eR,r2j,rnames,relemj,rlibrary)
    if isempty(Fs)
        fml1=eR("count~keyF")
        fml0=eR("count~1")
        L=callr("data.frame", keyF=keyF)
    else
        Flb=f"F$1".(1:length(Fs))
        t=join(Flb,"+")
        fml1=eR("count~$t+keyF")
        fml0=eR("count~$t")
        L=callr("data.frame", tb(tb(Flb,Fs),keyF=keyF))
    end
    
    rlibrary("DESeq")
    cds=callr("newCountDataSet",M,L)
    if isempty(sizeFactors)
        cds=callr("estimateSizeFactors",cds)
    else
        # callr("`sizeFactors<-`", cds, sizeFactors)
        cds=callr("function(x, y){sizeFactors(x)<-y;x}", cds, sizeFactors./mean(sizeFactors))
    end
    cds=callr("estimateDispersions",cds, method="pooled-CR", modelFormula=fml1)
    fit1=callr("fitNbinomGLMs", cds, fml1)
    fit0=callr("fitNbinomGLMs", cds, fml0)
    pval=callr("nbinomGLMTest",fit1,fit0)
    kys=filter(r"^keyF", rnames(fit1))
    T=rename!(relemj(fit1, kys, NA=NaN, keepvec=true), r"^keyF", "")
    (r2j(callr("p.adjust",pval,method="fdr"), keepvec=true), r2j(pval, keepvec=true), T)
end
function DESeqMF(M::AbstractMatrix{T1}, keyF::AbstractVector{Bool}, arg...; kw...) where {T1<:Integer}
    fdr, rawP, T=DESeqMF(M, ifelse.(keyF, "B", "A"), arg...; kw...)
    (fdr, rawP, T["B"])
end
#}}

#{{ DESeq2 DESeq2LRT
#[[ DESeq2 DESeq2LRT ]]
# Table = DESeq2(A_counts::Matrix|Vector, B_counts::Matrix|Vector)
# Or ...= DESeq2(counts::Matrix; isB=Bool_vec)
#     --Run DESeq2 to test different gene expression using Wald test.
# Table = DESeq2LRT(counts::Matrix, keyFactor[, otherFactor1, otherFactor2, ...]; levels=...)
#     --Run DESeq2 in likelihood ratio test.
# Function outputs below fields, all in Vector{Float64}. The positive direction of fold change is B-A.
#   baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
# In DESeq2LRT, keyFactor and otherFactors can be either categroicial (string or symbol) or continuous (numeric). For categroical factor, the order of keyFactor is assigned by `levels=...'. Assigning orders for otherFactors are not supported (and it is not necessary).
# See also: DESeq, DESeqMF, edgeR_anova
# 23 Feb 2022 > 13 Dec 2023 > 28 Feb 2024

export DESeq2, DESeq2LRT
function DESeq2(A::AbstractArray{<:Integer}, B::AbstractArray{<:Integer})
    #Based on the guides: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
    @pkgfun(myR, callrj, rlibrary)
    X=hcat(A, B)
    L=[fill(:A, size(A, 2)); fill(:B, size(B, 2))]
    rlibrary("DESeq2")
    callrj("""function(X, L){
dds <- DESeqDataSetFromMatrix(countData=X, colData=data.frame(condition=L), design=~condition)
dds <- DESeq(dds)
res <- results(dds)
as.list(res)
}""", X, L)
end
function DESeq2(M::AbstractMatrix{<:Integer}; isB::Union{Nothing, AbstractVector{Bool}}=nothing)
    @assert(!isnothing(isB), "Missed the 2nd dataset or missed isB=....")
    DESeq2(M[:, .!isB], M[:, isB])
end
function DESeq2LRT(M::AbstractMatrix{<:Integer}, keyF::AbstractVector{T}, Fs::AbstractVector...;
                   levels::Union{Nothing, AbstractVector{T}}=nothing) where T<:Union{String, Symbol, Real}
    @assert(size(M, 2)==length(keyF), "Factor length should equal to the column number of matrix.")
    @pkgfun(myR, callrj, rlibrary, callr, eR, j2r_na)
    factorKeyF=if isnothing(levels)
        if T<:Real
            j2r_na
        else
            x->callr("factor", x)
        end
    else
        T<:Real && error("For fool-proofing, key factor cannot be numeric when levels is given.")
        x->callr("factor", x, levels=levels)
    end
    factorF=x->isa(x, Real) ? j2r_na : callr("factor", x)
    if isempty(Fs)
        fml1=eR("~keyF")
        fml0=eR("~1")
        L=callr("data.frame", keyF=factorKeyF(keyF))
    else
        Flb=f"F$1".(1:length(Fs))
        t=join(Flb,"+")
        fml1=eR("~$t+keyF")
        fml0=eR("~$t")
        L=callr("data.frame", tb(ds(Flb, map(factorF, Fs)), keyF=factorKeyF(keyF)))
    end
    rlibrary("DESeq2")
    callrj("""function(M, L, fml1, fml0){
dds <- DESeqDataSetFromMatrix(countData=M, colData=L, design=fml1)
dds <- DESeq(dds, test="LRT", reduced=fml0)
res <- results(dds)
as.list(res)
}""", M, L, fml1, fml0)
end
#}}

#{{ edgeR_anova
#[[ edgeR_anova ]]
#Syntax: Table, comparison_group = edgeR_anova(Readnum_Matrix, Group; lib!size=vec, norm!factors=vec, remove!zeros=false)
# Table = edgeR_anova(Readnum_Matrix, Group.==comparison_group; ...)
#Run edgeR ANOVA mode to test different gene expression. 
#Each row is a gene and each column is a sample. Group length should equal to matrix column number. Group could be a vector of string, symbol or integer.
# lib!size: numeric vector giving the total count (sequence depth) for each library.
# norm!factors: numeric vector of normalization factors that modify the library sizes.
# remove!zeros: logical, whether to remove rows that have 0 total count.
# See ?DGEList in edgeR.
#Returns fields:  "logFC[_group...]", "LR", "PValue", "FDR"
# Updates in 5 Aug 2019:
# The 2nd output "comparison_group" is the "comparison" field of edgeR result. When Group has only two factors, it is the "higher" factor, in String or in Symbol depending on the input, showing the direction of logFCis this field - other filed. (I have no idea about this factor so far when Group has three or more factors.)
#See also: DESeq, DESeqMF, GOenrichFB, GOenrichGrp, GSVA
#Xiong Jieyi, Jun 7, 2016 > 2 Jun 2018 > 5 Aug 2019

export edgeR_anova
function edgeR_anova(M::AbstractMatrix{T1}, G::AbstractVector{T2}; wargs...) where {T1<:Real,T2<:Union{AbstractString,Symbol}}
    rnum(G)==size(M,2) || error("Group length should equal to matrix column number.")
    @pkgfun(myR, callr, eR, relem, relem1, j2r, r2j, rlibrary)
    rlibrary("edgeR")
    y=callr(:DGEList,counts=j2r(M,row!names=1:size(M,1)),
            group=callr(:factor,G); wargs...)
    design=callr(:model!matrix,eR("~group"),data=relem(y,"samples"))
    y=callr(:estimateGLMCommonDisp,y,design)
    y=callr(:estimateGLMTrendedDisp,y,design)
    y=callr(:estimateGLMTagwiseDisp,y,design)
    fit=callr(:glmFit,y,design)
    lrt=callr(:glmLRT,fit,coef=2:uninum(G))
    # tag=callr(:topTags,lrt,n=size(M,1))["table"]
    toptags=callr(:topTags,lrt,n=size(M,1))
    tag=relem(toptags, "table")
    T=r2j(tag, NA=NaN, keepvec=true) #Changed on 6 Nov 2018
    
    rowidx=int(r2j(callr(:row!names,tag), keepvec=true))

    compgrp=relem1(toptags, "comparison") #Now I have no idea what compgrp will be when facor number >=3. 5 Aug 2019
    if isa(compgrp, String)
        compgrp=replace(compgrp, r"^group"=>"")
        if eltype(G)<:Symbol
            compgrp=Symbol(compgrp)
        end
    end
    (dtshift(1:size(M,1),rowidx,T), compgrp)
end
function edgeR_anova(M::AbstractMatrix{T1}, G::AbstractVector{Bool}; kw...) where {T1<:Real}
    T, cpg=edgeR_anova(M, ifelse.(G, "B", "A"); kw...)
    @assert cpg=="B"
    T
end
#}}

#{{ GSVA
#[[ GSVA ]]
# sig1, sig2, ... = GSVA(Matrix, geneid, gene_list1[, gene_list2, ...]; counts=false, pn=1, method="gsva|ssgsea|zscore|plage")
# ... = GSVA(Matrix, Bool_vec1=ismbr(geneid, gene_list)[, Bool_vec2...]; ... )
# Calculate siguature value using R-GSVA pacakge.
# The matrix of first input: gene for row, sample for column. The output is one or more vectors with the length of sample number. If the input is read numbers rather than expression values, set counts=true.
# pn is the parallel computation core number.
# method="ssgsea" will be much faster than the default method "gsva".
# See the manual on https://www.rdocumentation.org/packages/GSVA/versions/1.20.0/topics/gsva
# `method`: Method to employ in the estimation of gene-set enrichment scores per sample. By default this is set to `gsva` (Hänzelmann et al, 2013) and other options are `ssgsea` (Barbie et al, 2009), `zscore` (Lee et al, 2008) or `plage` (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores over the samples and, in the case of zscore, it combines them together as their sum divided by the square-root of the size of the gene set, while in the case of plage they are used to calculate the singular value decomposition (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector as pathway activity profile.
# See also: DESeq, edgeR_anova, GOenrichFB, GOenrichGrp
# Jieyi Xiong, 13 Sep 2018 > 27 Feb 2019

export GSVA
function GSVA(count::Matrix, geneid::Vector{String}, gene_list1::Vector{String}, gene_lists::Vector{String}...; counts::Bool=false, pn::Int=1, kwargs...)
    if counts && !(eltype(count)<:Integer)
        error("Counts should be integers when counts=true.")
    end
    if eltype(count)<:Integer
        counts || @warn("Input count matrix is integer, but counts=false.")
    else
        counts && @warn("Input count matrix is not integer, but counts=true.")
    end
    @pkgfun(myR, rlibrary, callr, j2r, r2j)
    rlibrary("GSVA")
    if pn>1
        rlibrary("snow")
    end
    M=r2j(callr("gsva", j2r(count, rownames=geneid), callr("list", gene_list1, gene_lists...), parallel!sz=pn, kcdf=counts ? "Poisson" : "Gaussian", kwargs...), NA=NaN, keepvec=true)
    if isempty(gene_lists)
        vec(M)
    else
        tuple(eachr(M, as_vec=true)...)
    end
end
function GSVA(count::Matrix, isfg1::AbstractVector{Bool}, isfgs::AbstractVector{Bool}...; wargs...)
    geneid=string.(1:length(isfg1))
    gene_lists=map(x->geneid[x], [isfg1, isfgs...])
    GSVA(count, geneid, gene_lists...; wargs...)
end
#}}

#{{ GEOdownload
#[[ GEOdownload ]]
#smp_info = GEOdownload(c"GSMXXXXX, GSMXXXXX, ..."; label=GSM_IDs, path="_GEO_downloaded", trymax=10, paral=T|F, directly=false)
# ... = GEOdownload(c"SRRxxxxx, SRRxxxxx"; isSRR=true, ...) #Download by SRR ID.
#Download files from GEO and convert to fastq files. Multi-threads is supported. One GSM-id could links to multiple SRR files. GEOdownload will not merge these SRR fastq files since in some cases, their technical parameters (like whether pair-ended) are not identical.
#The finial files will be $path/$(label).$(SSR).SE|PE1|PE2.fastq.gz. Besides, under $path, SRA_file/ is the sra files; fastq/ is the fastq.gz files before merging; GEO_sample_info.jld is the sample info table the same as function output.
#If directly=true, function will download and convert fastq file using fastq-dump directly. Otherwise, a .sra file will be download by wget first, then converted to fastq using fastq-dump.
#If FTP download is failed, function will retry until get trymax times.
#If path exists and some finished fastq files are in it, these files will not be re-downloaded.
#Entrezdirect (ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/) and SRA Toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) need be installed.
#See also:
#Xiong Jieyi, 16 May 2017 > 7 Jun 2018 > 30 Oct 2019 > 13 May 2020

export GEOdownload
function GEOdownload(GSMs::Vector{String}; label::Vector{String}=GSMs, path::String="_GEO_downloaded", trymax::Int=10, paral::Bool=nprocs()>1, directly::Bool=false, isSRR::Bool=false)
    
    if !(isuniqr(GSMs) && isuniqr(label))
        error("GSM ID or labels is not unique.")
    end

    if isSRR
        smp=tb(SRR=GSMs, label=label)
        if isdir(path)
            println("Path $path has already exists. Check unfinished file...")
            smp1=readtb(joinpath(path, "GEO_sample_info.tsv"))
            if Set(eachr(d"SRR, label"smp))!=Set(eachr(d"SRR, label"smp1))
                error("Sample information doesn't match. Stop downloading.")
            end
            l=isfile.(f"$path/$1.$2.PE1.fastq.gz".(smp["label"], smp["SRR"])) .| isfile.(f"$path/$1.$2.SE.fastq.gz".(smp["label"], smp["SRR"]))
            if all(l)
                println("All files are found. Do nothing.")
                return zerorow(smp)
            else
                println(f"$1 finished fastq file(s) has(ve) found. Skip them."(sum(l)))
                smp=rec(smp, .!l)
            end
        else
            mkdir(path)
            writetb(joinpath(path, "GEO_sample_info.tsv"), smp)
        end
    else
        usmp=tb(GSM=GSMs, label=label)
        if isdir(path)
            println("Path $path has already exists. Check unfinished file...")
            smp=readtb(joinpath(path, "GEO_sample_info.tsv"))
            if Set(eachr(d"GSM, label"smp))!=Set(eachr(d"GSM, label"usmp))
                error("Sample information doesn't match. Stop downloading.")
            end
            l=isfile.(f"$path/$1.$2.PE1.fastq.gz".(smp["label"], smp["SRR"])) .| isfile.(f"$path/$1.$2.SE.fastq.gz".(smp["label"], smp["SRR"]))
            if all(l)
                println("All files are found. Do nothing.")
                return zerorow(smp)
            else
                println(f"$1 finished fastq file(s) has(ve) found. Skip them."(sum(l)))
                smp=rec(smp, .!l)
            end
        else
            println("Fetching sample info")
            smp, _=rowfunwith(usmp) do T
                print(T["GSM"])
                S=sfc"esearch -db sra -query \"$1\"| efetch -format docsum"(T["GSM"])|>readchomp
                isempty(S) && error("GEO ID is not found: "*T["GSM"])
                O=ParseStr(S)
                T["SRX"]=next(O,"<Experiment acc=\"",("\"",),(:no,"SRX"))[1] #SRX is no more used.
                # SRP=next(O,"<Study acc=\"",("\"",),(:no,"SRP"))[1]
                SRRs=String[]
                while next(O, (:test,"<Run acc=\""))[1]
                    push!(SRRs, next(O,("\"",))[1])
                end
                # (SRX,join(SRRs,','))
                println("\t"*join(SRRs, ","))
                tb(reprow(T, length(SRRs)), SRR=SRRs)
            end
            mkdir(path)
            writetb(joinpath(path, "GEO_sample_info.tsv"), smp)
        end
    end
    
    cd(path) do
        # jsave("GEO_sample_info", asmp)
        # for (fi,(lb,SRX,t)) in enumerate(eachrow((smp["label"],smp["SRX"],smp["SRRs"])))
        mapfun=paral ? pmapwd : map
        mapfun(1:rnum(smp), smp["label"], smp["SRR"]) do fi, lb, SRR
            println("$lb $SRR -- Downloading")
            if directly
                for tryi=1:trymax
                    try
                        sf"fastq-dump --outdir _fastq --gzip --split-files $1"(SRR)
                        break
                    catch err
                        if tryi==trymax
                            rethrow(err)
                        else
                            println(f"Download failed $1 times. Wait $2 second and try again...."(tryi, tryi*5))
                            sleep(tryi*5)
                        end
                    end
                end
            else
                for tryi=1:trymax
                    try
                        sf"prefetch -O _SRA_file $1"(SRR)
                        break
                    catch err
                        if tryi==trymax
                            rethrow(err)
                        else
                            println(f"Download failed $1 times. Wait $2 second and try again...."(tryi, tryi*5))
                            sleep(tryi*5)
                        end
                    end
                end
                println("Converting to Fastq")
                sf"fastq-dump --outdir _fastq --gzip --skip-technical  --readids --dumpbase --split-files --clip _SRA_file/$1/$1.sra"(SRR)
            end
            # if multiSRR
            #     #=
            #     13 May 2020
            #     In some rare cases, samples are mixed with PE and SE data.
            #     =#
            #     println("Merging")
            #     s1=f"_fastq/$(1)_1.fastq.gz".(SRRs)
            #     @assert all(isfile, s1)
            #     l=isfile.(f"_fastq/$(1)_2.fastq.gz".(SRRs))
            #     if any(l)
            #         s2=f"_fastq/$(1)_2.fastq.gz".(SRRs[l])
            #         sf"zcat $1|gzip -c > $(2)_PE2.fastq.gz"(join(s2[l], ' '), lb[l])
            #         sf"zcat $1|gzip -c > $(2)_PE1.fastq.gz"(join(s1[l], ' '), lb[l])
            #         rm.(s2)
            #     end
            #     if any(.!l)
            #         sf"zcat $1|gzip -c > $(2)_SE.fastq.gz"(join(s1[.!l], ' '), lb[.!l])
            #     end
            #     rm.(s1)
            # else
            if isfile(f"_fastq/$(1)_2.fastq.gz"(SRR))
                sf"mv '_fastq/$(1)_2.fastq.gz' '$(2).$(1).PE2.fastq.gz'"(SRR, lb)
                sf"mv '_fastq/$(1)_1.fastq.gz' '$(2).$(1).PE1.fastq.gz'"(SRR, lb)
            else
                sf"mv '_fastq/$(1)_1.fastq.gz' '$(2).$(1).SE.fastq.gz'"(SRR, lb)
            end
            # end
            if !directly
                rm("_SRA_file/$SRR/$SRR.sra")
                rm("_SRA_file/$SRR")
            end
            println("Done.")
        end
        directly || sh"rmdir _SRA_file"
        sh"rmdir _fastq"

        ## Check PE or SE
        T=readtb("GEO_sample_info.tsv")
        if !haskey(T, "PEorSE")
            T["PEorSE"]=map(T["label"], T["SRR"]) do lb, SRR
                t=matchfile("$lb.$SRR.*.fastq.gz", onlyone=false, hasone=true)
                c"SE, PE"[length(t)]
            end
            writetb("GEO_sample_info.tsv", T)
        end
    end
    println("Finished.")

    smp
end
#}}

#{{ survdiff, survcoxph
#[[ survdiff survcoxph ]]
# pvalue = survdiff(surv_time, status, group; vb=0)
# pvalue, coef, Robj = survcoxph(surv_time, status, mainFactor[, otherFactor1, otherfactor2, ...]; vb=0)
# Test survival differences among catalogus or to a numeric variable.
# `status` could be a vector of 0/1 or 1/2 to represent censored vs “dead”.
# verbose: Whether show model detail on screen.
# When mainFactor is numeric, the coefficients in a Cox regression relate to hazard; a positive coefficient indicates a worse prognosis and a negative coefficient indicates a protective effect of the variable with which it is associated.
# When mainFactor is factorial, the coeff is a Pair vector for each factor => coefficient.
# If calculation is failed, function return NaN.
# It needs survival package in R (see http://bioconnector.org/workshops/r-survival.html).
# (beta) The 3rd output of survcoxph() is the R object of coxph(...). 3 Aug 2022
# See also: survivalplot, survfit
#Xiong Jieyi, 15 Jun 2018

export survdiff, survcoxph
function survdiff(time::AbstractVector{Tt}, status::AbstractVector{Ts}, grp::AbstractVector; vb::Int=0) where {Tt<:Real, Ts<:Integer}
    0 < length(time) == length(status) == length(grp) || error("Inputs have different lengths or zero lengths.")
    @pkgfun(myR, rlibrary, callr, r2j)
    rlibrary("survival")
    rfun=k"""function(x, y, z, vb){
        sd<-survdiff(Surv(x, y)~z)
        if(vb==1 || vb==3){show(sd)}
        if(vb==2 || vb==3){show(summary(sd))}
        pchisq(sd$chisq, length(sd$n)-1, lower.tail=FALSE)
        }    
        """
    r2j(callr(rfun, time, status, grp, vb), keepvec=false, NA=NaN)
end
function survcoxph(time::AbstractVector{Tt}, status::AbstractVector{Ts}, X::AbstractVector, Fs::AbstractVector...; vb::Int=0) where {Tt<:Real, Ts<:Integer}
    0 < length(time) == length(status) == length(X) || error("Inputs have different lengths or zero lengths.")
    @pkgfun(myR, rlibrary, callr, r2j)
    rlibrary("survival")
    if isempty(Fs)
        f1=""
        f2=""
        f3="1"
    else
        f1=join(f", F$1".(1:length(Fs)))
        f2=join(f"F$1+".(1:length(Fs)))
        f3=join(f"F$1".(1:length(Fs)), "+")
    end
    rfun="""function(x, y, z, vb$(f1)){
        m1<-coxph(Surv(x, y)~$(f2)z)
        m0<-coxph(Surv(x, y)~$(f3))
        anv<-anova(m1, m0)
        if(vb==1 || vb==3){show(m1); show(anv)}
        if(vb==2 || vb==3){show(summary(m1)); show(summary(anv))}
        list(tail(anv[["Pr(>|Chi|)"]], n=1), m1\$coefficients, names(m1\$coefficients), m1)
        }    
        """
    C=callr(rfun, time, status, X, vb, Fs...)
    pval=r2j(@il(C[1]), keepvec=false, NA=NaN)
    coef, coefname =( r2j(@il(C[2]), keepvec=true), r2j(@il(C[3]), keepvec=true) )
    l=startswith.(coefname, "z")
    if sum(l)==1 && coefname[l][1]=="z"
        (pval, coef[l][1], @il(C[4]))
    else
        (pval, replace.(coefname[l], [r"^z"=>""]).=>coef[l], @il(C[4]))
    end
end
#}}

#{{ survfit
#[[ survfit ]]
# nothing = survfit(time, status[, group])
# Show mean and median survival time on screen. No output value.
# Equal to print(survfit(Surv(time, status) ~ group), print.rmean=TRUE) in R.
# See also: survdiff, survcoxph
# Jieyi Xiong, 2 Apr 2019

export survfit
function survfit(time::AbstractVector{Tt}, status::AbstractVector{Ts}, grp::AbstractVector) where {Tt<:Real, Ts<:Integer}
    @pkgfun(myR, rlibrary, callr)
    rlibrary("survival")
    rfun="function(time, status, group){print(survfit(Surv(time, status) ~ group), print.rmean=TRUE)}"
    callr(rfun, time, status, grp)
    nothing         
end
function survfit(time::AbstractVector{Tt}, status::AbstractVector{Ts}) where {Tt<:Real, Ts<:Integer}
    @pkgfun(myR, rlibrary, callr)
    rlibrary("survival")
    rfun="function(time, status){print(survfit(Surv(time, status) ~ 1), print.rmean=TRUE)}"
    callr(rfun, time, status)
    nothing         
end
#}}

#{{ survivalplot
#[[ survivalplot ]]
# handles = survivalplot(sfit; legend_inp=nothing|inp([labels;] ...) | legend_kw=(...))
# Or ...   = survivalplot(time, status[, group]; vb=0, ...) #Fit the model automatically.
# Draw the survival curves fitted by R package "survival".
# For R:survival package, see http://bioconnector.org/workshops/r-survival.html
# `status` could be a vector of 0/1 or 1/2 to represent censored vs “dead”.
# When legend_inp is nothing, no legend will be show.
# See also: survdiff, survcoxph
#Xiong Jieyi, 1 Jul 2017 > 3 Jan 2018 >15 Jun 2018

export survivalplot
function survivalplot(sfit; legend_inp::Union{Tuple, Void}=inp(), legend_kw=legend_inp[2])
    # myr=importpkg(:myR, preloaded=true)
    @pkgfun(myR, rlibrary, callr, callr1, relem, relemj, rnames)
    plt=importpkg(:PyPlot, preloaded=true)
    xx=relemj(sfit,"time", keepvec=true)
    yy=relemj(sfit,"surv", keepvec=true)
    ss=relemj(sfit,"n.censor", keepvec=true).>0
    num=relemj(sfit,"n", keepvec=true)
    
    is_strata=callr1("function(x,y) y %in% names(x)", sfit, "strata")
    if is_strata
        strata=[0;cumsum(relemj(sfit,"strata", keepvec=true))]
        labels=map(x->replace(x, r"^(z|group)\=" => ""), rnames(relem(sfit, "strata")))
        labels=f"$1 (N=$2)".(labels,num)
    else
        strata=[0,length(xx)]
    end
    
    hds=[]
    for i=1:length(strata)-1
        l=strata[i]+1:strata[i+1]
        x=[0;xx[l]]
        y=[1;yy[l]]
        s=[false;ss[l]]
        h=plt.step(x,y,where="post")[1] #Here draw step curves.
        plt.plot(x[s],y[s],"w+",mec=h.get_color()) #Here draw cencor dots.
        push!(hds,h)
    end
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    if is_strata && legend_inp!=nothing
        if !isempty(legend_inp[1])
            labels=legend_inp[1]
        end
        plt.legend(hds, labels; legend_kw...)
    end
    
    hds
end

function survivalplot(time::Vector{Tt}, status::Vector{Ts}, grp::Vector=[]; legend_inp::Union{Tuple, Void}=inp(), legend_kw=legend_inp[2], vb::Int=0) where {Tt<:Real, Ts<:Integer}
    ( (0 < length(time) == length(status)) && (isempty(grp) || length(time)==length(grp)) ) || error("Inputs have different lengths or zero lengths.")
    @pkgfun(myR, rlibrary, callr, callrw, eR, rshow)
    rlibrary("survival")
    if isempty(grp)
        data=callr("data.frame",time=time,status=status)
        sfit=callrw("survfit",eR("Surv(time, status)~1"),data)
    else
        data=callr("data.frame", group=grp,time=time,status=status)
        sfit=callrw("survfit",eR("Surv(time, status)~group"),data)
    end
    rshow(sfit, vb=vb)
    survivalplot(sfit; legend_inp=legend_inp, legend_kw=legend_kw)
end
#}}

#{{ liftover
#[[ liftover ]]
# chr|chrno, pos[, strand] = liftover(chr|chrno, pos, chainfile; strand=..., liftover_param="liftOver")
# Or ... = liftover((chr|chrno, pos[, strand]), chainfile; ...)
# Warped liftOver in Julia. The output always has the same row number as the input. The unliftoverable coordinates will be filled with empty values. The liftOver program and other parameters can be assigned in liftover_param.
# See also:
# Jieyi Xiong, 25 Apr 2018

export liftover
function liftover(chr::Vector{T1}, pos::AbstractMatrix{T2}, chainfile::AbstractString; strand::AbstractVector{T3}=Char[], liftover_param::AbstractString="liftOver") where T1<:AbstractString where T2<:Number where T3<:Char
    fn=tempname()
    fnout=tempname()
    open(fn,"w") do fid
        for i=1:length(chr)
            println(fid, join((chr[i],pos[i,1]-1,pos[i,2],i,255,isempty(strand) ? '+' : strand[i]), '\t'))
        end
    end
    liftover_param=isempty(liftover_param) ? "" : " "*liftover_param
    sf"$1 $2 $3 $4 /dev/null"(liftover_param, fn, chainfile, fnout)
    T=readtb(fnout, head=c"chr, N::pos, <, N::NO, X::, C::strand")
    rm(fn)
    rm(fnout)
    T["pos"][:,1]+=1
    fT=dtshift(1:length(chr), T["NO"], T, nothing)
    if isempty(strand)
        d"chr, pos"fT
    else
        d"chr, pos, strand"fT
    end
end
function liftover(chr::Vector{T1}, arg...; warg...) where T1<:Number
    C=liftover(no2chr(chr), arg...; warg...)
    tuple(chr2no(C[1]), C[2:end]...)
end
function liftover(chrpos::Tuple, arg...; warg...)
    if length(chrpos)==2
        liftover(chrpos..., arg...; warg...)
    else
        length(chrpos)==3 || error("Invalid tuple input.")
        liftover(chrpos[1], chrpos[2], arg...; strand=chrpos[3], warg...)
    end
end
#}}

#{{ findcloest
#[[ findcloest ]]
# [frontTagIdx backTagIdx]=findcloest(Base, Tag; eqin1=false)
# Find the cloest tags around each base. Base and Tag can be two numeric vectors or (chr, pos1) style. Note pos1 should be a vector rather than a N x 2 matrix. Output a 2-column index matrix with the same row number as Base. 0 in output means no any tag in thr front/back.
# When any Base and Tag are equal, eqin1 designed whether this index is reported as frontTagIdx (true) or backTagIdx (false).
# e.g. findcloest([1,3,4], [2,4]) => [0 1; 1 2; 1 2]
#      findcloest([1,3,4], [2,4]; eqin1=true) => [0 1; 1 2; 2 0]
#See also: genomemap
# Jieyi Xiong, 8 Jun 2021

export findcloest
function findcloest(A::AbstractVector{T1}, B::AbstractVector{T2}; eqin1::Bool=false) where {T1<:Real, T2<:Real}
    if eqin1
        I=sortri([B; A])
        K=[(1:length(B)); -(1:length(A))][I]
    else
        I=sortri([A; B])
        K=[-(1:length(A)); (1:length(B))][I]
    end
    O=fill(0, length(A), 2)
    p=1
    while !isnothing((ci=findnext(x->x<0, K, p);))
        fi=findprev(x->x>0, K, ci)
        bi=findnext(x->x>0, K, ci)
        O[-K[ci], :]=[isnothing(fi) ? 0 : K[fi], isnothing(bi) ? 0 : K[bi]]
        p=ci+1
    end
    return O
end
function findcloest((chrA, posA)::Tuple{Group, AbstractVector{T1}}, (chrB, posB)::Tuple{Group, AbstractVector{T2}}; eqin1::Bool=false) where {T1<:Real, T2<:Real}
    DposB=hashgrp(chrB, (posB, 1:length(posB)))
    grpfunexp(chrA, posA, id=true) do chr, cposA
        cposB, idxB=get(DposB, chr, (T2[], Int[]))
        if isempty(cposB)
            zeros(Int, length(cposA), 2)
        else
            M=findcloest(cposA, cposB; eqin1=eqin1)
            O=fill(0, size(M))
            l=M.>0
            O[l].=idxB[M[l]]
            O
        end
    end
end
#}}

# Code needs be further developed
#=
 #{{ segcloest
# #[[ segcloest ]]
#15 May 2017: Problem for this code: It's beheave is hard to be understand. The repeated segments is not be considered. It costs more coding time to use this function than do it manually by genomemap.
#
# Host_idx, Guest_idx, distance=segcloest(pos|(chr, pos); grp=group_ID, lim=(-XX, XX)|farest=XX, strand='+'s, uniqhost=false, dup=false) #A/B is the item with smaller/larger group ID in each pair.
# ...             =segcloest(Host_pos, Guest_pos; lim=(-XX, XX)|farest=XX, strand=CharVector_of_host, uniqhost=true)
#Choose the cloest pairs in any different groups. Their index are Ai and Bi. grp[Ai]!=grp[Bi]. The Ai is the group with smaller ID (host-group) and Bi is the group with larger ID (guest-group). The strand input are only for the host-group. If the guest group is in the downstream/upstream of host-group, distance is positive/negative. If host group are fully or partly overlapped with guest group, distance is 0. Only lim[1]<=distance<=lim[2] will be reported. If farest is assigned, lim=(-farest, farest). Note that if two segments are joined but not overlapped, the distance is 1/-1 rather than 0.
#In the pos-grp model, although all the strands are inputed, only the strand of host group is considered. In the pos1-pos2 model, only strand of host group is required.
#If grp is missing or be [], function will consider every row of pos are different, and the pos with smaller number are always be the host segment in each output pairs, irrelavent with the order of input pos.
#If uniqhost=true, function will try to uniquefy Ai at last and only report the pairs with the minimum abs(distance) for each Ai, unless there is more than one minimum value. In pos-grp module, this argument only work when grp is given. Note that the default values in pair mode and pos-grp mode are different.
#If dup is true, function will return ([Ai,Bi],[Bi,Ai],[D,+/-D]) rather than (Ai,Bi,D).
#See also: genomemap
#Xiong Jieyi, Sep 3, 2015

export segcloest
#<v0.6# @compat function segcloest{T<:Real}(pos::Matrix{T}; grp::Group=[], lim::Tuple{Real,Real}=(-Inf,Inf), farest::Real=NaN, strand::Vector{Char}=fill('+',size(pos,1)), uniqhost::Bool=false, dup::Bool=false)
@compat function segcloest(pos::Matrix{T}; grp::Group=[], lim::Tuple{Real,Real}=(-Inf,Inf), farest::Real=NaN, strand::Vector{Char}=fill('+',size(pos,1)), uniqhost::Bool=false, dup::Bool=false) where {T<:Real}
    if !isnan(farest)
        lim=(-farest,farest)
    end
    
    spos,idx=sortr(pos)
    sS=strand[idx]
    if isempty(grp)
        sG=1:size(pos,1)
        uniqhost=false #Not need to unique host.
    else
        sG=getrow(grp,idx)
    end

    rlt=rowpile()
    #Remove overlapped items
    len=size(spos,1)
    l=trues(len)
    op=1
    np=2
    while np<=len
        if spos[op,2]>spos[np,2]
            if getrow(sG,op)!=getrow(sG,np)
                if getrow(sG,op)>getrow(sG,np)
                    cAi=np
                    cBi=op
                else
                    cAi=op
                    cBi=np
                end
                addrow!(rlt,idx[cAi],idx[cBi],0)
            end
            l[np]=false
            np+=1
        else
            op=np
            np+=1
        end
    end
    sG=sG[l]
    sS=sS[l]
    spos=spos[l,:]
    idx=idx[l]
    
    for i=2:size(spos,1)
        Gp=getrow(sG,i-1)
        Gi=getrow(sG,i)
        if Gp==Gi
            continue
        end
        if spos[i-1,2]>=spos[i,1]
            D=0
        else
            D=spos[i,1]-spos[i-1,2]
        end
        if Gp<Gi
            cAi=i-1
            cBi=i
            rev=sS[i-1]=='-'
        else
            cAi=i
            cBi=i-1
            rev=sS[i]!='-'
        end
        if rev
            D=-D
        end
        if lim[1]<=D<=lim[2]
            addrow!(rlt,idx[cAi],idx[cBi],D)
        end
    end
    if isvirgin(rlt)
        return (Int[],Int[],T[])
    end
    Ai,Bi,D=value(rlt)
    if dup
        @assert(all((strand.=='+')|(strand.=='-')),"When dup is using, strand can only be + or -.")
        nD=D
        for i=1:length(D)
            if strand[Ai[i]]==strand[Bi[i]]
                nD[i]=-nD[i]
            end
        end
        (Ai,Bi,D)=([Ai,Bi],[Bi,Ai],[D,nD])
    end
    if uniqhost
        l=grpfunexp(x->x.==minimum(x),Ai,abs(D))
        Ai=Ai[l]
        Bi=Bi[l]
        D=D[l]
    # else
    #     if any((pos[Ai,1].==pos[Bi,1]) & (pos[Ai,2].==pos[Bi,2]))
    #         warn("Some host segment and guest segment are identical. In uniqhost=false, one of two pairs in these identical segments will lost in the report due to the algorithm limitation.")
    #     end
    end
    Ai,Bi,D
end
@compat function segcloest(chrpos::Tuple{Group,Matrix};grp=1:rownum(chrpos),strand::Vector{Char}=fill('+',rownum(chrpos)),wargs...)
    rlt,=grpfunwith(chrpos[1],chrpos[2], grp, strand, 1:rownum(chrpos), multi_output=true) do pos, cgrp, cstrand, l
        Ai,Bi,D=segcloest(pos;grp=cgrp,strand=cstrand,wargs...)
        (l[Ai],l[Bi],D)
    end
    rlt
end
@compat segcloest(chrpos::Tuple{Group,Matrix,Vector{Char}};wargs...)=segcloest(chrpos[1:2];strand=chrpos[3],wargs...)
function segcloest(pos1,pos2;grp=nothing,strand::Vector{Char}=fill('+',rownum(pos1)),uniqhost::Bool=true,dup::Bool=false,wargs...)
    @assert(grp==nothing,"If grp is given, function only accept one pos argument.")
    @assert(!dup,"dup=true is only supported in group mode so far.")
    l=[1:rownum(pos1);1:rownum(pos2)]
    grp,(pos,astrand)=vcatr_with(1:2,(pos1,strand),(pos2,fill('.',rownum(pos2))))
    Ai,Bi,D=segcloest(pos;grp=grp,strand=astrand,uniqhost=uniqhost,dup=false,wargs...)
    (l[Ai],l[Bi],D)
end
# function segcloestself(pos;grp::Void=nothing,uniqhost::Bool=false,dup::Bool=false,wargs...)
#     @assert(uniqhost==false,"In segcloestself, uniqhost can only be false.")
#     if isa(pos,Tuple) && length(pos)==3
#         strand=pos[3]
#         wargs=push!(wargs,(:strand,pos[3]))
#         pos=tuple(pos[1],pos[2])
#     end
#     Ai,Bi,D=segcloest(pos,pos;uniqhost=false,strand=strand,wargs...)
#     l=Ai.!=Bi
#     Ai=Ai[l]
#     Bi=Bi[l]
#     D=D[l]
#     if dup
#         nD=D
#         for i=1:length(D)
#             if strand[Ai[i]]==strand[Bi[i]]
#                 nD[i]=-nD[i]
#             end
#         end
#         ([Ai,Bi],[Bi,Ai],[D,nD])
#     else
#         (Ai,Bi,D)
#     end
# end
#}}
=#
