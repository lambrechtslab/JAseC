module JAseC
#{{ Preprocess
include("dataProcessKit/dataProcessKit.jl")
include("dataProcessKit/myBioXAM.jl")
using DelimitedFiles
using ArgParse
using Distributed
using .dataProcessKit
using .myBioXAM

function safesave(fn::AbstractString, arg...; kw...)
    @assert endswith(fn, ".jld2")
    jsave(fn*".part.jld2", arg...; kw...)
    mv(fn*".part.jld2", fn)
    nothing
end

global CMDARGS
global TMPDIR
global OUTDIR
global ANN
global SMP
global CELLBC
global ANN_DENOVO
# global IS10X

function parseUserInputs()
    # Parse inputted arguments
    global CMDARGS=let
        argparse_settings = ArgParseSettings()
        @add_arg_table! argparse_settings begin
            # "--mode", "-M"
            # help = "Whether handling 10x 5'-scRNAseq data or bulk RNA-seq data. Only support `10x' and `bulk' so far."
            # required = true
            # "--strand-lib"
            # help = "In bulk mode, this parameter assigns whether sequencing library is either fr-firststrand (--strand-lib=fr), fr-secondstrand (--strand-lib=rf), or non-stranded (--strand-lib=no)."
            # default = nothing
            "--ref-annotation", "-g"
            help = "GTF files of public reference annotation."
            required = true
            "--assembled-annotation"
            help = "GTF files of de-noval assembled annotation. When missing, JAseC will assemble it from bam files using StringTie."
            default=""
            "--stringtie"
            help = "Command of stringtie software."
            default=joinpath(@__DIR__, "../deps/stringtie/stringtie")
            "--temp-dir"
            help = "Tempary directionary."
            default = ""
            "--output-dir", "-o"
            help = "Directory for the outputs."
            default = "."
            "--samplesheet", "-s"
            help = "A tsv file with columns of sample name and bam files."
            required = true
            "--cell-barcode", "-c"
            help = "A tsv(.gz) file with columns of sample name, cell barcode and cell group."
            default = nothing
            "--thread", "-p"
            help = "Number of threads to run JAseC."
            arg_type = Int
            default = 1
        end
        args=parse_args(argparse_settings)
        if isempty(args["temp-dir"])
            args["temp-dir"]=joinpath(args["output-dir"], "_tmp")
        end

        ## Rerun check
        if isfile(joinpath(args["temp-dir"], "parameters.jld2"))
            lstargs=jload(joinpath(args["temp-dir"], "parameters.jld2"))
            for (k, v) in args
                if v!=lstargs[k]
                    error(f"Program detected an unfinished previous running in different parameters. To run JAseC in new parameters, please first remove the directionary $1"(args["temp-dir"]))
                end
            end
        else
            mkpath(args["temp-dir"])
            safesave(joinpath(args["temp-dir"], "parameters.jld2"), args)
        end
        args
    end

    # Add procs number
    pn=int(CMDARGS["thread"])
    println("Running JAseC in $pn cpu(s).")
    if pn>1
        @everywhere addprocs(pn) begin
	    Base.MainInclude.eval(include(@__FILE__))
	    Base.MainInclude.eval(using .JAseC)
	end
    end
    
    global TMPDIR=CMDARGS["temp-dir"]
    mkpath(TMPDIR)
    global OUTDIR=CMDARGS["output-dir"]
    mkpath(OUTDIR)
    global ANN=if endswith(CMDARGS["ref-annotation"], ".gz")
        parseGTF(sfc"zcat $1"(CMDARGS["ref-annotation"]))
    else
        parseGTF(CMDARGS["ref-annotation"])
    end
    
    #Check input and read sample sheet
    global SMP=let
        isfile(CMDARGS["ref-annotation"]) || @error(f"Cannot find ref-annotation file $1"(CMDARGS["ref-annotation"]))
        isfile(CMDARGS["samplesheet"]) || @error(f"Cannot find sample sheet file $1"(CMDARGS["samplesheet"]))
        smp=readtb(CMDARGS["samplesheet"], head=c"label, bam")
        isuniqr(smp["label"]) || @error("Sample name contains duplicates.")
        isuniqr(smp["bam"]) || @error("Bam files contain duplicates.")
        for fn in smp["bam"]
            isfile(fn) || @error("Cannot find bam file $fn")
        end
        println(f"$1 samples detected."(rnum(smp)))
        smp
    end

    #Load cell barcodes
    print("Reading cell barcodes....")
    global CELLBC=let
        isnothing(CMDARGS["cell-barcode"]) && return nothing
        isfile(CMDARGS["cell-barcode"]) || @error(f"Cannot find cell-barcode file $1"(CMDARGS["cell-barcode"]))
        fn=CMDARGS["cell-barcode"]
        T=readtb(fn, head=c"sample, cellBC, celltype")
        hashgrp(T["sample"], dd"cellBC, celltype"T)
    end
    println("Done.")
    
    #Assemble transcriptome if needed
    global ANN_DENOVO=let
        outfn=if isempty(CMDARGS["assembled-annotation"])
            joinpath(OUTDIR, "stringtie_assembled_annotation.gtf")
        else
            CMDARGS["assembled-annotation"]
        end
        if !isfile(outfn)
            println("Assembling annotation for each bam files in parallel....")
            mkpath(joinpath(TMPDIR, "stringtie_output"))
            # for (lb, bamfn) in zip(SMP["label"], SMP["bam"])
            outfnls=pmapwd(SMP["label"], SMP["bam"], env=(TMPDIR, CMDARGS)) do lb, bamfn, TMPDIR, CMDARGS
                coutfn=joinpath(TMPDIR, "stringtie_output", lb*".stringtie.gtf")
                isfile(coutfn) && return nothing
                # println("Assemble transcriptome in $lb")
                sf"$1 $2 -G $3 -o $4.part --rf"(CMDARGS["stringtie"], bamfn, CMDARGS["ref-annotation"], coutfn)
                mv(coutfn*".part", coutfn)
                coutfn
            end
            println("Done.")
            gtflistfn=joinpath(TMPDIR, "stringtie_output", "all_gtf_list.txt")
            writeall(gtflistfn, outfnls)
            print("Merging assembled transcriptome...")
            sf"$1 --merge -G $2 -o $3.part $4"(CMDARGS["stringtie"], CMDARGS["ref-annotation"], outfn, gtflistfn)
            mv(outfn*".part", outfn)
            println("Done.")
        end
        print("Reading GTF file....")
        gtf=parseGTF(outfn)
        println("Done.")
        gtf
    end
end
#}}
#{{ Main
export main
function main()
    parseUserInputs()
    
    ## Extract junctions
    println("Extracting junction reads from each bam file parallelly....")
    fetch_junction_reads()
    
    print("Unifying junction information....")
    union_junctions()
    add_unique_junction_id()
    println("Done.")

    ## Construct ASEs
    ase_info_file=joinpath(TMPDIR, "valid_splicingEvents_allAsTyp.jld2")
    T=if !isfile(ase_info_file)
        print("Detecting ASEs....")
        T=pairing_junction_into_ASE()
        T=add_gene_name(T, ANN)
        T=annotate_known_ASE(T, ANN)
        T=annotate_unknown_ASE(T, ANN, ANN_DENOVO)

        Tmxe=detect_mxe(ANN)
        if rnum(Tmxe)>0
            Tmxe=annotate_mxe(Tmxe, ANN)
            T=merge_mxe_into_ASE(T, Tmxe)
        else
            N=rnum(T)
            T["Ejc_pos_B"]=fill(0, N, 2)
            T["Ejcno_B"]=fill(0, N)
            T["Ejcno_B2"]=fill(0, N)
        end

        Tir=detect_intron_retention(ANN, ANN_DENOVO)
        if rnum(Tir)>0
            count_intron_retention_reads_in_each_cell(Tir)
            Tir=get_intron_counts(Tir)
            Tir=annotate_intron_retention(Tir, ANN)
            T=merge_intron_retention_into_ASE(T, Tir)
        end
        write_ase_info(T)
        safesave(ase_info_file, T)
        println("Done.")
        T
    else
        jload(ase_info_file)
    end
    print("Quantifying ASEs on each bam file parallelly....")
    count_reads_for_each_ASE_in_each_cell(T)
    println("Done.")

    ## Count ASE reads
    if !isnothing(CELLBC)
        print("Calculate ASEs for cell groups....")
        summarize_reads_for_each_cell_group(T)
        println("Done.")
    end

    rm(TMPDIR, recursive=true)
end
#}}
#{{ Fetch junction reads
function fetch_junction_reads()
    wkdir=joinpath(TMPDIR, "parse_junctions")
    if isdir(wkdir) && !isfile(joinpath(wkdir, ".fetch_junction_reads_unfinished"))
        return nothing
    end
    mkpath(wkdir)
    touch(joinpath(wkdir, ".fetch_junction_reads_unfinished"))
    touch(joinpath(wkdir, ".add_unique_junction_id_unfinished"))
    # for (lb, fn) in zip(d"label, bam"SMP...)
    pmapwd(d"label, bam"SMP...) do lb, fn
        outfnjc=joinpath(wkdir, lb*".junction_info.jld2")
        outfnrd=joinpath(wkdir, lb*".junction_readnum.jld2")
        if !(isfile(outfnjc) && isfile(outfnrd))
            # println("Fetch junction reads in sample $lb")
            juncT, readT=bamJuncCount_10x(fn)
            safesave(outfnjc, juncT)
            safesave(outfnrd, readT)
        end
    end
    rm(joinpath(wkdir, ".fetch_junction_reads_unfinished"))
end
function union_junctions()
    outfn=joinpath(TMPDIR, "unified_junction_info.jld2")
    if isfile(outfn)
        return
    end
    T=mapr(SMP["label"], mrows=true) do smp
        jload(joinpath(TMPDIR, "parse_junctions", "$(smp).junction_info"))
    end
    uT, _=grpfun(d"chrno, jc_pos, strand"T, T) do cT
        sT=getrow(dd"chrno, jc_pos, strand"cT, 1)
        sT["uni_readnum"]=sum(cT["uni_readnum"])
        sT
    end
    uT["uJcno"]=collect(1:rnum(uT))
    safesave(outfn, uT)
end
function add_unique_junction_id()
    wkdir=joinpath(TMPDIR, "parse_junctions")
    if isdir(wkdir) && !isfile(joinpath(wkdir, ".add_unique_junction_id_unfinished"))
        return nothing
    end
    touch(joinpath(wkdir, ".add_unique_junction_id_unfinished"))
    uT=jload(joinpath(TMPDIR, "unified_junction_info.jld2"))
    # for lb in SMP["label"]
    pmapwd(SMP["label"]) do lb
        fnjc=joinpath(wkdir, lb*".junction_info.jld2")
        fnrd=joinpath(wkdir, lb*".junction_readnum.jld2")
        Tj=jload(fnjc)
        Tj["uJcno"]=dtshift(d"chrno, jc_pos, strand"Tj, d"chrno, jc_pos, strand"uT, uT["uJcno"])
        
        Trd=jload(fnrd)
        Trd["uJcno"]=dtshift(Trd["i_jcno"], Tj["i_jcno"], Tj["uJcno"])

        jsave(fnjc, Tj)
        jsave(fnrd, Trd)
        nothing
    end
    rm(joinpath(wkdir, ".add_unique_junction_id_unfinished"))
end
#}}
#{{ Pairing junctions
function pairing_junction_into_ASE()
    jc=jload(joinpath(TMPDIR, "unified_junction_info"))
    rename!(jc, "uni_readnum"=>"num", "uJcno"=>"jcno")
    #Merge the antisense junction
    jc, _ =grpfun(d"chrno, jc_pos"jc, jc) do cjc
        if rownum(cjc)==1
            cjc["jcno2"]=0
            cjc
        else
            i=argmax(cjc["num"])
            T=getrow(cjc, i)
            T["num"]=sum(cjc["num"])
            T["jcno2"]=cjc["jcno"][3-i]
            T
        end
    end

    jcL=tbx(jc, jc, ((jc["chrno"], jc["jc_pos"][:, 1]), (jc["chrno"], jc["jc_pos"][:, 1])), aprefix="E", bprefix="I")
    jcL=rec(jcL, jcL["Ijc_pos"][:, 2].<jcL["Ejc_pos"][:, 2])

    jcR=tbx(jc, jc, ((jc["chrno"], jc["jc_pos"][:, 2]), (jc["chrno"], jc["jc_pos"][:, 2])), aprefix="E", bprefix="I")
    jcR=rec(jcR, jcR["Ijc_pos"][:, 1].>jcR["Ejc_pos"][:, 1])

    jcT=tbx(jcL, jcR, key=c"chrno, Ejcno", aprefix="L", bprefix="R", fillempty=:AB)
    jcT=rec(jcT, (jcT["LIjcno"].==0) .| (jcT["RIjcno"].==0) .| (jcT["LIjc_pos"][:, 2].<jcT["RIjc_pos"][:, 1]))

    #Assign strand for ASE
    jcT["strand"]=rowfun(jcT) do T
        S=@with T [$Estrand, $LIstrand, $RIstrand]
        N=@with T [$Enum, $LInum, $RInum]
        uN, uS=grpfun(sum, S, N)
        l=findall(uN.==maximum(uN))
        if length(l)==1
            uS[l[1]]
        else
            T["Estrand"]
        end
    end
    fd_i!(jcT, c"Estrand, LIstrand, RIstrand")
    
    jcT["Inum"]=max.(jcT["LInum"], jcT["RInum"])

    fstGrpV=grpfunexp(linkgrp, jcT["chrno"], [jcT["Ejc_pos"] jcT["LIjc_pos"][:, 2] jcT["RIjc_pos"][:, 1]])
    l=grpfunexp(
        fstGrpV,
        hcat(min.(jcT["Enum"], jcT["Inum"]),
             max.(jcT["Enum"], jcT["Inum"]),
         -abs.(jcT["LInum"].-jcT["RInum"])),
        [jcT["Ejc_pos"] jcT["LIjc_pos"][:, 2] jcT["RIjc_pos"][:, 1]]) do mrdnum, Mpos
            len=size(Mpos, 1)
            I=sortri(mrdnum)
            ckl=trues(len) #for I order
            kpL=falses(len) #for original order
            cgrpv=ones(Int, len)
            pi=len
            px=I[pi]
            while mrdnum[px, 1]>=5
                kpL[px]=true
                t=(Mpos[I, 1].==Mpos[px, 1]) .| (Mpos[I, 2].==Mpos[px, 2])
                if Mpos[px, 3]>0
                    t=t .| (Mpos[I, 2].==Mpos[px, 3])
                end
                if Mpos[px, 4]>0
                    t=t .| (Mpos[I, 1].==Mpos[px, 4])
                end
                ckl[t].=false
                pi=findlast(ckl)
                isnothing(pi) && break
                px=I[pi]
            end
            kpL
        end
    jcT=rec(jcT, l)

    jcT=rec(jcT, (jcT["Enum"].>=5) .& (jcT["Inum"].>=5) .& (0.01.<(jcT["Inum"]./(jcT["Enum"].+jcT["Inum"])).<0.99))
    
    # safesave("valid_splicingEvents.jld2", jcT)
    return jcT
end
#}}
#{{ Add gene name
function add_gene_name(jc, exn)
    # jc=jload("valid_splicingEvents.jld2")
    # exn=jload(exnfn)

    #Remove non-spliced transcripts
    l=grpfunexp(isrowsame, exn["transid"], exn["exn_pos"])
    exn=rec(exn, .!l)
    
    #Add geneid
    texn=rec(exn, .!isempty.(exn["geneid"]))
    Bl=hashgrp((texn["chrno"], texn["exn_pos"][:, 1]), texn["geneid"])
    Br=hashgrp((texn["chrno"], texn["exn_pos"][:, 2]), texn["geneid"])
    jc["geneid"]=rowfun(jc) do T
        Gel=get(Br, (T["chrno"], T["Ejc_pos"][1]), String[])
        Ger=get(Bl, (T["chrno"], T["Ejc_pos"][2]), String[])
        # Gil=get(Bl, (T["chrno"], T["LIjc_pos"][1]), String[]) #bug?? Should be T["LIjc_pos"][2]
        Gil=get(Bl, (T["chrno"], T["LIjc_pos"][2]), String[])
        # Gir=get(Br, (T["chrno"], T["RIjc_pos"][2]), String[]) #bug?? Should be T["RIjc_pos"][1]
        Gir=get(Br, (T["chrno"], T["RIjc_pos"][1]), String[])
        setno, gels=vcatr_with(1:4, Gel, Ger, Gil, Gir)
        if isempty(gels)
            ""
        else
            setfq=idennum(setno)
            gefq=idennum(gels)
            gels[sortri((setfq, gefq), rev=true)[1]]
        end
    end

    #Using gene strand will correct 0.81% ASE.
    d"gene_name, gene_typ, strand"jc=dtshift(jc["geneid"], exn["geneid"], d"gene_name, gene_typ, strand"exn, (fill("", rnum(jc)), fill("", rnum(jc)), jc["strand"]))
    
    # safesave(outfn, jc)
    return jc
end
#}}
#{{ Annotate known AS events
function annotate_known_ASE(jc, exn)
    # outfn="valid_splicingEvents_withAsType_step1.jld2"
    # isfile(outfn) && return nothing
    
    #Preprocess annotation
    # exn=jload(exnfn)
    exn["exn_rvrank"]=grpfunexp(exn["transid"], exn["exn_rank"]) do x
        (length(x).-x).+1
    end
    exn["isend"]=(exn["exn_rank"].==1) .| (exn["exn_rvrank"].==1)
    exn, _=grpfunwith(exn["transid"], exn) do ex
        ex=rec(ex, sortri(ex["exn_rank"]))
        P=0
        flag=false
        ex["loc_typ"], ex["cds_loc"]=if all(ex["cds_pos"][:, 1].==0)
            reprow(("NC", [0 0]), rownum(ex))
        else
            rowfun(ex["cds_pos"], as_vec=true) do x
                if x[1]==0
                    (flag ? "UTR3" : "UTR5", [0 0])
                else
                    cdslen=x[2]-x[1]+1
                    out=("CDS", [P+1 P+cdslen])
                    P=P+cdslen
                    flag=true
                    out
                end
            end
        end
        ex
    end
    
    trEndL=Dict{Tuple{String, Int}, Int}()
    grploop((exn["transid"], exn["exn_pos"][:, 1]), exn["exn_pos"][:, 2], id=true) do ky, x
        trEndL[ky]=maximum(x)
    end
    
    trEndR=Dict{Tuple{String, Int}, Int}()
    grploop((exn["transid"], exn["exn_pos"][:, 2]), exn["exn_pos"][:, 1], id=true) do ky, x
        trEndR[ky]=minimum(x)
    end

    # D_transToExnPos=hashgrp(exn["transid"], exn["exn_pos"]) #@1
    
    exn=rec(exn, sortri(d"transid, exn_pos"exn))
    d"cdsLenToLeft, cdsLenToRight"exn=grpfunexp(exn["transid"], seglen(exn["cds_pos"])) do x
        (cumsum(x), reverse(cumsum(reverse(x))))
    end
    exn["hasCdsDwstream"]=grpfunexp(exn["transid"], exn["cds_pos"], exn["strand"]) do cP, strand
        k=findfirst(cP[:, 1].>0)
        if isnothing(k)
            false
        elseif strand[1]=='+'
            falsesbut(length(strand), 1:k)
        else
            @assert strand[1]=='-'
            falsesbut(length(strand), k:length(strand))
        end
    end
    DL=hashgrp((exn["chrno"], exn["exn_pos"][:, 2]))
    DR=hashgrp((exn["chrno"], exn["exn_pos"][:, 1]))

    empexn=zerorow(exn)
    
    pool=Tuple{Dict, Float64}[]
    
    # jc=jload("valid_splicingEvents_withGene")
    jc=rowfun(jc, no=true) do loopi, cjc
        empty!(pool)
        
        if mod(loopi, 100)==0
            if mod(loopi, 1000)==0
                # println(loopi)
            else
                # print(".")
            end
        end
        getexn=(D, chrpos)->begin
            l=get(D, chrpos, nothing)
            if l==nothing
                empexn
            else
                getrow(exn, l)
            end
        end
        exnEL=getexn(DL, (cjc["chrno"], cjc["Ejc_pos"][1, 1]))
        exnER=getexn(DR, (cjc["chrno"], cjc["Ejc_pos"][1, 2]))
        exnIL=getexn(DR, (cjc["chrno"], cjc["LIjc_pos"][1, 2]))
        exnIR=getexn(DL, (cjc["chrno"], cjc["RIjc_pos"][1, 1]))
        #CE
        let trS=intersect(exnEL["transid"], exnER["transid"], exnIL["transid"], exnIR["transid"])
            if !isempty(trS)
                cjc=deepcopy(cjc)
                cjc["AS_typ"]="CE"
                t=vcat(exnEL["geneid"], exnER["geneid"], exnIL["geneid"], exnIR["geneid"])
                cdslen=if any(exnIL["loc_typ"].=="CDS") || any(exnIR["loc_typ"].=="CDS")
                    exnIL, exnIR=tbx(exnIL, exnIR, (exnIL["transid"], exnIR["transid"]), nomerge=true)
                    #() cdsCgLen=seglen([exnIL["cds_loc"][:, 1] exnIR["cds_loc"][:, 2]]) #A bug: cdsCgLen could be negative. 24 May 2020
                    t=[exnIL["cds_loc"] exnIR["cds_loc"]]
                    M=t[anyh(t.>0), :]
                    if size(M, 1)>0
                        cdsCgLen=rowfun(M, as_vec=true) do x
                            t=extrema(x[x.>0])
                            t[2]-t[1]+1
                        end
                        if any((l=mod.(cdsCgLen, 3).==0;))
                            maximum(cdsCgLen[l])
                        else
                            maximum(cdsCgLen)
                        end                    
                    else
                        0
                    end
                else
                    0
                end
                cjc["cds_len_EI"]=[0 cdslen]

                trES=intersect(exnEL["transid"], exnER["transid"])
                cjc["isKnownEI"]=if isempty(setdiff(trES, trS))
                    # "reComb"
                    [false true]
                else
                    # "known"
                    [true true]
                end
                cjc["cds_cg"]=if cdslen==0
                    ""
                elseif any(exnEL["cds_pos"][:, 1].>0) && any(exnER["cds_pos"][:, 1].>0)
                    if mod(cdslen, 3)==0
                        "in-frame"
                    else
                        "frame-shift"
                    end
                elseif (cjc["strand"]=='+' && any(exnEL["cds_pos"][:, 1].>0)) || (cjc["strand"]=='-' && any(exnER["cds_pos"][:, 1].>0))
                    "tail-cut"
                elseif (cjc["strand"]=='+' && any(exnER["hasCdsDwstream"])) || (cjc["strand"]=='-' && any(exnEL["hasCdsDwstream"]))
                    "head-cut"
                else
                    "nonsense"
                end
                push!(pool, (cjc, 10*min(cjc["LInum"], cjc["RInum"])/max(cjc["LInum"], cjc["RInum"])))
            end
        end

        trE=intersect(exnEL["transid"], exnER["transid"])
        isKnownE=!isempty(trE)
        #L-pairs
        let trI=setdiff(intersect(exnEL["transid"], exnIL["transid"]), trE)
            # if !isempty(trI) #@1
            #     l=map(trI) do trid
            #         exnpos=D_transToExnPos[trid]
            #         any((cjc["Ejc_pos"][1] .< exnpos[:, 1]) .& (exnpos[:, 2] .< cjc["LIjc_pos"][2]))
            #     end
            #     deleteat!(trI, l)
            # end
            if !isempty(trI) && cjc["LInum"]>cjc["RInum"]
                cjc=deepcopy(cjc)
                # cjc["knowledge"]="known"
                cjc["isKnownEI"]=[isKnownE true]
                eIL=getrow(exnIL, ismbr(exnIL["transid"], trI))
                if any(eIL["exn_pos"][:, 2].>cjc["Ejc_pos"][2]) #ASS
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "ASS5" : "ASS3"
                    exnIL=rec(exnIL, ismbr(exnIL["transid"], trI) .& (exnIL["loc_typ"].=="CDS"))
                    cdslen=if rownum(exnIL)>0
                        if any((l=0 .< exnIL["cds_pos"][:, 1] .< cjc["Ejc_pos"][2];))
                            t=cjc["Ejc_pos"][2].-exnIL["cds_pos"][l, 1]
                            if any((l=mod.(t, 3)==0;))
                                maximum(t[l])
                            else
                                maximum(t)
                            end
                        else
                            0
                        end
                    else
                        0
                    end
                    cjc["cds_len_EI"]=[0 cdslen]
                    cjc["cds_cg"]=if cdslen==0
                        ""
                    elseif any(exnEL["cds_pos"][:, 1].>0) && any(exnEL["cds_pos"][:, 1].>0)
                        if mod(cdslen, 3)==0
                            "in-frame"
                        else
                            "frame-shift"
                        end
                    elseif (cjc["strand"]=='+' && any(exnEL["cds_pos"][:, 1].>0)) || (cjc["strand"]=='-' && any(exnER["cds_pos"][:, 1].>0))
                        "tail-cut"
                    elseif (cjc["strand"]=='+' && any(exnER["hasCdsDwstream"])) || (cjc["strand"]=='-' && any(exnEL["hasCdsDwstream"]))
                        "head-cut"
                    else
                        "nonsense"
                    end
                elseif maximum(map((x, y)->trEndL[(x, y)], eIL["transid"], eIL["exn_pos"][:, 1]))<cjc["Ejc_pos"][2]
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "altF" : "altL"
                    # t=eIL["cdsLenToRight"]
                    # if any((l=(t.>0) .& (mod.(t, 3).==0);))
                    #     t=t[l]
                    # end
                    # Icdslen=maximum(t)
                    # t=exnER["cdsLenToRight"]
                    # if any((l=(t.>0) .& (mod.(t, 3).==0);))
                    #     t=t[l]
                    # end
                    # Ecdslen=isempty(t) ? -1 : maximum(t) # exnER could be empty.
                    tI=eIL["cdsLenToRight"]
                    mI=map(x->x>0 ? mod(x, 3) : -1, tI)
                    tE=exnER["cdsLenToRight"]
                    mE=map(x->x>0 ? mod(x, 3) : -1, tE)
                    mU=intersect(mI[mI.>=0], mE[mE.>=0])
                    if !isempty(mU)
                        tI=tI[ismbr(mI, mU)]
                        tE=tE[ismbr(mE, mU)]
                    end
                    Icdslen=maximum(tI)
                    Ecdslen=isempty(tE) ? 0 : maximum(tE) # exnER could be empty.
                    # Ecdslen=maximum(tE)
                    
                    cjc["cds_len_EI"]=[Ecdslen Icdslen]
                    cjc["cds_cg"]=if Ecdslen==0 && Icdslen==0
                        ""
                    elseif !any(exnEL["cds_pos"][:, 1].>0)
                        "nonsense"
                    elseif cjc["AS_typ"]=="altF"
                        if Ecdslen>0 && Icdslen>0
                            if mod(abs(Ecdslen-Icdslen), 3)==0
                                "in-frame"
                            else
                                "frame-shift"
                            end
                        elseif any(exnEL["hasCdsDwstream"])
                            "head-cut"
                        else
                            "nonsense"
                        end                        
                    else
                        if Ecdslen>0 && Icdslen>0
                            if mod(abs(Ecdslen-Icdslen), 3)==0
                                "in-frame"
                            else
                                "frame-shift"
                            end
                        else
                            "tail-cut"
                        end                        
                    end
                    # return cjc
                else
                    cjc["AS_typ"]="Mix"
                    cjc["cds_len_EI"]=[-1 -1]
                    cjc["cds_cg"]=""
                    # return cjc
                end
                push!(pool, (cjc, 1-cjc["RInum"]/cjc["LInum"]))
            end
        end
        
        #R-pairs
        let trI=setdiff(intersect(exnER["transid"], exnIR["transid"]), trE)
            # if !isempty(trI) #@1
            #     l=map(trI) do trid
            #         exnpos=D_transToExnPos[trid]
            #         any((cjc["RIjc_pos"][1] .< exnpos[:, 1]) .& (exnpos[:, 2] .< cjc["Ejc_pos"][2]))
            #     end
            #     deleteat!(trI, l)
            # end
            if !isempty(trI) && cjc["LInum"]<cjc["RInum"]
                cjc=deepcopy(cjc)
                cjc["isKnownEI"]=[isKnownE true]
                eIR=getrow(exnIR, ismbr(exnIR["transid"], trI))
                if any(eIR["exn_pos"][:, 1].<cjc["Ejc_pos"][1])
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "ASS3" : "ASS5"
                    exnIR=rec(exnIR, ismbr(exnIR["transid"], trI) .& (exnIR["loc_typ"].=="CDS"))
                    cdslen=if rownum(exnIR)>0
                        if any((l=eIR["cds_pos"][:, 2].>cjc["Ejc_pos"][1];))
                            t=eIR["cds_pos"][l, 2].-cjc["Ejc_pos"][1]
                            if any((l=mod.(t, 3)==0;))
                                maximum(t[l])
                            else
                                maximum(t)
                            end
                        else
                            0
                        end
                    else
                        0
                    end
                    cjc["cds_len_EI"]=[0 cdslen]
                    cjc["cds_cg"]=if cdslen==0
                        ""
                    elseif any(exnEL["cds_pos"][:, 1].>0) && any(exnER["cds_pos"][:, 1].>0)
                        if mod(cdslen, 3)==0
                            "in-frame"
                        else
                            "frame-shift"
                        end
                    elseif (cjc["strand"]=='+' && any(exnEL["cds_pos"][:, 1].>0)) || (cjc["strand"]=='-' && any(exnER["cds_pos"][:, 1].>0))
                        "tail-cut"
                    elseif (cjc["strand"]=='+' && any(exnER["hasCdsDwstream"])) || (cjc["strand"]=='-' && any(exnEL["hasCdsDwstream"]))
                        "head-cut"
                    else
                        "nonsense"
                    end
                    # return cjc
                elseif minimum(map((x, y)->trEndR[(x, y)], eIR["transid"], eIR["exn_pos"][:, 2]))>cjc["Ejc_pos"][1]
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "altL" : "altF"
                    # t=eIR["cdsLenToLeft"]
                    # if any((l=(t.>0) .& (mod.(t, 3).==0);))
                    #     t=t[l]
                    # end
                    # Icdslen=maximum(t)
                    # t=exnEL["cdsLenToLeft"]
                    # if any((l=(t.>0) .& (mod.(t, 3).==0);))
                    #     t=t[l]
                    # end
                    # Ecdslen=isempty(t) ? -1 : maximum(t)
                    
                    tI=eIR["cdsLenToLeft"]
                    mI=map(x->x>0 ? mod(x, 3) : -1, tI)
                    tE=exnEL["cdsLenToLeft"]
                    mE=map(x->x>0 ? mod(x, 3) : -1, tE)
                    mU=intersect(mI[mI.>=0], mE[mE.>=0])
                    if !isempty(mU)
                        tI=tI[ismbr(mI, mU)]
                        tE=tE[ismbr(mE, mU)]
                    end
                    Icdslen=maximum(tI)
                    Ecdslen=isempty(tE) ? 0 : maximum(tE) # exnER could be empty.
                    # Ecdslen=maximum(tE)
                    
                    cjc["cds_len_EI"]=[Ecdslen Icdslen]
                    cjc["cds_cg"]=if Ecdslen==0 && Icdslen==0
                        ""
                    elseif !any(exnER["cds_pos"][:, 1].>0)
                        "nonsense"
                    elseif cjc["AS_typ"]=="altF"
                        if Ecdslen>0 && Icdslen>0
                            if mod(abs(Ecdslen-Icdslen), 3)==0
                                "in-frame"
                            else
                                "frame-shift"
                            end
                        elseif any(exnER["hasCdsDwstream"])
                            "head-cut"
                        else
                            "nonsense"
                        end                        
                    else
                        if Ecdslen>0 && Icdslen>0
                            if mod(abs(Ecdslen-Icdslen), 3)==0
                                "in-frame"
                            else
                                "frame-shift"
                            end
                        else
                            "tail-cut"
                        end                        
                    end
                else
                    cjc["AS_typ"]="Mix"
                    cjc["cds_len_EI"]=[-1 -1]
                    cjc["cds_cg"]=""
                    # return cjc
                end
                push!(pool, (cjc, 1-cjc["LInum"]/cjc["RInum"]))
            end
        end
        
        if isempty(pool)
            cjc["AS_typ"]="UN"
            cjc["isKnownEI"]=[false false]
            cjc["cds_len_EI"]=[-1 -1]
            cjc["cds_cg"]=""
            return cjc
        else
            t1, t2=disent(pool)
            return t1[argmax(t2)]
        end
    end
    
    # jsave(outfn, jc)
    return jc
end
#}}
#{{ Annotate unknonwn AS events
function annotate_unknown_ASE(jc0, exn_known, exn)
    #Preprocess annotation
    # exn=jload(annfn)
    exn["exn_rvrank"]=grpfunexp(exn["transid"], exn["exn_rank"]) do x
        (length(x).-x).+1
    end
    exn["isend"]=(exn["exn_rank"].==1) .| (exn["exn_rvrank"].==1)
    
    DL=hashgrp((exn["chrno"], exn["exn_pos"][:, 2], exn["strand"]))
    DR=hashgrp((exn["chrno"], exn["exn_pos"][:, 1], exn["strand"]))
    
    trEndL=Dict{Tuple{String, Int}, Int}()
    grploop((exn["transid"], exn["exn_pos"][:, 1]), exn["exn_pos"][:, 2], id=true) do ky, x
        trEndL[ky]=maximum(x)
    end
    
    trEndR=Dict{Tuple{String, Int}, Int}()
    grploop((exn["transid"], exn["exn_pos"][:, 2]), exn["exn_pos"][:, 1], id=true) do ky, x
        trEndR[ky]=minimum(x)
    end
    
    empexn=zerorow(exn)
    pool=Tuple{Dict, Float64}[]
    
    # jc0=jload("valid_splicingEvents_withAsType_step1")
    jc=rec(jc0, (jc0["AS_typ"].!="CE") .& .!allh(jc0["isKnownEI"]))

    
    jc=rowfun(jc, no=true) do loopi, cjc
        empty!(pool)
        if mod(loopi, 100)==0
            if mod(loopi, 1000)==0
                println(loopi)
            else
                print(".")
            end
        end
        getexn=(D, chrpos)->begin
            l=get(D, chrpos, nothing)
            if l==nothing
                empexn
            else
                getrow(exn, l)
            end
        end
        exnEL=getexn(DL, (cjc["chrno"], cjc["Ejc_pos"][1, 1], cjc["strand"]))
        exnER=getexn(DR, (cjc["chrno"], cjc["Ejc_pos"][1, 2], cjc["strand"]))
        exnIL=getexn(DR, (cjc["chrno"], cjc["LIjc_pos"][1, 2], cjc["strand"]))
        exnIR=getexn(DL, (cjc["chrno"], cjc["RIjc_pos"][1, 1], cjc["strand"]))

        #CE
        let trS=intersect(exnEL["transid"], exnER["transid"], exnIL["transid"], exnIR["transid"])
            if !isempty(trS)
                cjc=deepcopy(cjc)
                cjc["AS_typ"]="CE"
                push!(pool, (cjc, 10*min(cjc["LInum"], cjc["RInum"])/max(cjc["LInum"], cjc["RInum"])))
            end
        end

        trE=intersect(exnEL["transid"], exnER["transid"])
        #L-pairs
        let trI=setdiff(intersect(exnEL["transid"], exnIL["transid"]), trE)
            if !isempty(trI)
                cjc=deepcopy(cjc)
                eIL=getrow(exnIL, ismbr(exnIL["transid"], trI))
                if any(eIL["exn_pos"][:, 2].>cjc["Ejc_pos"][2]) #ASS
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "ASS5" : "ASS3"
                    # return cjc
                elseif maximum(map((x, y)->trEndL[(x, y)], eIL["transid"], eIL["exn_pos"][:, 1]))<cjc["Ejc_pos"][2]
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "altF" : "altL"
                    # return cjc
                else
                    cjc["AS_typ"]="Mix"
                    # return cjc
                end
                push!(pool, (cjc, 1-cjc["RInum"]/cjc["LInum"]))
            end
        end
        
        #R-pairs
        let trI=setdiff(intersect(exnER["transid"], exnIR["transid"]), trE)
            if !isempty(trI)
                cjc=deepcopy(cjc)
                eIR=getrow(exnIR, ismbr(exnIR["transid"], trI))
                if any(eIR["exn_pos"][:, 1].<cjc["Ejc_pos"][1])
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "ASS3" : "ASS5"
                    # return cjc
                elseif minimum(map((x, y)->trEndR[(x, y)], eIR["transid"], eIR["exn_pos"][:, 2]))>cjc["Ejc_pos"][1]
                    cjc["AS_typ"]=cjc["strand"]=='-' ? "altL" : "altF"
                    # return cjc
                else
                    cjc["AS_typ"]="Mix"
                    # return cjc
                end
                push!(pool, (cjc, 1-cjc["LInum"]/cjc["RInum"]))
            end
        end

        if isempty(pool)
            cjc["AS_typ"]="UN"
            return cjc
        else
            t1, t2=disent(pool)
            return t1[argmax(t2)]
        end
    end
    jc=rec(jc, jc["AS_typ"].!="UN")
    t=copy(jc0["AS_typ"])
    jc0["AS_typ"]=dtshift!(jc0["Ejcno"], jc["Ejcno"], jc["AS_typ"], jc0["AS_typ"], safe=true)
    jc0["cds_cg"][t.!=jc0["AS_typ"]].="" #If AS type changed, remove previous cds_cg data. 28 Jul 2021
    jc0["cds_len_EI"][t.!=jc0["AS_typ"], :].=0
    
    jc=rec(jc0, ismbr(jc0["AS_typ"], c"CE, ASS3, ASS5") .& allh((jc0["cds_len_EI"].==0)))
    
    #Check if noval exon are CDS
    # exn=jload(exnfn)
    cds=rec(exn_known, exn_known["cds_pos"][:, 1].>0)

    l1=ismbr((jc["chrno"], jc["Ejc_pos"][:, 1], jc["strand"]), (cds["chrno"], cds["exn_pos"][:, 2], cds["strand"]))
    l2=ismbr((jc["chrno"], jc["Ejc_pos"][:, 2], jc["strand"]), (cds["chrno"], cds["exn_pos"][:, 1], cds["strand"]))
    l=l1 .& l2
    jc["cds_len_EI"][l, :]=rowfun(jc["Ejc_pos"][l, :], jc["LIjc_pos"][l, :], jc["RIjc_pos"][l, :], as_vec=true) do e, l, r
        P1=l[2]==0 ? e[1] : max(e[1], l[2])
        P2=r[1]==0 ? e[2] : min(e[2], r[1])
        [0 P2-P1+1]
    end
    jc0["cds_len_EI"]=dtshift!(jc0["Ejcno"], jc["Ejcno"], jc["cds_len_EI"], jc0["cds_len_EI"], safe=true)
    
    #For ASS, altF/L, remove the junction in one side
    (ff!)=(T, f)->begin
        T["$(f)Ijc_pos"]=[0 0]
        T["$(f)Ijcno"]=0
        T["$(f)Ijcno2"]=0
        T["$(f)Inum"]=0
        T
    end
    jc0=rowfun(jc0) do T
        if T["strand"]=='+'
            if T["AS_typ"]=="ASS3" || T["AS_typ"]=="altL"
                T["LInum"]>T["RInum"] || print("A")
                ff!(T, "R")
            elseif T["AS_typ"]=="ASS5" || T["AS_typ"]=="altF"
                T["LInum"]<T["RInum"] || print("B")
                ff!(T, "L")
            end
        else
            @assert T["strand"]=='-'
            if T["AS_typ"]=="ASS3" || T["AS_typ"]=="altL"
                T["LInum"]<T["RInum"] || print("C")
                ff!(T, "L")
            elseif T["AS_typ"]=="ASS5" || T["AS_typ"]=="altF"
                T["LInum"]>T["RInum"] || print("D")
                ff!(T, "R")
            end
        end
        T
    end
    
    # jsave("valid_splicingEvents_withAsType", jc0)
    return jc0
end
#}}
#{{ Define MXE, method
  #{{ Detected MXE
function detect_mxe(exn)
    # exn=jload(annfn)
    # exnitn=rowpile()
    iexn, _=grpfunwith(exn["transid"], exn) do T
        if rownum(T)<3
            return nothing
        end
        T=rec(T, sortri(T["exn_pos"]))
        T["itnrg_pos"]=[[0; T["exn_pos"][1:end-1, 2]] [T["exn_pos"][2:end, 1]; 0]]
        getrow(T, 2:rownum(T)-1)
    end
    l=grpfunexp(d"geneid, itnrg_pos"iexn, iexn["exn_pos"]) do x
        size(segunion(x)[1], 1)>=2
    end
    iexn=rec(iexn, l)
    iexn=tbuniq(dd"chrno, geneid, strand, exn_pos, itnrg_pos"iexn)

    jc=jload(joinpath(TMPDIR, "unified_junction_info"))
    ujc=tb()
    d"jcno, num"ujc, d"chrno, jc_pos"ujc=grpfun(d"chrno, jc_pos"jc, jc["uJcno"], jc["uni_readnum"]) do jcno, rdnum
        if length(jcno)==1
            ([jcno 0], rdnum)
        else
            t=argmax(rdnum)
            ([jcno[t] jcno[3-t]], sum(rdnum))
        end
    end

    t1, iexn["Ljcno"]=dtshift((iexn["chrno"], [iexn["itnrg_pos"][:, 1] iexn["exn_pos"][:, 1]]),
                              d"chrno, jc_pos"ujc, d"num, jcno"ujc, (0, [0 0]))
    t2, iexn["Rjcno"]=dtshift((iexn["chrno"], [iexn["exn_pos"][:, 2] iexn["itnrg_pos"][:, 2]]),
                              d"chrno, jc_pos"ujc, d"num, jcno"ujc, (0, [0 0]))
    iexn["LRjc_num"]=[t1 t2]
    mxe, _=grpfunwith(d"geneid, itnrg_pos"iexn, iexn) do T
        ri=sortri([minh(T["LRjc_num"]) maxh(T["LRjc_num"])])[1:2]
        ai, bi=T["exn_pos"][ri[1], 1]<T["exn_pos"][ri[2], 1] ? ri : (ri[2], ri[1])
        tb(chrno=T["chrno"][ai],
           strand=T["strand"][ai],
           Ejcno=T["Ljcno"][ai, 1],
           Ejcno2=T["Ljcno"][ai, 2],
           Enum=T["LRjc_num"][ai, 1],
           Enum_B=T["LRjc_num"][ai, 2],
           Ejcno_B=T["Rjcno"][ai, 1],
           Ejcno_B2=T["Rjcno"][ai, 2],
           LIjcno=T["Ljcno"][bi, 1],
           LIjcno2=T["Ljcno"][bi, 2],
           RIjcno=T["Rjcno"][bi, 1],
           RIjcno2=T["Rjcno"][bi, 2],
           LInum=T["LRjc_num"][bi, 1],
           RInum=T["LRjc_num"][bi, 2],
           Ejc_pos=[T["itnrg_pos"][ai, 1] T["exn_pos"][ai, 1]],
           Ejc_pos_B=[T["exn_pos"][ai, 2] T["itnrg_pos"][ai, 2]],
           LIjc_pos=[T["itnrg_pos"][bi, 1] T["exn_pos"][bi, 1]],
           RIjc_pos=[T["exn_pos"][bi, 2] T["itnrg_pos"][bi, 2]],
           AS_typ="MXE")
    end

    #Remove fake MXE
    l=mxe["Ejc_pos_B"][:, 1] .< mxe["LIjc_pos"][:, 2]
    mxe=rec(mxe, l)

    #Remove the MXE which both two exons can be included in one transcript.
    # intron=tb(c"transid, itn_pos"=>findintron(exn["transid"], exn["exn_pos"]))
    # d"chrno, strand"intron=dtshift(intron["transid"], exn["transid"], d"chrno, strand"exn)
    # l=ismbr((mxe["chrno"], [mxe["Ejc_pos_B"][:, 1] mxe["LIjc_pos"][:, 2]]),
    #         (intron["chrno"], intron["itn_pos"]|>x->[x[:, 1].-1 x[:, 2].+1]))
    # mxe=rec(mxe, .!l)
    # 30 Apr 2020 Test below: >>>
    mxe["nonMxJcNum"]=dtshift((mxe["chrno"], [mxe["Ejc_pos_B"][:, 1] mxe["LIjc_pos"][:, 2]]), d"chrno, jc_pos"ujc, ujc["num"], 0)
    l=min.(mxe["Enum_B"], mxe["LInum"]).>=mxe["nonMxJcNum"]
    mxe=rec(mxe, l)
    #<<<
    
    mxe=rec(mxe, (mxe["Enum"].+mxe["Enum_B"].>=5) .& (mxe["LInum"].+mxe["RInum"].>=5) .&
            (0.01 .< (mxe["LInum"].+mxe["RInum"])./(mxe["Enum"].+mxe["Enum_B"].+mxe["LInum"].+mxe["RInum"]) .< 0.99) .&
            (max.(mxe["Enum"], mxe["Enum_B"]) .< 20*min.(mxe["Enum"], mxe["Enum_B"])) .&
            (max.(mxe["LInum"], mxe["RInum"]) .< 20*min.(mxe["LInum"], mxe["RInum"])))
    
    # jsave("valid_MXE_splicingEvents", mxe)
    return mxe
end
#}}
  #{{ Annotate MXE
function annotate_mxe(mxe, exn)
    # exn=jload(exnfn)
    iexn, _=grpfunwith(exn["transid"], exn) do T
        if rownum(T)<3
            return nothing
        end
        T=rec(T, sortri(T["exn_pos"]))
        T["itnrg_pos"]=[[0; T["exn_pos"][1:end-1, 2]] [T["exn_pos"][2:end, 1]; 0]]
        getrow(T, 2:rownum(T)-1)
    end
    l=grpfunexp(x->falsesbut(length(x), argmax(x)),
                d"geneid, itnrg_pos, exn_pos"iexn, seglen(iexn["cds_pos"]))
    iexn=rec(iexn, l)

    # mxe=jload("valid_MXE_splicingEvents")
    l1=memberr((mxe["chrno"], hcat(d"Ejc_pos, Ejc_pos_B"mxe...)),
               (iexn["chrno"], [iexn["itnrg_pos"][:, 1] iexn["exn_pos"] iexn["itnrg_pos"][:, 2]]))
    l2=memberr((mxe["chrno"], hcat(d"LIjc_pos, RIjc_pos"mxe...)),
               (iexn["chrno"], [iexn["itnrg_pos"][:, 1] iexn["exn_pos"] iexn["itnrg_pos"][:, 2]]))
    # mxe["knowledge"]=ifelse.((l1.>0) .& (l2.>0), "known", "novel")
    mxe["isKnownEI"]=[l1.>0 l2.>0]
    l3=memberr((mxe["chrno"], [mxe["Ejc_pos"][:, 1] mxe["Ejc_pos_B"][:, 2]]), d"chrno, itnrg_pos"iexn)
    mxe["geneid"]=fill("", rnum(mxe))
    mxe["geneid"][l3.>0].=iexn["geneid"][l3[l3.>0]]
    mxe["geneid"][l2.>0].=iexn["geneid"][l2[l2.>0]]
    mxe["geneid"][l1.>0].=iexn["geneid"][l1[l1.>0]]

    ll=(l1.>0) .& (l2.>0)
    mxe["cds_len_EI"]=fill(0, rnum(mxe), 2)
    mxe["cds_len_EI"][ll, 1]=seglen(iexn["cds_pos"][l1[ll], :])
    mxe["cds_len_EI"][ll, 2]=seglen(iexn["cds_pos"][l2[ll], :])
    mxe["cds_cg"]=if rownum(mxe)==0
        String[]
    else
        rowfun(mxe["cds_len_EI"]) do (x, y)
            if x==0 && y==0
                ""
            elseif mod(abs(x-y), 3)==0
                "in-frame"
            else
                "frame-shift"
            end
        end
    end
    d"gene_name, gene_typ"mxe=dtshift(mxe["geneid"], iexn["geneid"], d"gene_name, gene_typ"iexn, ("", ""))
    # jsave("valid_MXE_splicingEvents_withAnn", mxe)
    return mxe
end
   #}}
  #{{ Merge other splicing events
function merge_mxe_into_ASE(ase, mxe)
    # mxe=jload("valid_MXE_splicingEvents_withAnn")
    mxe=@rec(mxe, ($Ejcno .> 0) .& ($Ejcno_B .> 0) .& ($LIjcno .> 0) .& ($RIjcno .> 0))
    mxe["Inum"]=mxe["LInum"].+mxe["RInum"]
    mxe["Enum"]=mxe["Enum"].+mxe["Enum_B"]
    delete!(mxe, "Enum_B")

    #Filter nested mxe
    G=linkgrp(hcat(d"Ejcno, Ejcno2, LIjcno, LIjcno2, RIjcno, RIjcno2"mxe...))
    l=grpfunexp(G, mxe["Inum"], mxe["Enum"]) do x, y
        falsesbut(length(x), sortri((min.(x, y), max.(x, y)), rev=true)[1])
    end
    mxe=rec(mxe, l)
    
    # ase=jload("valid_splicingEvents_withAsType")
    N=rnum(ase)
    ase["Ejc_pos_B"]=fill(0, N, 2)
    ase["Ejcno_B"]=fill(0, N)
    ase["Ejcno_B2"]=fill(0, N)

    t=ismbr((ase["chrno"], ase["Ejc_pos"], ase["LIjc_pos"]), (mxe["chrno"], mxe["Ejc_pos"], mxe["LIjc_pos"]))
    @assert all(.!t)
    l1=memberr((ase["chrno"], ase["Ejc_pos"], ase["LIjc_pos"]), (mxe["chrno"], mxe["LIjc_pos"], mxe["Ejc_pos"]))
    delasei=Int[]
    delmxei=Int[]
    for (ai, mi) in zip(findall(l1.>0), l1[l1.>0])
        if min(ase["Enum"][ai], ase["Inum"][ai])<=min(mxe["Enum"][mi], mxe["Inum"][mi])
            push!(delasei, ai)
        else
            push!(delmxei, mi)
        end
    end
    ase=rec_i(ase, delasei)
    mxe=rec_i(mxe, delmxei)

    t=ismbr((ase["chrno"], ase["Ejc_pos"], ase["LIjc_pos"]), (mxe["chrno"], mxe["Ejc_pos"], mxe["LIjc_pos"]))
    @assert all(.!t)
    l2=memberr((ase["chrno"], ase["Ejc_pos"], ase["RIjc_pos"]), (mxe["chrno"], mxe["Ejc_pos_B"], mxe["RIjc_pos"]))
    delasei=Int[]
    delmxei=Int[]
    for (ai, mi) in zip(findall(l2.>0), l2[l2.>0])
        if min(ase["Enum"][ai], ase["Inum"][ai])<=min(mxe["Enum"][mi], mxe["Inum"][mi])
            push!(delasei, ai)
        else
            push!(delmxei, mi)
        end
    end
    ase=rec_i(ase, delasei)
    mxe=rec_i(mxe, delmxei)

    # Solve the Ejcno conflict problem
    delasei=Int[]
    delmxei=Int[]
    for (mi, ai) in zip(setxri(mxe["Ejcno"], ase["Ejcno"])...)
        mV=(mxe["Enum"][mi], mxe["LInum"][mi]+mxe["RInum"][mi])
        aV=(ase["Enum"][ai], ase["LInum"][ai]+ase["RInum"][ai])
        t=sortri([min(mV...) max(mV...); min(aV...) max(aV...)], rev=true)[1]
        if t==1
            push!(delasei, ai)
        else
            push!(delmxei, mi)
        end
    end
    mxe=rec_i(mxe, delmxei)
    ase=rec_i(ase, delasei)
    @assert(!any(ismbr(mxe["Ejcno"], ase["Ejcno"])) && all(mxe["Ejcno"].>0))
    
    T=vcatr(ase, mxe)
    # jsave("valid_splicingEvents_withAsTypeAndMXE", T)
    return T
end
#}}
#}}
#{{ Define intron retention
  #{{ Detect intron retention
function detect_intron_retention(kexn, exn)
    # exn=jload(annfn)
    eie, _=grpfunwith(exn["transid"], exn) do T
        if rownum(T)<2
            return nothing
        end
        T=rec(T, sortri(T["exn_pos"]))
        t=T["exn_pos"]
        T=rec(T, 1:rownum(T)-1)
        # T["eie_pos"]=[t[1:end-1, 1] t[2:end, 2]]
        T["itn_pos"]=[t[1:end-1, 2] t[2:end, 1]]
        fd(T, c"chrno, itn_pos, strand")
    end
    eie=tbuniq(eie)
    # l=ismbr(d"chrno, eie_pos"eie, d"chrno, exn_pos"exn)
    # eie=rec(eie, l)
    # l=grpfunexp(x->falsesbut(length(x), argmin(x)),
    #             d"chrno, itn_pos, strand"eie, seglen(eie["eie_pos"]))
    # eie=rec(fd_i(eie, c"transid, trans_name"), l)
    _, ei=genomemap(d"chrno, exn_pos, strand"exn, d"chrno, itn_pos, strand"eie)
    eie=rec(eie, unival(ei))
    
    # Remove intron touched with known exon boundary
    # kexn=jload(exnfn)
    ei, _=genomemap(d"chrno, itn_pos, strand"eie,
                    (kexn["chrno"]|>x->[x; x],
                     kexn["exn_pos"]|>x->[x[:, 1].+1; x[:, 2].-1]|>x->[x x],
                     kexn["strand"]|>x->[x; x]))
    eie=rec_i(eie, ei)

    # Add geneid
    d"gene_len, cds_len"kexn=grpfunexp(kexn["transid"], kexn["exn_pos"], kexn["cds_pos"]) do p1, p2
        (sum(seglen(p1)), sum(seglen(p2)))
    end
    intron=tb(c"transid, itn_pos"=>findintron(kexn["transid"], kexn["exn_pos"]))
    d"chrno, strand"intron=dtshift(intron["transid"], kexn["transid"], d"chrno, strand"kexn, safe=true)
    intron["itn_pos"]=intron["itn_pos"]|>x->[x[:, 1].-1 x[:, 2].+1]
    d"geneid, gene_len, cds_len, chrno, strand"intron=dtshift(intron["transid"], kexn["transid"], d"geneid, gene_len, cds_len, chrno, strand"kexn)
    eie["geneid"]=grpfun(fastgrp(d"chrno, itn_pos, strand"intron, d"chrno, itn_pos, strand"eie), intron, default="") do T
        T["geneid"][sortri((T["cds_len"], T["gene_len"]), rev=true)[1]]
    end
    l=findall(.!startswith.(eie["geneid"], "ENSG")) #check from here 30 Mar 2020
    ki, ii=genomemap(d"chrno, exn_pos, strand"kexn, getrow(d"chrno, itn_pos, strand"eie, l))
    geid, uii=grpfun(ii, rec(kexn, ki)) do T
        T["geneid"][sortri((T["cds_len"], T["gene_len"]), rev=true)[1]]
    end
    eie["geneid"][l[uii]]=geid
    
    d"gene_name, gene_typ"eie=dtshift(eie["geneid"], kexn["geneid"], d"gene_name, gene_typ"kexn, ("", ""))
    
    jc=jload(joinpath(TMPDIR, "unified_junction_info"))
    d"num, jcno"eie=dtshift(d"chrno, itn_pos, strand"eie, d"chrno, jc_pos, strand"jc, d"uni_readnum, uJcno"jc, (0, 0))
    l=grpfunexp(x->falsesbut(length(x), argmax(x)), d"chrno, itn_pos"eie, eie["num"])
    eie=rec(eie, l)
    
    #Remove nested IR
    ai, bi=genomemap(d"chrno, itn_pos"eie, d"chrno, itn_pos"eie)
    eie=rec_i(eie, ai[ai.!=bi])
    
    # jsave("all_intron_retention_ASE", eie)
    return eie
end
  #}}
  #{{ count intron retention reads eachcell
# KK(path="intron_retention_unspliced_readcount")do
function count_intron_retention_reads_in_each_cell(itn)
    itn=rec(itn, itn["jcno"].>0)
    itn["chr"]=no2chr(itn["chrno"], chr=false)
    mkpath(joinpath(TMPDIR, "intron_retention_unspliced_readcount"))
    # for (smp, bamfn) in zip(SMP["label"], SMP["bam"])
    pmapwd(SMP["label"], SMP["bam"], env=(TMPDIR,)) do smp, bamfn, TMPDIR
        outfn=joinpath(TMPDIR, "intron_retention_unspliced_readcount", "$(smp)_intronRet_unspliced_rdnum.jld2")
        isfile(outfn) && return nothing
        
        t=bamReadCount(bamfn, (itn["chr"], itn["itn_pos"][:, 1]|>x->[x x.+3], itn["strand"]), UMI="UB", smpBC="CB", revstrand=true, read_span_site=true)
        T1=tb(c"LInum, site_idx, CB"=>t[[1, 3, 4]])
        
        t=bamReadCount(bamfn, (itn["chr"], itn["itn_pos"][:, 2]|>x->[x.-3 x], itn["strand"]), UMI="UB", smpBC="CB", revstrand=true, read_span_site=true)
        T2=tb(c"RInum, site_idx, CB"=>t[[1, 3, 4]])
        T=tbx(T1, T2, key=c"site_idx, CB", auniq=true, buniq=true, fillempty=:AB, emptyA=ds(LInum=0), emptyB=ds(RInum=0))
        T["uJcno"]=itn["jcno"][T["site_idx"]]
        delete!(T, "site_idx")
        
        safesave(outfn, T)
    end
end
#}}
  #{{ Get intron counts
function get_intron_counts(itr)
    # itr=jload("all_intron_retention_ASE")
    fd_i!(itr, c"eie_pos, exn_pos, exn_rank, exnid, protid")
    rename!(itr, "num"=>"Enum", "jcno"=>"Ejcno")
    itr["LInum"]=fill(0, rnum(itr))
    itr["RInum"]=fill(0, rnum(itr))
    for lb in SMP["label"]
        rd=jload(joinpath(TMPDIR, "intron_retention_unspliced_readcount", "$(lb)_intronRet_unspliced_rdnum.jld2"))
        urd=tbuniq(rd, "uJcno", "LInum"=>sum, "RInum"=>sum, trim=true)
        t=dtshift(itr["Ejcno"], urd["uJcno"], d"LInum, RInum"urd, (0, 0), safe=true)
        itr["LInum"].+=t[1]
        itr["RInum"].+=t[2]
    end
    # mxe=rec(mxe, (mxe["Enum"].+mxe["Enum_B"].>=5) .& (mxe["LInum"].+mxe["RInum"].>=5) .& (0.01 .< (mxe["LInum"].+mxe["RInum"])./(mxe["Enum"].+mxe["Enum_B"].+mxe["LInum"].+mxe["RInum"]) .< 0.99))
    itr["Inum"]=max.(itr["LInum"], itr["RInum"])
    #In the version before 20200915, above sentence has a bug as: itr["Inum"]=max(itr["LInum"], itr["RInum"])
    itr=rec(itr, (itr["Enum"].>=5) .& (itr["Inum"].>=5) .& (0.01 .< itr["Inum"]./(itr["Inum"].+itr["Enum"]) .< 0.99))
    # jsave("valid_intron_retention_ASE", itr)
    return itr
end
  #}}
  #{{ Annotate intron retention
function annotate_intron_retention(itr, exn)
    # exn=jload(exnfn)

    exn=rec(exn, sortri(d"transid, exn_pos"exn))
    exn["hasCdsDwstream"]=grpfunexp(exn["transid"], exn["cds_pos"], exn["strand"]) do cP, strand
        k=findfirst(cP[:, 1].>0)
        if isnothing(k)
            false
        elseif strand[1]=='+'
            falsesbut(length(strand), 1:k)
        else
            @assert strand[1]=='-'
            falsesbut(length(strand), k:length(strand))
        end
    end

    cds=rec(exn, exn["cds_pos"][:, 1].>0)
    intron=tb(c"transid, itn_pos"=>findintron(exn["transid"], exn["exn_pos"]))
    d"chrno, strand"intron=dtshift(intron["transid"], exn["transid"], d"chrno, strand"exn, safe=true)
    intron["itn_pos"]=intron["itn_pos"]|>x->[x[:, 1].-1 x[:, 2].+1]
    cdsitn=tb(c"transid, itn_pos"=>findintron(cds["transid"], cds["cds_pos"]))
    d"chrno, strand"cdsitn=dtshift(cdsitn["transid"], cds["transid"], d"chrno, strand"cds, safe=true)
    cdsitn["itn_pos"]=cdsitn["itn_pos"]|>x->[x[:, 1].-1 x[:, 2].+1]

    # itr=jload("valid_intron_retention_ASE")
    for x in c"AS_typ, cds_cg"
        itr[x]=fill("", rnum(itr))
    end
    itr["isKnownEI"]=falses(rnum(itr), 2)
    is_jc_anned=ismbr(d"chrno, itn_pos, strand"itr, d"chrno, itn_pos, strand"cdsitn)
    _, t1=genomemap(d"chrno, cds_pos, strand"cds,
                    (itr["chrno"], itr["itn_pos"][:, 1]|>x->[x x], itr["strand"]))
    _, t2=genomemap(d"chrno, cds_pos, strand"cds,
                    (itr["chrno"], itr["itn_pos"][:, 2]|>x->[x x], itr["strand"]))
    is_jc_in_cds=falsesbut(rnum(itr), intersect(t1, t2))
    _, t=genomemap(d"chrno, cds_pos, strand"cds, d"chrno, itn_pos, strand"itr)
    is_itn_in_cds=falsesbut(rnum(itr), t)
    _, t=genomemap(d"chrno, exn_pos, strand"exn, d"chrno, itn_pos, strand"itr)
    is_itn_in_exn=falsesbut(rnum(itr), t)

    t=grpfun(cds["transid"], cds["cds_pos"], cds["strand"]) do x, s
        L=minimum(x[:, 1])
        R=maximum(x[:, 2])
        if s[1]=='+'
            ([L L], [R R])
        else
            @assert s[1]=='-'
            ([R R], [L L])
        end
    end
    trans=tb((c"TSS, TES", "transid")=>t)
    d"chrno, strand"trans=dtshift(trans["transid"], cds["transid"], d"chrno, strand"cds)

    t, _=genomemap(d"chrno, itn_pos, strand"itr, d"chrno, TSS, strand"trans)
    is_itn_TSS=falsesbut(rnum(itr), t)
    t, _=genomemap(d"chrno, itn_pos, strand"itr, d"chrno, TES, strand"trans)
    is_itn_TES=falsesbut(rnum(itr), t)

    # is_itn_termcode=getseq(d"chrno, itn_pos, strand"itr, path="/user/leuven/323/vsc32366/projects/biodata/genome/hg38/hg38_jld") do seq
    #     in(AA_Term, translate(LongDNASeq(seq[1:div(length(seq), 3)*3])))
    # end

    ei, ii=genomemap(d"chrno, exn_pos, strand"exn, d"chrno, itn_pos, strand"itr)
    t, uii=grpfun(ii, ei) do l
        any(exn["hasCdsDwstream"][l])
    end
    is_cds_downstream=falses(rnum(itr))
    is_cds_downstream[uii]=t

    # d"isKnownEI, cds_cg"itr=mapr(is_jc_anned, is_jc_in_cds, is_itn_in_cds, is_itn_in_exn, is_itn_termcode, is_itn_TSS, is_itn_TES, is_cds_downstream, rem.(seglen(itr["itn_pos"]), 3).==0) do jc_anned, jc_in_cds, itn_in_cds, itn_in_exn, itn_termcode, itn_TSS, itn_TES, cds_downstream, div3
    d"isKnownEI, cds_cg"itr=mapr(is_jc_anned, is_jc_in_cds, is_itn_in_cds, is_itn_in_exn, is_itn_TSS, is_itn_TES, is_cds_downstream, rem.(seglen(itr["itn_pos"]), 3).==0) do jc_anned, jc_in_cds, itn_in_cds, itn_in_exn, itn_TSS, itn_TES, cds_downstream, div3
        if jc_anned
            if itn_in_cds
                ([true true], "in-frame", [0 0])
            else
                ck=itn_in_exn ? [true true] : [true false]
                if jc_in_cds
                    if div3
                        if itn_TES
                            (ck, "tail-cut")
                        elseif itn_TSS
                            (ck, "head-cut")
                        # elseif itn_termcode
                        #     (ck, "pre-term")
                        # else
                            #     (ck, "in-frame")
                        else
                            (ck, "in-frame-or-pre-term")
                        end
                    else
                        (ck, "frame-shift")
                    end
                elseif itn_TSS
                    if cds_downstream #here is a leak: if a second CDS starts in the "half-exon" just downstream of the retented-trion, this situation will be misassigned to 'nonsense' but it is actually be a "head-cut".
                        (ck, "head-cut")
                    else
                        (ck, "nonsense")
                    end
                elseif itn_TES
                    (ck, "tail-cut")
                else
                    (ck, "")
                end
            end
        else
            if itn_in_cds
                if div3
                    ([false true], "in-frame")
                else
                    ([false true], "frame-shift")
                end
            elseif itn_TSS
                if cds_downstream #here is a leak: if a second CDS starts in the "half-exon" just downstream of the retented-trion, this situation will be misassigned to 'nonsense' but it is actually be a "head-cut".
                    ([false true], "head-cut")
                else
                    ([false true], "nonsense")
                end
            elseif itn_TES
                ([false true], "tail-cut")
            elseif itn_in_exn
                ([false true], "")
            else
                ([false false], "")
            end
        end
    end
    # jsave("valid_intron_retention_ASE_annotated", itr)
    return itr
end
  #}}
  #{{ Merge with other splicing events
function merge_intron_retention_into_ASE(ase, itr)
    # ase=jload("valid_splicingEvents_withAsTypeAndMXE")
    # itr=jload("valid_intron_retention_ASE_annotated")
    rename!(itr, "itn_pos"=>"Ejc_pos")
    itr=rec(itr, itr["Ejcno"].>0)
    delete!(itr, "AS_typ")
    l=ismbr(itr["Ejcno"], @with(ase, vcat($Ejcno, $Ejcno2, $Ejcno_B, $Ejcno_B2, $LIjcno, $LIjcno2, $RIjcno, $RIjcno2)))
    itr=rec(itr, .!l)
    t=ds(AS_typ="IR", RIjcno2=0, LIjcno=0, LIjc_pos=[0 0], cds_len_EI=[0 0], Ejcno_B2=0, Ejcno_B=0, RIjcno=0, Ejc_pos_B=[0 0], LIjcno2=0, RIjc_pos=[0 0], Ejcno2=0)
    dictconc!(itr, reprow(t, rnum(itr)))
    T=vcatr(ase, itr)
    # jsave("valid_splicingEvents_allAsTyp", T)
    return T
end
  #}}
#}}
#{{ Write ASE info to tsv
function write_ase_info(T)
    outfn=joinpath(OUTDIR, "detected_ASE_info.tsv")
    # isfile(outfn) && return nothing
    ff=X->begin
        rowfun(X, as_vec=true) do x
            if x[1]==0
                "NA"
            else
                f"$1-$2"(x[1], x[2])
            end
        end
    end
    M=@with T asmx("ASE_no"=>f"ASE$1".($Ejcno),
                   "gene_EnsemblID"=>$geneid,
                   "gene_name"=>$gene_name,
                   "AS_type"=>$AS_typ,
                   "coding_change"=>isempty($cds_cg) ? "NA" : $cds_cg,
                   "chromosome"=>no2chr($chrno),
                   "strand"=>$strand,
                   "E_junc_A"=>ff($Ejc_pos),
                   "E_junc_B"=>ff($Ejc_pos_B),
                   "I_junc_A"=>ff($LIjc_pos),
                   "I_junc_B"=>ff($RIjc_pos),
                   "is_E_junc_known"=>$isKnownEI[:, 1],
                   "is_I_junc_known"=>$isKnownEI[:, 2])
    writedlm(outfn*".part", M)
    mv(outfn*".part", outfn)
end
#}}
#{{ Inclusive and exclusive read number per cell
function count_reads_for_each_ASE_in_each_cell(se)
    # se=jload("valid_splicingEvents_allAsTyp")
    # for lb in SMP["label"]
    pmapwd(SMP["label"], env=(OUTDIR, TMPDIR)) do lb, OUTDIR, TMPDIR
        outfn=joinpath(OUTDIR, "$(lb)_count_perASE_perCell.tsv")
        isfile(outfn*".gz") && return nothing

        jc=jload(joinpath(TMPDIR, "parse_junctions", "$lb.junction_readnum.jld2"))
        jcD=hashgrp(jc["uJcno"], jc["cellBC"])
        rdD=Dict(eachr((jc["uJcno"], jc["cellBC"])).=>jc["uni_readnum"])
        irrd=jload(joinpath(TMPDIR, "intron_retention_unspliced_readcount", "$(lb)_intronRet_unspliced_rdnum.jld2"))
        irrd["CB"]=cut"-1".(irrd["CB"])
        irD=Dict(zip(irrd["uJcno"], irrd["CB"]).=>zip(irrd["LInum"], irrd["RInum"]))
        # rl=rowpile()
        open(outfn, "w") do fid
            println(fid, join(c"# ASE_no, cell_barcode, E_junc_count, I_junc_A_count, I_junc_B_count", "\t"))
            for T in eachr(se)
                cellBC=unival([get(jcD, T["Ejcno"], []);
                               get(jcD, T["Ejcno2"], []);
                               get(jcD, T["Ejcno_B"], []);
                               get(jcD, T["Ejcno_B2"], []);
                               get(jcD, T["LIjcno"], []);
                               get(jcD, T["LIjcno2"], []);
                               get(jcD, T["RIjcno"], []);
                               get(jcD, T["RIjcno2"], [])])
                for cc in cellBC
                    O=tb(Ejcno=T["Ejcno"],
                         ASEid=f"ASE$1"(T["Ejcno"]),
                         LIjcno=T["LIjcno"],
                         RIjcno=T["RIjcno"],
                         cellBC=cc,
                         Enum=get(rdD, (T["Ejcno"], cc), 0)+get(rdD, (T["Ejcno2"], cc), 0)+get(rdD, (T["Ejcno_B"], cc), 0)+get(rdD, (T["Ejcno_B2"], cc), 0))
                    if T["AS_typ"]=="IR"
                        d"LInum, RInum"O=get(irD, (T["Ejcno"], cc), (0, 0))
                    else
                        O["LInum"]=get(rdD, (T["LIjcno"], cc), 0)+get(rdD, (T["LIjcno2"], cc), 0)
                        O["RInum"]=get(rdD, (T["RIjcno"], cc), 0)+get(rdD, (T["RIjcno2"], cc), 0)
                    end
                    # addrow!(rl, O)
                    println(fid, join(d"ASEid, cellBC, Enum, LInum, RInum"O, '\t'))
                end
            end
        end
        sf"gzip $1"(outfn)
        # value(rl)
    end
    # T=vcatr_with("patient"=>smp["patient"], C...)
    # jsave("uniReadnum_perEvent_perCell", T)
end
#}}
#{{ Inclusive and exclusive read number per cell type
function summarize_reads_for_each_cell_group(ase)
    # for lb in SMP["label"]
    pmapwd(SMP["label"], env=(OUTDIR, CELLBC)) do lb, OUTDIR, CELLBC
        outfn=joinpath(OUTDIR, "$(lb)_count_perASE_perCellGroup.tsv")
        isfile(outfn*".gz") && return nothing
        
        cel=CELLBC[lb]
        fn=joinpath(OUTDIR, "$(lb)_count_perASE_perCell.tsv.gz")
        @assert isfile(fn)
        T=readtb(`zcat $fn`, head=c"ASEid, cellBC, N::Enum, N::LInum, N::RInum")
        T["celltype"]=dtshift(T["cellBC"], cel["cellBC"], cel["celltype"], "")
        T=rec(T, .!isempty.(T["celltype"]))
        uT=tbuniq(T, c"ASEid, celltype", "Enum"=>sum, "LInum"=>sum, "RInum"=>sum, trim=true)
        
        D=Dict(f"ASE$1".(ase["Ejcno"]).=>ase["AS_typ"])
        uT["Inum"]=map(uT["LInum"], uT["RInum"], uT["ASEid"]) do a, b, x
            if D[x]=="MXE"
                # print(".")
                a+b
            else
                max(a, b)
            end
        end
        @with uT $PSI=$Inum./($Inum.+$Enum)
        # writedlm(outfn, asmx(c"ASEid, Enum, Inum, PSI"=>uT))
        writetb(outfn, uT, "ASEid=>ASE_no, celltype=>cell_group, Enum=>E_junc_count, Inum=>I_junc_count, PSI", with_type=false)
        sf"gzip $1"(outfn)
    end
end
#}}
end
