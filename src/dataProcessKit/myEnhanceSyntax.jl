#{{ importpkg importpy @pkgfun @pypkgfun
#[[ importpkg importpy @pkgfun @pypkgfun ]]
# @pkgfun(Julia_Pkg, fun1[=>rename][, fun2[=>rename], ]...)
# @pypkgfun(Python.package, fun1[=>rename][, fun2[=>rename], ]...)
#Model=importpkg(:Julia_Pkg|"Julia_Pkg"; preloaded=false)
#Model=importpy("Python.package"; preloaded=false)
#In importpkg function, symbol input have 6 time faster speed than string input when package is loaded, but it will not hide private variable. AbstractString input will hide private variable but slower speed.
#Import the Julia or python module only at running time, and it will not reload the same package, which is different from Module=PyCall.pywrap(PyCall.pyimport("...")). If the package name is given as symbol, dot "." should be replaced by "!".
#When preloaded is true, function will just throw an error if this package is not loaded previously.
#Note that if any other module is occupied already with the same module name, these function cannot distinglish the real module but just return the wrong module. In this case, try to set another model_name to avoid confliction.
#@pkgfun|@pypkgfun imports some functions from a Julia|Python package. importpkg||importpy can only import a preloaded package if it work in non-REPL, otherwise an error will throw while using its functions. (See https://docs.julialang.org/en/stable/manual/methods/#Redefining-Methods-1). The @pkgfun|@pypkgfun will solve this problem.
#See also: @pyhelp, pyfun, callr
#Xiong Jieyi, January 24, 2015 >Aug 14, 2015>Aug 26, 2015>Aug 28, 2015, 2 Jul 2017>13 May 2019

export importpkg, importpy, @pkgfun, @pypkgfun
function importpkg(X::AbstractString; preloaded::Bool=false)
    mdnm=Symbol("warpedpkg_$X")
    if !isdefined(Main, mdnm)
        if !preloaded && !isdefined(Main, Symbol(X))
            eval(Main,Meta.parse("import $X"))
        end
        mth=join(names(getfield(Main, Symbol(X))),",")
        eval(Main,Meta.parse("""
            module $mdnm
              import $X:$mth
              export $mth
            end"""))
    end
    getfield(Main,mdnm)::Module
end
function importpkg(X::Symbol; preloaded::Bool=false)
    if !preloaded && !isdefined(Main, X)
        # preloaded && error("Package $X has not been loaded.")
        eval(Main,Meta.parse("import $X"))
    end
    getfield(Main, X)::Module
end
__loaded_python_package__=Dict{String, Any}()
function importpy(X::AbstractString; preloaded::Bool=false)
    # # >>> Until Julia 1.1.1, a few python packages (such as scanpy) cannot be loaded by pyimport directly due to the LLVM version confliction. When this problem has solved in Julia, this step should be removed. 13 Jun 2019
    # if !preloaded && !isdefined(Main, :Libdl)
    #     Libdl=importpkg(:Libdl)
    #     invokelatest(Libdl.dlopen, ENV["HOME"]*"/.julia/conda/3/lib/python3.7/site-packages/llvmlite/binding/libllvmlite.so", Libdl.RTLD_DEEPBIND)
    # end
    # # <<<
    
    if !preloaded && !isdefined(Main, :PyCall)
        eval(Main, Meta.parse("using PyCall")) #Seems only using PyCall can completely avoid core crash in o.pymember("..."). Aug 26, 2015
    end
    global __loaded_python_package__
    O=get(__loaded_python_package__, X, nothing)
    if O==nothing
        O=invokelatest(Main.PyCall.pyimport, X)
        __loaded_python_package__[X]=O
    end
    O
end
importpy(X::Symbol; kw...)=importpy(replace(string(X), "!"=>"."); kw...)
# const invokelatest=Base.invokelatest
macro pkgfun(pkg, funs...)
    C=map(funs) do eachfun
        fun=string(eachfun)
        t=split(fun," => ")
        if length(t)>1
            sfun=t[1]
            tfun=t[2]
        else
            sfun=fun
            tfun=fun
        end
        if VERSION>=v"0.7.0"
            "$tfun(args...; kwargs...)=Base.invokelatest(importpkg(:$pkg).$sfun, args...; kwargs...)"
        elseif VERSION>=v"0.6.0-rc1"
            "$tfun(args...; kwargs...)=my.invokelatest(importpkg(:$pkg).$sfun, args...; kwargs...)"
        else
            "$tfun(args...; kwargs...)=importpkg(:$pkg).$sfun(args...; kwargs...)"
        end
    end
    txt=if length(C)==1
        C[1]
    else
        "begin\n"*join(C,'\n')*"\nend"
    end

    esc(Meta.parse(txt))
end
macro pypkgfun(pkg, funs...)
    C=map(funs) do eachfun
        fun=string(eachfun)
        t=split(fun," => ")
        if length(t)>1
            sfun=t[1]
            tfun=t[2]
        else
            sfun=fun
            tfun=fun
        end
        "$tfun(args...; kwargs...)=my.invokelatest(my.invokelatest(Main.PyCall.getproperty, pYpKG_$pkg, Symbol(\"$sfun\")), args...; kwargs...)"
    end
    # txt=if length(C)==1
    #     C[1]
    # else
    txt="begin\n(pYpKG_$pkg) = importpy(\"$pkg\")\n"*join(C,'\n')*"\nend"
    # end
    esc(Meta.parse(txt))
end
#}}

#{{ int
#[[ int ]]
# Int64 = int( 1.0|"1"|"1.0"|Symbol("1")|Array )
# Convert data to Int64 format.
# It used to be a standrad function in pervious julia version but be removed later. However, other similar function float() is kept. Here I kept this function in library for the convenience.
#See also: float
#Xiong Jieyi, 22 Feb 2018 > 2 Oct 2018

int(X::AbstractArray{T}) where {T<:AbstractFloat} = AbstractArray{Int}(X)
int(X::AbstractArray{T}) where {T<:Integer} = AbstractArray{Int}(X)
int(X::Integer)=Int(X)
int(X::AbstractArray{T}) where {T<:AbstractString} = int.(X)
int(X::AbstractArray{T}) where {T<:Symbol} = int.(X)

#Below will tragger a warning.
int(X::AbstractFloat)=Int(X)
function int(X::AbstractString)
    if in('.', X)
        Int(float(X))
    else
        Base.parse(Int, X)
    end
end
int(X::Symbol)=int(string(X))
export int
#}}

#{{ @k_str protect string

#[[ @k_str @k_mstr ]]
#Syntax: k"String_without_any_interpretation" or k"""..."""
#See also: @c_str @d_str
#Xiong Jieyi, December 11, 2014

macro k_str(x)
    x
end
macro k_mstr(x)
    x
end
export @k_str, @k_mstr

#}}

#{{ ismatche startswithe etc.

#[[ containssub startswithsub endswithsub occursinfirst occursinsub ismatchfirst containsfirst startswithfirst endswithfirst ]]
# MatchedItemVector = occurinsub|startswithsub|endswithsub|containssub(StringVector, FindString)
# (MatchedFirstString|nothing, Index) = occurinfirst|ismatchfirst(Regex, StringVector)
# (MatchedFirstString|nothing, Index) = occurinfirst|containsfirst|startswithfirst|endswithfirst(StringVector, FindString)
#Return all matched strings or the first matched string for a string-vector input.
#See also: filter
#Xiong Jieyi, 9 Jul 2014 >Aug 10, 2015 >1 Nov 2018

# ismatchl(R::Regex,V::Vector{T}) where {T<:AbstractString} = BitArray(map(x->ismatch(R,x), V))
# function ismatchvec(R::Regex,V::Vector{T}) where {T<:AbstractString}
#     warn("ismatchvec will be obolished. Use ismatchl or containsl instead.")
#     ismatchl(R,V)
# end
# function ismatchvec(R::AbstractString,V::Vector{T}) where {T<:AbstractString}
#     warn("ismatchvec will be obolished. Use ismatchl or containsl instead.")
#     containsl(V,R)
# end

# containsl(V::Vector{T},R::AbstractString) where {T<:AbstractString} = BitArray(map(x->contains(x,R), V))
# beginswithl(V::Vector{T},R::AbstractString) where {T<:AbstractString} = BitArray(map(x->startswith(x,R), V)) #For compatiable
# startswithl(V::Vector{T},R::AbstractString) where {T<:AbstractString} = BitArray(map(x->startswith(x,R), V))
# endswithl(V::Vector{T},R::AbstractString) where {T<:AbstractString} = BitArray(map(x->endswith(x,R), V))
function containssub(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    O=T[]
    for x in V
        if contains(x,R)
            push!(O,x)
        end
    end
    O
end
function beginswithsub(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    O=T[]
    for x in V
        if startswith(x, R)
            push!(O,x)
        end
    end
    O
end
function startswithsub(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    O=T[]
    for x in V
        if startswith(x, R)
            push!(O,x)
        end
    end
    O
end
function endswithsub(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    O=T[]
    for x in V
        if endswith(x, R)
            push!(O,x)
        end
    end
    O
end
function ismatchfirst(R::Regex,V::AbstractVector{T}) where {T<:AbstractString}
    for (i,x) in enumerate(V)
        if ismatch(R, x)
            return (x,i)
        end
    end
    return (nothing, 0)
end
function containsfirst(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    for (i,x) in enumerate(V)
        if contains(x,R)
            return (x,i)
        end
    end
    return (nothing, 0)
end
function beginswithfirst(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    for (i,x) in enumerate(V)
        if startswith(x, R)
            return (x,i)
        end
    end
    return (nothing, 0)
end
function startswithfirst(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    for (i,x) in enumerate(V)
        if startswith(x, R)
            return (x,i)
        end
    end
    return (nothing, 0)
end
function endswithfirst(V::AbstractVector{T},R::Union{AbstractString, AbstractChar}) where {T<:AbstractString}
    for (i,x) in enumerate(V)
    if endswith(x, R)
            return (x,i)
        end
    end
    return (nothing, 0)
end
function occursinsub(A, B::AbstractVector{T}) where {T<:AbstractString}
    l=occursin.(A, B)
    B[l]
end
function occursinfirst(A, B::AbstractVector{T}) where {T<:AbstractString}
    for (i, x) in enumerate(B)
        if occursin(x, B)
            return (x, i)
        end
    end
    return (nothing, 0)
end
export containssub, beginswithsub, startswithsub, endswithsub, ismatchfirst, containsfirst, beginswithfirst, startswithfirst, endswithfirst, occursinsub, occursinfirst

#}}

#{{ @c_str
#[[ @c_str ]]
#c" Astr , Bstr , Cstr " or c"""..."""
#Return ["Astr","Bstr","Cstr"]. $var can also be replaced.
#c"Astr#Bstr#Cstr#"c have the same result. Here, always use the last charactor as seperator.
# c"""
#   line1
#   line2
#   ...
# """c  #Return ["line1", "line2", ...] 
#The blanks at the flanks of each item will be striped.
#See also: my.char my.ismatch, @d_str, @k_str, @f_str, @cut_str
#Xiong Jieyi, 13 May 2014>9 Sep 2014>December 11, 2014>January 6, 2015>Oct 19, 2015>Oct 19, 2015>May 13, 2016>Aug 22, 2016>5 Sep 2022

export @c_str
macro c_str(S,P="")
    if isempty(P)
        D=','
    else
        @assert(P=="c", "Only allow c\"Astr,Bstr\" or c\"Astr#Bstr#\"c.")
        if ismatch(r"\n[\ \t]+$",S)
            D='\n'
            S=rstrip(rstrip(S,' '),'\t')
        end
        D=S[end]
        S=String(S[1:end-1])
    end
    C=map(split(S,D)) do x
        replace(String(strip(x)), '"'=>"\\\"")
    end

    esc(Meta.parse("[\""*join(C,"\",\"")*"\"]"))
    
    ## Replaced Meta.parse() to Expr: 2 Sep 2022
    # C=map(C) do S
    #     p=1
    #     K=Any[]
    #     while !isnothing((m=match(r"(?<!\\)\$(?:([A-Za-z\_]\w*)|\(([A-Za-z\_]\w*)\))", S, p);))
    #         d1, d2=m.captures
    #         d=!isnothing(d1) ? d1 : d2
            
    #         push!(K, S[p:(m.offset-1)])
    #         push!(K, Symbol(d))
    #         p=m.offset+length(m.match)
    #     end
    #     if isempty(K)
    #         S
    #     else
    #         push!(K, S[p:end])
    #         Expr(:string, K...)
    #     end
    # end
    # esc(Expr(:vect, C...))
end
#}}

#{{ @cols_str
#[[ @cols_str ]]
# [A1, A2, ...], [B1, B2, ...], ... = cols"""[;; control_params...;]A1, B1; A2, B2; ... """
# A grammar candy to write aligned vectors. Input will be parsed by readdlm(..., ',', String, ';') and then be stripped.
# The control_params can be:
#  hcat( #col1:#col2, $col3:$col4, ... ) #Merge columns. It must be at the beginning of control_params. The rest control parameters is based on the merged column index rather than the original one.
#  #col[:#col]=Type e.g. 2=Int, 3:5=Float64. To assigning data type of given column. Supports all numeric, Char, Bool, Symbol. The unassigned data are strings.
# Below changes are work only in text but not in control_params.
#  c=','|' '|\t|TAB      To assign column delimeter.
#  r=;|\n|ENT            To assign end-of-line delimeter.
#  q='|\0|NO             To assign quote charactor. Quotes are only used for including blanks at the begin or end, not for escaping delimeter charactors. Quote charactors inside text will not be stripped. e.g. ' Rock'n'roll '. Note that to assign comma or blank, need quoted by single quotation mark. NO means no quotes.
# This macro parse data using readdlm(;quotes=true), followed by strip() for each item. Therefore ,"...", can be used to escape delimeters but cannot protect the flanking blanks unless ,"'...'", is nested (e.g., ," a,b ",=>"a,b", ,"' a;b' ",=>" a;b"). However, when q=", readdlm(;quotes=false) will be performed and delimeters inside "..." will no longer be escaped.
# Jieyi Xiong, 6 Oct 2022 > 11 Jul 2023 >20 Mar 2024>20 Mar 2024

export @cols_str
macro cols_str(X::String)
    catn1=Int[]
    catn2=Int[]
    kys=Int[]
    funs=Function[]
    dm=','
    eof=';'
    qut='\''
    if length(X)>=2 && X[1:2]==";;"
        i=findnext(';', X, 3)
        P=X[3:i-1]
        X=X[i+1:end]
        P=lstrip(P)
        
        if startswith(P, "hcat(")
            pn=findfirst(P, ')')
            if isnothing(pn)
                error("Invalid ;;hcat(1:2, 3:4, ...); format.")
            end
            for p in split(P[5:pn-1], ',')
                m=match(r"\s*(\d+)\s*\:\s*(\d+)\s*", p)
                if isnothing(m)
                    error("Invalid ;;hcat(1:2, 3:4, ...); format.")
                end
                push!(catn1, int(m[1]))
                push!(catn2, int(m[2]))
            end
            P=strip(lstrip(P[pn+1:end], ','))
        end
        if !isempty(P)
            for p in split(P, r"\,(?!\')")
                ky, sval=split(p, '=')
                ky=strip(ky, ' ')
                sval=unescape_string(strip(sval, ' '))
                if length(sval)>2 && sval[1]=='\'' && sval[end]=='\''
                    sval=sval[2:end-1]
                end

                kym=match(r"^(\d+)(?:\s*\:(\d+)\s*)?$", ky)
                if !isnothing(kym)
                    cfunt=eval(Meta.parse(sval))::Type
                    cfun=if cfunt<:Number
                        x->parse(cfunt, x)
                    elseif cfunt<:Char
                        only
                    else
                        x->cfunt(x)
                    end
                    cky1=parse(Int, kym[1])
                    cky2=isnothing(kym[2]) ? cky1 : parse(Int, kym[2])
                    for i=cky1:cky2
                        push!(kys, i)
                        push!(funs, cfun)
                    end
                else
                    cc::Union{Char, Nothing}=if sval=="ENT"
                        '\n'
                    elseif sval=="TAB"
                        '\t'
                    elseif sval=="NO"
                        nothing
                    elseif length(sval)==1
                        sval[1]
                    else
                        Meta.parse("'$sval'")
                    end
                    if ky=="c"
                        dm=cc::Char
                    elseif ky=="r"
                        eof=cc::Char
                    elseif ky=="q"
                        qut=cc
                    else
                        error("Unsupported command `$ky'.")
                    end
                end
            end
        end
    end
    X=replace(unescape_string(X), r"\;\s*$"=>"") #Remove the last semicolon.
    O=map(eachcol(readdlm(IOBuffer(X), dm, String, eof; quotes=qut!='"'))) do xs
        map(xs) do x
            x=String(strip(x))
            if !isnothing(qut) && startswith(x, qut) && endswith(x, qut)
                x=x[2:end-1]
            end
            x
        end
    end
    
    if !isempty(catn1)
        ldel=falses(length(O))
        O=convert(Vector{Any}, O)
        for (i, j) in zip(catn1, catn2)
            if any(ldel[i:j])
                error("Overlapped hcat columns are not allowed.")
            end
            O[i]=hcat(O[i:j]...)
            ldel[i+1:j].=true
        end
        deleteat!(O, ldel)
    end
    
    if !isempty(kys)
        O=convert(Vector{Any}, O)
        for (ky, cfun) in zip(kys, funs)
            O[ky]=cfun.(O[ky])
        end
    end
    Tuple(O)
end
#}}

#{{ @f_str @F_str @cf_str @cF_str
#[[ @f_str @F_str @cf_str @cF_str ]]
# f"\$1=$1 and $(2)-$3 $(4=.3g|x|X) $(5=h)"("aa",100,200,pi,1000) will return string
#  $1=aa and 100-200 3.14 1k
# f"...$(1=...[<|>][n-width|$n][m(def)|t|e|w[fill_char=' ']])..." : Fixed width replacement **
# f"..." actually return a function. For example,
# map(f"$1-$2", 1:3, 4:6) or f"$1-$2".(1:3, 4:6) outputs:
# 3-element Array{ASCIIString,1}:
#  "1-4"
#  "2-5"
#  "3-6"
# A clear way to form a string with variable. f"""..."""(...) are also acceptable.
# f"$(N=...)" allows to specify a printf-like format after = ('%' mark should not be included), or human-readable number ($(N=h), transfer-by saynum()), but format a variable one more time is not supported. (Updated in 5 Aug 2021: f"$(1=.3g)"(NaN) return "NaN" instead of triggering an error like @sprintf(), so does Inf.)
# In format $(1=.3x|X), number will be transfered to LaTex way of scientific notation. The difference between `x' and `X' is only `X' will warp string inside $...$.
# In f"$1-$x" but not in F"$1-$x", the $x will also be replaced.
#
# ** Fixed width replacement: fillin string (or other data) in fixed width.
# e.g. f"|$(1=<5)|"("A") => "|A    |"
#      f"|$(1=>5m&)|"("A") => "|&&&&A|"
# The width can also be an argument. Note only `$2' but `$(2)' is supported.
# e.g. f"d$(1=>$2m0)b"(3, 4) => "d0003b"
#m|t|e|w: In case string is longer than width, m--ignore; t--trim the end; e--throw error; w--trim the end with a warning. "m" is the default value.
# Consider using fontname="monospace" in this mode to visually align the charactors.
#
# cf"..."[c] and cF"..."[c] : combine f"..." and c"..."[c] string macros.
# See also: @c_str, @cut_str, @k_str
# Xiong Jieyi, Aug 23, 2015 >Nov 4, 2015 >Aug 4, 2016 >17 May 2018 > 23 Aug 2022>28 May 2024

export @f_str, @F_str, @cf_str, @cF_str
# function _f_str(T::String, M::String; keep_dollor::Bool=false)
#     M1, M2=if M=="then_c"
#         ("c", "")
#     elseif M=="then_cc"
#         ("c", "c")
#     else
#         @assert isempty(M)
#         ("", "")
#     end
#     D=Set{String}() #Just for error promote.
#     C=map(eachmatch(r"(?<!\\)\$\((\d+)\=([\d\w\.\-\+\#\*]+)\)", T)) do m
#         d,f=m.captures
#         if in(d, D)
#             error("Token #$d has been repeatly formated, which is not supported so far.")
#         end
#         push!(D, d)
#         if f=="h"
#             "__insidE_valuE_dO_noT_usE__[$d]=saynum(__insidE_valuE_dO_noT_usE__[$d])"
#         elseif length(f)>1 && (f[end]=='x' || f[end]=='X')
#             ff=f[1:end-1]
#             if f[end]=='X'
#                 "__insidE_valuE_dO_noT_usE__[$d]=(__insidE_valuE_dO_noT_usE__[$d]|>x->isnan(x) || isinf(x) ? x : replace(@sprintf(\"%$(ff)g\", x), r\"(\\-?\\d+(?:\\.\\d+)?)e\\+?(\\-)?(?:0*)(\\d+)\"=>s\"\$\\1\\\\times10^{\\2\\3}\$\"))"
#             else
#                 "__insidE_valuE_dO_noT_usE__[$d]=(__insidE_valuE_dO_noT_usE__[$d]|>x->isnan(x) || isinf(x) ? x : replace(@sprintf(\"%$(ff)g\", x), r\"(\\-?\\d+(?:\\.\\d+)?)e\\+?(\\-)?(?:0*)(\\d+)\"=>s\"\\1\\\\times10^{\\2\\3}\"))"
#             end
#         else
#             "__insidE_valuE_dO_noT_usE__[$d]=(__insidE_valuE_dO_noT_usE__[$d]|>x->isnan(x) || isinf(x) ? x : @sprintf(\"%$f\", x))"
#         end
#     end
#     I=isempty(C) ? "" : "(__insidE_valuE_dO_noT_usE__=Any[__insidE_valuE_dO_noT_usE__...];"*join(C,";")*";__insidE_valuE_dO_noT_usE__)|>__insidE_valuE_dO_noT_usE__->"
    
#     T=replace(T, r"(?<!\\)\$(\d+)" => s"$(__insidE_valuE_dO_noT_usE__[\1])")
#     T=replace(T, r"(?<!\\)\$\((\d+)\)" => s"$(__insidE_valuE_dO_noT_usE__[\1])")
#     T=replace(T, r"(?<!\\)\$\((\d+)\=[\d\w\.\-\+\#\*]+\)" => s"$(__insidE_valuE_dO_noT_usE__[\1])")

#     T=replace(T, '"' => "\\\"")
#     if keep_dollor
#         T=replace(T, r"(?<!\\)\$(?!\(__insidE_valuE_dO_noT_usE__\[\d+\]\))" => "\\\$")
#     end
#     Meta.parse("(__insidE_valuE_dO_noT_usE__::Union{AbstractString, Symbol, AbstractChar, Number, Type}...)->$I$M1\"$T\"$M2")
# end

#isesc: whether to replace $var
#isparse: whether to parse \t, \n, etc. When isparse=false, only \$ will be parsed.
#halfout: output x->"...$x..." or only "...$x...".
function _f_str(S::AbstractString; isesc::Bool=false, halfout::Bool=false, isparse::Bool=true)
    p=1
    K=Any[]

    function addto(part::AbstractString, isTheLast::Bool=false)
        if ismatch(r"(?:[^\\]|^)(?:\\\\)*\$\([\d\.\s]*\)", part)
            error("Invalid format inside parentheses of \$(1).")
        elseif ismatch(r"(?:[^\\]|^)(?:\\\\)*\$(?:\([^\)]*)?\b__insidE_valuE_dO_noT_usE__\b", part)
            error("Containing \$__insidE_valuE_dO_noT_usE__ in @f_str. This is a coding trick.")
        end        
        if isesc
            @assert(isparse)
            t=Meta.parse('"'*replace(part, '"'=>"\\\"")*'"')
            if isa(t, AbstractString)
                push!(K, t)
            else
                @assert(t.head==:string)
                append!(K, t.args)
            end
        else
            part=replace(part, r"([^\\]|^)((?:\\\\)*)\\\$(?=(?:\d|\(\d))"=>s"\1\2$") #\$ will be parsed in any case.
            if isparse
                part=unescape_string(part)
            else
                @assert(!occursin(r"\\\$(?:\d|\(\d)", part) && (isTheLast || isempty(part) || part[end]!='\\'), "So far multiple anti-slashes before \\\$1 is not supported. In the future it will have the same beheave as raw\"\\\\\\\"\".")
            end
            push!(K, part)
        end
    end

    # while !isnothing((m=match(r"(?:[^\\]|^|\b)(?:\\\\)*(\$)(?:(\d+)|\((\d+)\)|\((\d+)\=([\d\w\.\-\+\#\*]+)\))", S, p);))
    while !isnothing((mSp=match(r"(?:[^\\]|^|\b)(?:\\\\)*(\$)(?:(\d+)|\((\d+)\)|\((\d+)\=([^\)]+)\))", S[p:end]);)) #If using match(r"...", S, p), there will be a bug: when parse f"$(1)$2", the "$2" under p=5 will not match this Regex. 28 May 2024
        d0, d1, d2, d3, f0=mSp.captures
        @assert d0=="\$"
        dollaridx=p+mSp.offsets[1]-1
        d=!isnothing(d1) ? d1 : !isnothing(d2) ? d2 : d3
        if r"^\d+$"(d)
            d=parse(Int, d)
        end

        addto(S[p:(dollaridx-1)], false)

        f, f1, isleftali=if !isnothing(f0)
            if in('<', f0)
                (split(f0, '<', limit=2)..., true)
            elseif in('>', f0)
                (split(f0, '>', limit=2)..., false)
            else
                (f0, nothing, false)
            end
        else
            (f0, nothing, false)
        end
        cK=if isnothing(f) || isempty(f)
            if isa(d, Int)
                :(__insidE_valuE_dO_noT_usE__[$d])
            else
                @assert(d!="__insidE_valuE_dO_noT_usE__")
                Symbol(d)
            end
        elseif f=="h"
            :(saynum(__insidE_valuE_dO_noT_usE__[$d]))
        elseif length(f)>1 && (f[end]=='x' || f[end]=='X')
            spf="%"*f[1:end-1]*"g"
            if f[end]=='X'
                :(__insidE_valuE_dO_noT_usE__[$d] |> x->isnan(x) || isinf(x) ? x : replace(@sprintf($spf, x), r"(\-?\d+(?:\.\d+)?)e\+?(\-)?(?:0*)(\d+)"=>s"$\1\\times10^{\2\3}$"))
            else
                :(__insidE_valuE_dO_noT_usE__[$d] |> x->isnan(x) || isinf(x) ? x : replace(@sprintf($spf, x), r"(\-?\d+(?:\.\d+)?)e\+?(\-)?(?:0*)(\d+)"=>s"\1\\times10^{\2\3}"))
            end
        else
            spf="%$f"
            :(__insidE_valuE_dO_noT_usE__[$d] |> x->isnan(x) || isinf(x) ? x : @sprintf($spf, x))
        end
        if !isnothing(f1)
            # m2=match(r"^(\d+)([mtwe](\S)?)?$", f1)
            m2=match(r"^(?:(\d+)|\$(\d+))([mtwe](\S)?)?$", f1)
            isnothing(m2) && error("Invalid format: `$f1'")
            swidth=if !isnothing(m2[1])
                int(m2[1])
            else
                t=int(m2[2])
                :(__insidE_valuE_dO_noT_usE__[$t]::Integer)
            end
            typ=Symbol(something(m2[3], "m")[1])
            fillchr=only(something(m2[4], " "))
            chkfun=if typ==:m
                :(x->x)
            elseif typ==:t
                :(x->length(x)>$swidth ? x[1:$swidth] : x)
            elseif typ==:w
                quote
                    x->begin
                        swidth=$swidth
                        if length(x)>swidth
                            @warn("String is longer than $swidth: $x. It was trimmed.")
                            x[1:swidth]
                        else
                            x
                        end
                    end
                end
            elseif typ==:e
                quote
                    x->begin
                        swidth=$swidth
                        if length(x)>swidth
                            error("String is longer than $swidth: $x.")
                        end
                        x
                    end
                end
            else
                error("WTF!")
            end
            fillfun=if isleftali
                :(x->x*(length(x)<$swidth ? repeat($fillchr, $swidth-length(x)) : ""))
            else
                :(x->(length(x)<$swidth ? repeat($fillchr, $swidth-length(x)) : "")*x)
            end
            cK=:($fillfun($chkfun(string($cK))))
        end
        
        push!(K, cK)
        p=p+mSp.offset+length(mSp.match)-1
    end

    addto(S[p:end], true)
    
    hfexpr=Expr(:string, K...)
    if halfout
        hfexpr
    else
        :((__insidE_valuE_dO_noT_usE__...)->$hfexpr)
    end
end
function  _f_str(V::AbstractVector{<:AbstractString}; isesc::Bool=false, halfout::Bool=false)
    C=_f_str.(V, isesc=isesc, halfout=true)
    hfexpr=Expr(:vect, C...)
    if halfout
        hfexpr
    else
        :((__insidE_valuE_dO_noT_usE__...)->$hfexpr)
    end
end
macro f_str(T)
    esc(_f_str(T, isesc=true))
end
macro F_str(T)
    _f_str(T, isesc=false)
end
macro cF_str(S, P="")
    if isempty(P)
        D=','
    else
        @assert(P=="c", "Only allow c\"Astr,Bstr\" or c\"Astr#Bstr#\"c.")
        if ismatch(r"\n[\ \t]+$",S)
            D='\n'
            S=rstrip(rstrip(S,' '),'\t')
        end
        D=S[end]
        S=String(S[1:end-1])
    end
    C=map(x->String(strip(x)),split(S,D))
    _f_str(C, isesc=false)
end
macro cf_str(S, P="")
    if isempty(P)
        D=','
    else
        @assert(P=="c", "Only allow c\"Astr,Bstr\" or c\"Astr#Bstr#\"c.")
        if ismatch(r"\n[\ \t]+$",S)
            D='\n'
            S=rstrip(rstrip(S,' '),'\t')
        end
        D=S[end]
        S=String(S[1:end-1])
    end
    C=map(x->String(strip(x)),split(S,D))
    esc(_f_str(C, isesc=true))
end
#}}

#{{ @cut_str
#[[ @cut_str ]]
# cut"DNJ"(string) = join( split(string, D)[N], J)
# cut"DN"(string) = join( split(string, D)[N], D)
# cut"N"(string) = split(string)[N]
# Get a part of string. N could be a number, end, range (with end), and a vector.
# e.g. cut"-2"("a-b-c-d") = "b"
#  cut"-end"("a-b-c-d") = "d"
#  cut"-2:end"("a-b-c-d") = "b-c-d"
#  cut"-2:end+"("a-b-c-d") = "b+c+d"
#  cut"-1, 3"("a-b-c-d") = "a-c"
# D cannot be any digits or "end". J can not be any digits or alphabates. Both D and J support transfer charactors like \t, but $ should not be wrote as \$.
# See also: @c_str, @cc_str, @f_str, @F_str, @k_str
# Xiong Jieyi, 15 Oct 2018

export @cut_str
# macro cut_str(S)
#     m=match(r"([^\d]+)(\d+)$", S)
#     m==nothing && error("Invalid cut\"DN\" grammar.")
#     D=m.captures[1]
#     N=int(m.captures[2])
#     if N==0
#         x->split(x, D)[end]
#     else
#         x->split(x, D)[N]
#     end
# end
macro cut_str(S)
    S=Meta.parse("\""*replace(S, "\$" => "\\\\\\\$")*"\"")
    m=match(r"\d|end", S)
    m==nothing && error("Invalid cut\"DNJ\" grammar.")
    J=""
    if m.offset>1
        D=S[1:m.offset-1]
        J=D
    end
    R=S[m.offset:end]
    mR=match(r"[^\d\w]+$", R)
    if mR!=nothing
        J=mR.match
        R=R[1:mR.offset-1]
    end
    if in(',', R)
        R="[$R]"
    end
    o=if m.offset>1
        if in(':', R) || in(',', R)
            "x->join(split(x, \"$D\")[$R], \"$J\")"
        else
            "x->String(split(x, \"$D\")[$R])"
        end
    else
        if in(':', R) || in(',', R)
            isempty(J) && error("The join element is not specified.")
            "x->join(split(x)[$R], \"$J\")"
        else
            "x->String(split(x)[$R])"
        end
    end
    Meta.parse(o)
end

#}}

#{{ str2char

#[[ str2char ]]
#Char_vector=str2char("AbstractString"[, width])
#Char_matrix=str2char(String_Vector[, width])
#Convert string to char vector. The end will be filled with ' '.
#See also: ascii, utf8, rstrip, cellstr, @k_str
#Xiong Jieyi, 9 Jul 2014

str2char(X::AbstractString)=[x for x in X]
function str2char(X::AbstractString,N::Integer)
    C=fill(' ',N)
    for i=1:length(X)
        C[i]=X[i]
    end
end
#<v0.6# function str2char{T<:AbstractString}(X::Vector{T})
function str2char(X::Vector{T}) where {T<:AbstractString}
    N=maximum([length(x) for x in X])
    str2char(X,N)
end
#<v0.6# function str2char{T<:AbstractString}(X::Vector{T},N::Int)
function str2char(X::Vector{T},N::Int) where {T<:AbstractString}
    C=fill(' ',length(X),N)
    for ri=1:length(X)
        ctxt=X[ri]
        for ci=1:length(ctxt)
            C[ri,ci]=ctxt[ci]
        end
    end
    return C
end
export str2char

# #[[ char2str ]]
# #AbstractString=char2str(Char_vector) %Blanks in the end will be trimmed.
# #String_matrix=char2str(Char_matrix)
# #See also: cellstr, str2char, @k_str, rstrip
# #Xiong Jieyi, 11 Jul 2014

# function char2str(X::Vector{Char})
#     i=length(X)
#     while i>0 && X[i]==' '
#         i-=1
#     end
#     return string(X[1:i]...)
# end
# function char2str(X::Matrix{Char})
#     S=fill("",size(X,1))
#     for i=1:size(X,1)
#         S[i]=char2str(X[i,:])
#     end
#     return S
# end
# char2str

#}}

#{{ cellstr

#[[ cellstr ]]
# StringVector=cellstr(Char_Mx|Dict)
#The same function as cellstr in Matlab.
#If input is a Dict, run cellstr to every char_matrix.
#See also: ascii, utf8, rstrip, str2char, tidy, tidytype
#Xiong Jieyi, 3 Sep 2014

#<v0.6# function cellstr{T<:Char}(X::Matrix{T})
function cellstr(X::Matrix{T}) where {T<:Char}
    S=String(X'[:])
    L=size(X,2)
    return [rstrip(S[((i-1)*L+1):i*L]) for i=1:size(X,1)]
end
cellstr(X)=X
cellstr(X::Dict)=dictfun(cellstr,X)
export cellstr

#}}

#{{ num2str
#[[ num2str ]]
# AbstractString|StringArray|Tuple = num2str(Real|RealArray|Tuple; mindigit=N|fixdigit=false, maxdec=M, fillblanks=false, base=10)
#Convert number matrix to string matrix. The ending 0 in decimal part will be removed (i.e., 1.0=>"1" rather than "1.0").
#If mindigit is given, the number-string with length less than N will be filled with 0 or blanks on the left, depending on fillblanks=T|F.
#If maxdec is given, the decimal part will be trimed to M digitals.
#If fixdigit is ture (First input must be a array), mindigit will be set by the digital of maximum number automatically.
#base=... is the same as string(::Int; base=...). Only avaiable when input is a integer (array).
#See also: str2char, cellstr
#Xiong Jieyi, 27 Aug 2014>January 12, 2015>Feb 8, 2016>8 Nov 2018>30 Jan 2024

export num2str
function num2str(X::Real; mindigit::Integer=1, maxdec::Integer=-1, fillblanks::Bool=false, kw...)
    if maxdec>=0
        X=round(X, digits=maxdec)
    end
    if isinteger(X)
        S=string(int(X); kw...)
    else
        S=string(X; kw...)
    end
    if mindigit>1 && length(S)<mindigit
        if fillblanks
            S=repeat(' ', mindigit-length(S))*S
        else
            if startswith(S, '-')
                S=S[1:1]*repeat('0', mindigit-length(S))*S[2:end]
            else
                S=repeat('0', mindigit-length(S))*S
            end
        end
    end
    return S
end
function num2str(X::AbstractArray{<:Real}; fixdigit::Bool=false, fillblanks::Bool=false, kw...)
    isempty(X) && return Array{String}(undef, size(X))
    # if fixdigit
    #     mindigit=floor(Int, log(base, maximum(X))+1)
    # end
    C=map(x->num2str(x; kw...), X)
    if fixdigit
        maxlen=maximum(length, C)
        for i in 1:length(C)
            S=C[i]
            if length(S)<maxlen
                if fillblanks
                    C[i]=repeat(' ', maxlen-length(S))*S
                else
                    if startswith(S, '-')
                        C[i]=S[1:1]*repeat('0', maxlen-length(S))*S[2:end]
                    else
                        C[i]=repeat('0', maxlen-length(S))*S
                    end
                end
            end
        end
    end
    return C
end
num2str(X::Tuple; wargs...)=tuple(num2str([X...]; wargs...)...)
#}}

#{{ saynum
#[[ saynum ]]
#Syntax: [AbstractString] = mysay(Any_number)
#Or  String_arr|tuple = mysay(array|tuple of number)
#Change the N to human readable format.
#N should be a scalar. If no input and output is given, ans will be as input.
#For the SI prefix, please see http://en.wikipedia.org/wiki/SI_prefix
#5 Aug 2021: NaN and Inf will no longer trigger an error.
#See also: num2str
#Xiong Jieyi, February 8, 2015>Jun 1, 2015>5 Aug 2021

export saynum
function saynum(N::Real)
    if isnan(N) || isinf(N)
        return string(N)
    end
    if N==0
	return "0"
    elseif N<0
	N=-N
	sg="-"
    else
	sg=""
    end
    k=floor(log10(N))
    if k<3 && k>=0
	S=string(sg,num2str(N))
    elseif k>0
	unin="kMGTPEZY"
	cu=floor(Int,k/3)
	N=round(Int,N/10^(k-3+1))
	N=N/10^(2-mod(k,3))
	S=*(sg,num2str(N),unin[[cu]])
    else
	unin="munpfazy"
	cu=floor(Int,k/3)
	N=round(Int,N/10^(k-3+1))
	N=N/10^(2-mod(k,3))
	S=*(sg,num2str(N),unin[[-cu]])
    end
    S
end
saynum(V::Union{AbstractArray,Tuple})=map(saynum,V)
#}}

#{{ truesbut falsesbut
#[[ truesbut falsesbut ]]
# V=truesbut|falsesbut(N|(Nrow,Ncol), M...)
# Return a BitArray V, where V=trues|falses(N) but V[M...]=false|true.
#See also:
#Xiong Jieyi, March 27, 2015 > 8 Feb 2019

export truesbut, falsesbut
function truesbut(N, M...)
    O=trues(N)
    O[M...].=false
    O
end
function truesbut(N, M::Integer...)
    O=trues(N)
    O[M...]=false
    O
end
function falsesbut(N, M...)
    O=falses(N)
    O[M...].=true
    O
end
function falsesbut(N, M::Integer...)
    O=falses(N)
    O[M...]=true
    O
end
#}}

#{{ typein typein_dict
#typein(::Array{T})=T
#typein_dict(Dict)=(T_key, T_value)
#Return the type of the element of array or dict.
#See also: deany
#Xiong Jieyi, 2 Sep 2014

typein(X::Array{T}) where {T} =T
typein_dict(X::Dict{T1,T2}) where {T1,T2} =(T1,T2)
# export typein, typein_dict
#}}

#{{ tidynumbit, tidystr, tidyboolarr, tidytype, tidy

#[[ tidynumbit tidystr tidyboolarr tidytype tidy ]]
#tidied_X=tidynumbit|tidystr|tidyboolarr|tidy(X)
#Make tidy the types.
#tidynumbit: Convert any integer to Int, any float number to AbstractFloat.
#tidystr: Convert SubString to AbstractString, Convert @compat(Union{T<:AbstractString...}) to a proporate AbstractString type. It will dramatically enhance jldsave speed and it is default for jldsave.
#tidyboolarr: bitpack bool array.
#tidytype(X) = tidyboolarr(tidynumbit(tidystr(X)))
#tidy(X) = tidytype(col2vec(X)) #Also convert column-matrix to vector.
#All these function will do recurly for Tuple and Dict input, and output the input if the convertion is not needed.
#See also: col2vec, cellstr, retype, jldsave
#Xiong Jieyi, 28 Aug 2014 >May 20, 2015

export tidy, tidytype, tidynumbit, tidystr, tidyboolarr

tidynumbit(X::Bool)=X
tidynumbit(X::Char)=X
tidynumbit(X::Integer)=int(X)
tidynumbit(X::Array{Bool})=X
tidynumbit(X::Array{Char})=X
tidynumbit(X::Array{T}) where {T<:Integer} =int(X)
tidynumbit(X::AbstractFloat)=float(X)
tidynumbit(X::Array{T}) where {T<:AbstractFloat} =float(X)
tidynumbit(X::Tuple)=map(tidynumbit,X)
tidynumbit(X::Dict)=dictfun(tidynumbit,X)
tidynumbit(X)=X

tidystr(X::Array{T}) where {T<:String} =X

function tidystr(X::Array{T}) where {T<:AbstractString}
    try
        convert(Array{ASCIIString},X)
    catch
        try
            convert(Array{UTF8String},X)
        catch
            try
                convert(Array{UTF16String},X)
            catch
                convert(Array{UTF32String},X)
            end
        end
    end
end
tidystr(X::SubString{T}) where {T} =tidystr(convert(T,X))
tidystr(X::AbstractArray{SubString{T}}) where {T} =tidystr(convert(AbstractArray{T},X))
tidystr(X::Tuple)=map(tidystr,X)
tidystr(X::Dict)=dictfun(tidystr,X)
tidystr(X)=X

# tidytype(X::AbstractString)=ascii(X)
# # tidytype(X::@compat(Union{Array,Tuple}))=map(tidytype,X)
# tidytype{T<:AbstractString}(X::Array{T})=convert(Array{ASCIIString},X)
# tidytype(X::Tuple)=map(tidytype,X)
# tidytype(X::Dict)=dictfun(tidytype,X)
tidyboolarr(X::Array{Bool})=BitArray(X)
tidyboolarr(X::Tuple)=map(tidyboolarr,X)
tidyboolarr(X::Dict)=dictfun(tidyboolarr,X)
tidyboolarr(X)=X

tidytype(X)=tidyboolarr(tidynumbit(tidystr(X)))
tidy(X)=tidytype(col2vec(X))

#}}

#{{ desubstr
#[[ desubstr ]]
# String = desubstr(SubString)
# String_Array = desubstring(SubString_Array)
#Convert SubString to normal string type.
#See also: jldsave, tidystr
#Xiong Jieyi, Nov 29, 2015

export desubstr
desubstr(X::SubString{T}) where {T}=convert(T,X)
desubstr(X::AbstractArray{SubString{T}}) where {T}=convert(AbstractArray{T},X)
#}}

#{{ jldsave jldload jldlook jldreadme jldmeta
#[[ jldsave jldload jldload2 jldlook jldreadme jldmeta ]]
#jldsave(filename[, Dict][, "name1"=>val1, "name2"=>val2,...];
#        allow_any = false,
#        meta=Dict("mvar1"=>val1, "mvar2"=>val2,...), readme="...",
#        del=c"var1,...", istidystr=true, append=false, rewrite=false)
#
#Dict=jldload(filename[, varnames]; showinfo=isinteractive())
#(Value1,Value2,...)=jldload(filename,varname1,varname2,...; ...)
#
#Dict=jldmeta(filename[, meta_varnames_vec]; noreadme=false)
#Value1, Value2, ...=jldmeta(filename, meta_varname1, meta_varname2,...)
#Dict, MetaDict = jldload2(filename; noreadme=false) #No meta-tip show in screen.
#Dict, MetaVal1, MetaVal2, ... = jldload2(filename, meta_varname1, meta_varname2, ...)
#jldload_output, MetaVal1, MetaVal2, ... = jldload2((filename, ...), meta_varname1, meta_varname2, ...)
#
#varnames = jldlook(filename; showinfo=isinteractive())
#readme_str = jldreadme(filename)
#
#Saving and load data in JLD format. Filename will be added .jld automatically.
#When allow_any = false, saving Array{Any} data is not allowed (It's slow to save a Array{Any}, and in most cases it is not necessary.).
#When istidystr is true, tidystr will be used to convert the input data.
#When append is true, new variable or meta could be added in the old file.
#When rewrite is true, new variable with the same name of old variable will be rewritted. Otherwise a error is occur. rewrite=true imply append=true.
#When del or delmeta is given, the given variable will be removed before saving.
#For large ASCIIString data, JLD format are much faster than MAT format. Note that string should be non-Union and non-SubString type.
#Update in Jun 14, 2015: Support "/" in keyname.
#If want to save in JLD2 format, try jsave jload jload2 jlook jreadme jmeta.
#See also: ds, tidy, tidytype, matload, matsave, matlook
#Xiong Jieyi,November 23, 2014>December 11, 2014>May 23, 2015>May 27, 2015>Jun 14, 2015>Nov 29, 2015>Feb 3, 2016>Feb 3, 2016>8 Mar 2017>31 Dec 2017

export jldsave, jldload, jldmeta, jldlook, jldreadme, jldload2

__trslash(S)=replace(S, r"/" => "\\%")
__reslash(S)=replace(S, r"\\%" => "/")
function jldsave(filename::AbstractString, data::Dict{T1,}=Dict{ByteString,Any}(); istidystr::Bool=true, readme::Union{AbstractString, Nothing}=nothing, meta::Dict{T2,}=Dict{ByteString,Any}(),append::Bool=false,rewrite::Bool=false,del::Vector{T3}=ASCIIString[],delmeta::Vector{T4}=ASCIIString[],allow_any::Bool=false) where {T1<:ByteString,T2<:ByteString,T3<:ByteString,T4<:ByteString}
    isempty(delmeta) || error("jldsave(...; delmeta=...) is no more supported.") #Added in 8 Sep 2020
    
    @pkgfun(JLD, jldopen=>JLD_jldopen, o_delete=>JLD_o_delete, exists=>JLD_exists, g_create=>JLD_g_create)
    if !allow_any
        is_any=isa.(collect(values(data)), Array{Any})
        if any(is_any)
            error(f"Variable $1 is Array{Any}. Set allow_any=true."(collect(keys(data))[findfirst(is_any)]))
        end
        is_any=isa.(collect(values(meta)), Array{Any})
        if any(is_any)
            error(f"Meta variable $1 is Array{Any}. Set allow_any=true."(collect(keys(meta))[findfirst(is_any)]))
        end
    end
    
    
    if !ismatch(r"\.jld$",filename)
        filename*=".jld"
    end
    if istidystr
        data=tidystr(data)
    end
    if !isnothing(readme) && !isempty(readme)
        meta=convert(Dict{ByteString,Any},meta)
        meta["readme"]=readme
    # elseif rewrite && readme==""
    #     delmeta=[delmeta; "readme"]
    end
    if !isempty(del) || !isempty(delmeta)
        rewrite=true
    end
    if rewrite
        append=true
    end
    JLD_jldopen(filename,append ? "r+" : "w") do fid
        for k in del
            JLD_o_delete(fid,__trslash(k))
        end
        for (k,v) in data
            k=__trslash(k)
            if rewrite && JLD_exists(fid,k)
                JLD_o_delete(fid,k)
            end
            fid[k]=v
        end
        if !isempty(meta)
            if !JLD_exists(fid,"__metadata__")
                mtgrp=JLD_g_create(fid,"__metadata__")
            else
                mtgrp=fid["__metadata__"]
                for k in delmeta
                    JLD_o_delete(mtgrp,__trslash(k))
                end
            end
            for (k,v) in meta
                k=__trslash(k)
                if rewrite && JLD_exists(mtgrp,k)
                    JLD_o_delete(mtgrp,k)
                end
                mtgrp[k]=v
            end
        end
    end
end
jldsave(filename::AbstractString,T::Dict;wargs...)=jldsave(filename,retype_key(T,:check,ByteString);wargs...)#For [Any=>Any] table input.
jldsave(filename::AbstractString,arg1::Dict,args::Pair...;wargs...)=jldsave(filename,ds(arg1,args...);wargs...)
jldsave(filename::AbstractString,args::Pair...;wargs...)=jldsave(filename,ds(args...);wargs...)

function jldload(filename::AbstractString,varnames::Vector{Tv}=ASCIIString[]; showinfo::Bool=isinteractive()) where {Tv<:ByteString}
    @pkgfun(JLD, jldopen=>JLD_jldopen, exists=>HDF5_exists)
    # @pkgfun(HDF5, exists=>HDF5_exists)
    
    if !ismatch(r"\.jld$",filename)
        filename*=".jld"
    end
    readme=""
    T=JLD_jldopen(filename, "r") do fid
        if showinfo && HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD.JldGroup)
                if HDF5_exists(fid,"__metadata__/readme")
                    readme=read(fid,"__metadata__/readme")
                end
                mtnm=setdiff(names(fid["__metadata__"]),c"readme")
                if !isempty(mtnm)
                    readme=readme*"\n[Meta data] "*join(mtnm,",")
                end
            else
                readme=read(fid,"__metadata__")["readme"]
            end
        else
            showinfo=false
        end
        if isempty(varnames)
            varnames=names(fid)
            varnames=varnames[varnames.!="__metadata__"]
            # Dict{ByteString,Any}(__reslash(var) => _read(fid, var) for var in varnames)
            TT=Dict{ByteString,Any}()
            for var in varnames
                V=try
                    read(fid, var)
                catch err
                    @warn(string(err))
                    println("Field $var is ignored.")
                    undef
                end
                if !(V===undef)
                    TT[__reslash(var)]=V
                end
            end
            TT
        else
            Dict{ByteString,Any}(var => read(fid, __trslash(var)) for var in varnames)
        end
    end
    if showinfo
        @info(readme)
    end

    if VERSION<v"0.7"
        T=dictfun!(T, replace=true) do v
            if isa(v, Char) && v!='\0' && (UInt32(v) & 0x00ffffff)==0
                Char(UInt32(v) >>> 24)
            elseif isa(v, Array) && eltype(v)<:Char && all(x->(UInt32(x) & 0x00ffffff)==0, v)
                map(x->Char(UInt32(x) >>> 24), v)
            else
                v
            end
        end
    end
    T
end
function jldload(filename::AbstractString, varname1::ByteString, varname2::ByteString, varnames::ByteString...)
    T=jldload(filename,[varname1, varname2, varnames...])
    map(x->T[x],tuple(varname1, varname2, varnames...))
end
jldload(filename::AbstractString, varname::ByteString)=jldload(filename,[varname])[varname]
function jldlook(filename::ByteString; showinfo::Bool=isinteractive())
    @pkgfun(JLD, jldopen=>JLD_jldopen, exists=>HDF5_exists)
    # @pkgfun(HDF5, exists=>HDF5_exists)
    
    if !ismatch(r"\.jld$",filename)
        filename*=".jld"
    end    
    F=JLD_jldopen(filename, "r") do fid
        if showinfo && HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD.JldGroup)
                if HDF5_exists(fid,"__metadata__/readme")
                    readme=read(fid,"__metadata__/readme")
                else
                    readme=""
                end
                mtnm=setdiff(names(fid["__metadata__"]),c"readme")
                if !isempty(mtnm)
                    readme=readme*"\n[Meta data] "*join(mtnm,",")
                end
            else
                readme=read(fid,"__metadata__")["readme"]
            end
            @info(readme)
        end
        map(__reslash,names(fid))
    end
    F[F.!="__metadata__"]
end
function jldmeta(filename::AbstractString, varnames::Vector{Ts}=ByteString[]; noreadme::Bool=false) where {Ts<:ByteString}
    @pkgfun(JLD, jldopen=>JLD_jldopen, exists=>HDF5_exists)
    # @pkgfun(HDF5, exists=>HDF5_exists)
    
    if !ismatch(r"\.jld$",filename)
        filename*=".jld"
    end
    T=JLD_jldopen(filename, "r") do fid
        if HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD.JldGroup)
                mtgrp=fid["__metadata__"]
                if isempty(varnames)
                    varnames=names(mtgrp)
                    if noreadme
                        varnames=varnames[varnames.!="readme"]
                    end
                    Dict{ByteString,Any}(__reslash(k)=>read(mtgrp,k) for k in varnames)
                else
                    Dict{ByteString,Any}(k=>read(mtgrp,__trslash(k)) for k in varnames)
                end

            else
                meta=read(fid,"__metadata__")
                if !isempty(varnames)
                    fd(meta,varnames)
                end
            end
        else
            Dict{ByteString,Any}()
        end
    end
    
    if VERSION<v"0.7"
        T=dictfun!(T, replace=true) do v
            if isa(v, Char) && v!='\0' && (UInt32(v) & 0x00ffffff)==0
                Char(UInt32(v) >>> 24)
            elseif isa(v, Array) && eltype(v)<:Char && all(x->(UInt32(x) & 0x00ffffff)==0, v)
                map(x->Char(UInt32(x) >>> 24), v)
            else
                v
            end
        end
    end
    T
end
function jldmeta(filename::AbstractString, varname1::ByteString, varname2::ByteString, varnames::ByteString...)
    T=jldmeta(filename,[varname1, varname2, varnames...])
    map(x->T[x],tuple(varname1, varname2, varnames...))
end
jldmeta(filename::AbstractString, varname::ByteString)=jldmeta(filename,[varname])[varname]
function jldreadme(filename)
    @pkgfun(JLD, jldopen=>JLD_jldopen, exists=>HDF5_exists)
    # @pkgfun(HDF5, exists=>HDF5_exists)
    if !ismatch(r"\.jld$",filename)
        filename*=".jld"
    end
    T=JLD_jldopen(filename, "r") do fid
        if HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD.JldGroup)
                mtgrp=fid["__metadata__"]
                if HDF5_exists(mtgrp, "readme")
                    mtgrp["readme"]
                else
                    nothing
                end
            else
                nothing
            end
        else
            nothing
        end
    end
end
function jldload2(filename::AbstractString, args...; kw...)
    O1=jldload(filename, showinfo=false)
    O2=jldmeta(filename, args...; kw...)
    if length(args)>1
        tuple(O1, O2...)
    else
        (O1, O2)
    end
end
function jldload2(arg1::Tuple, args...; kw...)
    O1=jldload(arg1..., showinfo=false)
    O2=jldmeta(arg1[1], args...; kw...)
    if length(args)>1
        tuple(O1, O2...)
    else
        (O1, O2)
    end
end
#}}

#{{ jsave jload jlook jreadme jmeta
#[[ jsave jload jload2 jlook jreadme jmeta ]]
#jsave(filename[, Dict][, "name1"=>val1, "name2"=>val2,...];
#        allow_any = false,
#        meta=Dict("mvar1"=>val1, "mvar2"=>val2,...) | :unchanged, readme="...",
#        istidystr=true, append=false)
#jsave(...; del=c"var1,...", delmeta=c"mvar1,...", rewrite=T|F, ...) #Slow!!
#Dict=jload(filename[, varnames]; showinfo=isinteractive())
#(Value1,Value2,...)=jload(filename,varname1,varname2,...; ...)
#Dict=jmeta(filename[, meta_varnames_vec]; noreadme=false)
#Value1, Value2, ...=jmeta(filename, meta_varname1, meta_varname2,...)
#Dict, MetaDict = jload2(filename; noreadme=false) #No meta-tip show in screen.
#Dict, MetaVal1, MetaVal2, ... = jload2(filename, meta_varname1, meta_varname2, ...)
#jldload_output, MetaVal1, MetaVal2, ... = jload2((filename, ...), meta_varname1, meta_varname2, ...)
#varnames = jlook(filename; showinfo=isinteractive())
#readme_str|nothing = jreadme(filename)
#
#Saving and load data in JLD2 format. Filename will be added .jld2 automatically.
#In jsave, when meta=:unchange, the meta will be firstly read in the same file by jmeta.
#When allow_any = false, saving Array{Any} data is not allowed (It's slow to save a Array{Any}, and in most cases it is not necessary.).
#When istidystr is true, tidystr will be used to convert the input data.
#When append is true, new variable or meta could be added in the old file.
#When rewrite is true, new variable with the same name of old variable will be rewritted. When del=c"..." or delmeta=c"..." is given, these data or metadata will be removed from file. Since JLD2 do not support remove or rewrite fields, actually all the data will be loaded and resaved again, so thess funcitons cost time when file is large.
#When del or delmeta is given, the given variable will be removed before saving.
#See also: jldsave jldload jldlook jldreadme jldmeta, ds, tidy, tidytype, matload, matsave, matlook
#Xiong Jieyi, 3 Sep 2018 > 8 Sep 2020

export jsave, jload, jlook, jmeta, jload2, jreadme

#__trslash(S)=replace(S, r"/" => "\\%") #Already defined in jld... functions.
#__reslash(S)=replace(S, r"\\%" => "/")
function jsave(filename::AbstractString, data::Dict{T1,}=Dict{String,Any}(); istidystr::Bool=true, readme::Union{AbstractString, Nothing}=nothing, meta::Union{Dict{T2,}, Symbol}=Dict{String,Any}(),append::Bool=false,rewrite::Bool=false,del::Vector{T3}=String[],delmeta::Vector{T4}=String[],allow_any::Bool=false) where {T1<:String,T2<:String,T3<:String,T4<:String}
    if VERSION<v"0.7.0" && any(x->isa(x, Char), values(data))
        @warn("JLD2 saved Char data in Julia 0.6 will not be properly loaded in Julia 0.7.")
    end
    @pkgfun(JLD2, jldopen=>JLD_jldopen, haskey=>JLD_exists, Group=>JLD_g_create)
    JLD_o_delete(t, x)=error("Delete variable $x is not supported by JLD2 so far.") # This function could be redefined in the future once JLD2 supports delete fields. 8 Sep 2020
    if !ismatch(r"\.jld2$",filename)
        filename*=".jld2"
    end
    if isa(meta, Symbol)
        meta==:unchanged || error("meta should be a Dict or :unchanged.")
        meta=jmeta(filename)::Dict
    end
    if !allow_any
        is_any=isa.(collect(values(data)), Array{Any})
        if any(is_any)
            error(f"Variable $1 is Array{Any}. Set allow_any=true."(collect(keys(data))[findfirst(is_any)]))
        end
        is_any=isa.(collect(values(meta)), Array{Any})
        if any(is_any)
            error(f"Meta variable $1 is Array{Any}. Set allow_any=true."(collect(keys(meta))[findfirst(is_any)]))
        end
    end
    if istidystr
        data=tidystr(data)
    end
    if !isnothing(readme) && !isempty(readme)
        meta=convert(Dict{String,Any},meta)
        meta["readme"]=readme
    end
    if !isempty(del) || !isempty(delmeta)
        rewrite=true
    end
    if rewrite
        # append=true
        println("jsave rewrite mode will first load out data then resave it.")
        println("Loading...")
        oldfds=jlook(filename, showinfo=false)
        olddata, oldmeta=jload2((filename, setdiff!(oldfds, [del..., keys(data)...])))
        fd_i!(oldmeta, delmeta)
        if readme==""
            delete!(oldmeta, "readme")
        end
        print("Saving...")
        O=jsave(filename, merge!(olddata, data); meta=merge!(oldmeta, meta), istidystr=istidystr, readme=readme, allow_any=allow_any)
        println(" Done.")
        return O
    end
    JLD_jldopen(filename,append ? "r+" : "w") do fid
        for k in del
            JLD_o_delete(fid,__trslash(k))
        end
        for (k,v) in data
            k=__trslash(k)
            if rewrite && JLD_exists(fid,k)
                JLD_o_delete(fid,k)
            end
            fid[k]=v
        end
        if !isempty(meta)
            if !JLD_exists(fid,"__metadata__")
                mtgrp=JLD_g_create(fid,"__metadata__")
            else
                mtgrp=fid["__metadata__"]
                for k in delmeta
                    JLD_o_delete(mtgrp,__trslash(k))
                end
            end
            for (k,v) in meta
                k=__trslash(k)
                if rewrite && JLD_exists(mtgrp,k)
                    JLD_o_delete(mtgrp,k)
                end
                mtgrp[k]=v
            end
        end
    end
end
jsave(filename::AbstractString,T::Dict;wargs...)=jsave(filename,retype_key(T,:check,String);wargs...)#For [Any=>Any] table input.
jsave(filename::AbstractString,arg1::Dict,args::Pair...;wargs...)=jsave(filename,ds(arg1,args...);wargs...)
jsave(filename::AbstractString,args::Pair...;wargs...)=jsave(filename,ds(args...);wargs...)

function jload(filename::AbstractString,varnames::Vector{Tv}=String[]; showinfo::Bool=isinteractive()) where {Tv<:String}
    if !ismatch(r"\.jld2$",filename)
        filename*=".jld2"
    end
    isfile(filename) || error("Cannot find file: $filename")
    @pkgfun(JLD2, jldopen=>JLD_jldopen, keys=>JLD_keys, haskey=>HDF5_exists)
    readme=""
    T=JLD_jldopen(filename, "r") do fid
        if showinfo && HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD2.Group)
                if HDF5_exists(fid,"__metadata__/readme")
                    readme=read(fid,"__metadata__/readme")
                end
                mtnm=setdiff(JLD_keys(fid["__metadata__"]),c"readme")
                if !isempty(mtnm)
                    readme=readme*"\n[Meta data] "*join(mtnm,",")
                end
            else
                readme=read(fid,"__metadata__")["readme"]
            end
        else
            showinfo=false
        end
        if isempty(varnames)
            varnames=JLD_keys(fid)
            varnames=varnames[varnames.!="__metadata__"]
            Dict{String,Any}(__reslash(var) => read(fid, var) for var in varnames)
        else
            Dict{String,Any}(var => read(fid, __trslash(var)) for var in varnames)
        end
    end
    if showinfo
        @info(readme)
    end
    T
end
function jload(filename::AbstractString, varname1::String, varname2::String, varnames::String...; kw...)
    T=jload(filename,[varname1, varname2, varnames...]; kw...)
    map(x->T[x],tuple(varname1, varname2, varnames...))
end
jload(filename::AbstractString, varname::String)=jload(filename,[varname])[varname]
function jlook(filename::String; showinfo::Bool=isinteractive())
    if !ismatch(r"\.jld2$",filename)
        filename*=".jld2"
    end
    isfile(filename) || error("Cannot find file: $filename")
    @pkgfun(JLD2, jldopen=>JLD_jldopen, keys=>JLD_keys, haskey=>HDF5_exists)
    F=JLD_jldopen(filename, "r") do fid
        if showinfo && HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD2.Group)
                if HDF5_exists(fid,"__metadata__/readme")
                    readme=read(fid,"__metadata__/readme")
                else
                    readme=""
                end
                mtnm=setdiff(JLD_keys(fid["__metadata__"]),c"readme")
                if !isempty(mtnm)
                    readme=readme*"\n[Meta data] "*join(mtnm,",")
                end
            else
                readme=read(fid,"__metadata__")["readme"]
            end
            @info(readme)
        end
        map(__reslash, JLD_keys(fid))
    end
    F[F.!="__metadata__"]
end
function jmeta(filename::AbstractString, varnames::Vector{T}=String[]; noreadme::Bool=false) where {T<:String}
    if !ismatch(r"\.jld2$",filename)
        filename*=".jld2"
    end
    isfile(filename) || error("Cannot find file: $filename")
    @pkgfun(JLD2, jldopen=>JLD_jldopen, keys=>JLD_keys, haskey=>HDF5_exists)
    JLD_jldopen(filename, "r") do fid
        if HDF5_exists(fid,"__metadata__")
            if isa(fid["__metadata__"], Main.JLD2.Group)
                mtgrp=fid["__metadata__"]
                if isempty(varnames)
                    varnames=JLD_keys(mtgrp)
                    if noreadme
                        varnames=varnames[varnames.!="readme"]
                    end
                    Dict{String,Any}(__reslash(k)=>mtgrp[k] for k in varnames)
                else
                    Dict{String,Any}(k=>mtgrp[__trslash(k)] for k in varnames)
                end

            else
                meta=read(fid,"__metadata__")
                if !isempty(varnames)
                    fd(meta,varnames)
                end
            end
        else
            Dict{String,Any}()
        end
    end
end
function jmeta(filename::AbstractString, varname1::String, varname2::String, varnames::String...)
    T=jmeta(filename,[varname1, varname2, varnames...])
    map(x->T[x],tuple(varname1, varname2, varnames...))
end
jmeta(filename::AbstractString, varname::String)=jmeta(filename,[varname])[varname]
function jreadme(filename)
    # jmeta(filename, "readme")
    if !ismatch(r"\.jld2$",filename)
        filename*=".jld2"
    end
    isfile(filename) || error("Cannot find file: $filename")
    @pkgfun(JLD2, jldopen=>JLD_jldopen, keys=>JLD_keys, haskey=>HDF5_exists)
    JLD_jldopen(filename, "r") do fid
        if HDF5_exists(fid,"__metadata__") && isa(fid["__metadata__"], Main.JLD2.Group)
            mtgrp=fid["__metadata__"]
            if HDF5_exists(mtgrp, "readme")
                mtgrp["readme"]
            else
                nothing
            end
        else
            nothing
        end
    end
end
function jload2(filename::AbstractString, args...; kw...)
    O1=jload(filename, showinfo=false)
    O2=jmeta(filename, args...; kw...)
    if length(args)>1
        tuple(O1, O2...)
    else
        (O1, O2)
    end
end
function jload2(arg1::Tuple, args...; kw...)
    O1=jload(arg1..., showinfo=false)
    O2=jmeta(arg1[1], args...; kw...)
    if length(args)>1
        tuple(O1, O2...)
    else
        (O1, O2)
    end
end
#}}

#{{ pmapwd
#[[ pmapwd ]]
# ...=pmapwd((P1,P2,...,env1,evn2,...)->fun, P1, P2, ...; env=(env1, env2, ...), dir=pwd(), wargs...)
#Run pmap with the repeated env variables.
# ...=pmapwd(...; debug=true) #Run single process model for debug.
# ...=pmapwd(...; tryone=true) #Run one process first, and then run remains parallelly.
# nothing =pmapwd(...; hasout=false, ...) #Give up results collection.
#Note that env should always be a tuple. For the one input, write env=(V,) rather env=(V).
#The env=(...) parameter is not necessary since pmap can directly visit environment variables in Julia 1.0.
#The different beheaves between pmapwd() and pmap(): pmapwd() will cd to the "latest" current path when pmapwd() was called, but pmap() will stick on the "original" current path when the whole script was started.
#See also: @nohup, retype
#Xiong Jieyi, 3 Sep 2014>October 24, 2014>November 11, 2014>January 9, 2015>January 11, 2015>March 3, 2015>March 3, 2015>March 19, 2015>19 Feb 2019>13 Mar 2024

#The reason I don't export __pmapwith__(): this function have unstable beheave about the current path, depending on the "real" or "fake" multi-core mode.
function __pmapwith__(fun::Function, args...; env::Tuple=(), tryone::Bool=false, debug::Bool=false, wargs...)
    if debug
        @warn("pmapwd: Run single core in debug model.")
    elseif nprocs()==1
        @warn("Only one Julia core is running.")
        debug=true
    end
    pnum=length(args[1])
    repenv=map(x->fill(x,pnum),env)
    if debug || pnum<=1
        map(fun,args...,repenv...)
    elseif tryone
        println("Run first task on main CPU.")
        arg1=map(first, args)
        repenv1=map(first, repenv)
        out1=fun(arg1...,repenv1...)
        println("Run other $(pnum-1) tasks on $(nprocs()-1) CPUs parallelly.")
        argr=map(x->Iterators.drop(x, 1), args)
        repenvr=map(x->Iterators.drop(x, 1), repenv)
        outr=pmap(fun,argr...,repenvr...; wargs...)
        [out1, outr...]
    else
        pmap(fun,args...,repenv...; wargs...)
    end
end
function pmapwd(fun::Function, args...; env::Tuple=(), dir::AbstractString=pwd(), hasout::Bool=true, wargs...)
    pfun=if hasout
        (x...)->cd(()->fun(x[1:end-1]...),x[end])
    else
        (x...)->(cd(()->fun(x[1:end-1]...),x[end]); nothing)::Nothing
    end
    pout=__pmapwith__(pfun, args...; env=tuple(env..., dir), wargs...)
    if hasout
        pout
    else
        nothing
    end
end
export pmapwd
#}}

#{{ cmd runstr successstr
#[[ cmd runstr successstr ]]
#Cmd = cmd("linux_stript_as_str")
#runstr("linux_stript_as_str")
#T|F=successstr("linux_stript_as_str")
#Run linux script as string. Pipeline "|" in the str is allowed.
#See also: prun, flow
#Xiong Jieyi, January 9, 2015 > Aug 10, 2015

export runstr, successstr, cmd
cmd(S::AbstractString)=`bash -c $S`
function runstr(S::AbstractString)
    #method 1
    # writeall(`bash`|>stdout,S)

    #method 2
    # fn=tempname()
    # writeall(fn,S)
    # run(`chmod u+x $fn`)
    # run(`$fn`)
    # rm(fn)

    #method 3
    run(cmd(S))
    nothing    
end
function successstr(S::AbstractString)
    # fn=tempname()
    # writeall(fn,S)
    # run(`chmod u+x $fn`)
    # rlt=success(`$fn`)
    # rm(fn)
    # rlt
    success(cmd(S))
end

#}}

#{{ @sh_str @sf_str @shc_str @sfc_str

#[[ @sh_str @sf_str @shc_str @sfc_str ]]
#sh|shc"Astr,Bstr,Cstr" or sh|shc"""..."""
#sf|sfc"... $1 ... $2 ..."(Val1, Val2, ...)
#Run the string using runstr (sh and sf), or return a command (shc and sfc).
#sh-shell; sf-shell function; shc-shell command; sfc-shell function command.
#The "\" charactor will not be translated, different from Julia string.
#sf|sfc"..." will keep $ untouched except it following a number, but f"..." will parse other $val also.
#See also: runstr, prun, qsub_sh, @f_str
#Xiong Jieyi, May 15, 2015 > May 8, 2016 >May 15, 2016

export @sh_str, @shc_str, @sf_str, @sfc_str
macro sh_str(S)
    :(runstr($S))
end
macro shc_str(S)
    `bash -c $S`
end
macro sf_str(T)
    # T=replace(T, r"\\(?!\$(\d|\(\d+\)))" => "\\\\")
    # T=replace(T, '"' => "\\\"")
    # T=replace(T, r"(?<![^\\]\\)\$\d+" => x->"\$(__insidE_valuE_dO_noT_usE__[$(x[2:end])])")
    # T=replace(T, r"(?<![^\\]\\)\$\(\d+\)" => x->"\$(__insidE_valuE_dO_noT_usE__[$(x[3:end-1])])")
    # T=replace(T, r"(?<![^\\]\\)\$(?!\(__insidE_valuE_dO_noT_usE__)" => "\\\$")
    # Meta.parse("(__insidE_valuE_dO_noT_usE__...)->runstr(\"$T\")")
    t=_f_str(T, isesc=false, halfout=true, isparse=false)
    :((__insidE_valuE_dO_noT_usE__...)->runstr($t))
end
macro sfc_str(T,__insidE_valuE_dO_noT_usE__...)
    # T=replace(T, r"\\(?!\$(\d|\(\d+\)))" => "\\\\")
    # T=replace(T, '"' => "\\\"")
    # T=replace(T, r"(?<![^\\]\\)\$\d+" => x->"\$(__insidE_valuE_dO_noT_usE__[$(x[2:end])])")
    # T=replace(T, r"(?<![^\\]\\)\$\(\d+\)" => x->"\$(__insidE_valuE_dO_noT_usE__[$(x[3:end-1])])")
    # T=replace(T, r"(?<![^\\]\\)\$(?!\(__insidE_valuE_dO_noT_usE__)" => "\\\$")
    # esc(Meta.parse("(__insidE_valuE_dO_noT_usE__...)->`bash -c \$(\"$T\")`"))
    t=_f_str(T, isesc=false, halfout=true, isparse=false)
    # @assert t.head==:string
    # tt=Expr(:string, "bash -c \$(\"", t.args..., "\")")
    :((__insidE_valuE_dO_noT_usE__...)->cmd($t))
end

#}}

#{{ dict2tuple
#[[ dict2tuple ]]
#Tuple=dict2tuple(Dict[, FieldsIterater=keys(Dict)])
# (V1, (V2, V3)) = dict2tuple(Dict, ("F1", ("F2", "F3"))) #Recurser supported.
#Convert the values of a Dict to tuple. Only string-key Dict is supported.
#Recursed fields is supported since 27 Jan 2020.
#See also: fd, fd_i
#Xiong Jieyi, 7 Sep 2014 > 22 Nov 2018 > 27 Jan 2020

# dict2tuple(D::Dict, Fs)=tuple([D[f] for f in Fs]...)
function dict2tuple(D::Dict{Tk,}, Fs::Tuple) where {Tk<:AbstractString}
    map(Fs) do f
        if isa(f, AbstractString)
            D[f]
        else
            dict2tuple(D, f)
        end
    end
end
dict2tuple(D::Dict{Tk,}, Fs::AbstractVector{Tv}) where {Tk<:AbstractString, Tv<:AbstractString}=dict2tuple(D, tuple(Fs...))
dict2tuple(D::Dict{Tk,}) where {Tk<:AbstractString} = tuple(values(D)...)
export dict2tuple
#}}

#{{ promote_array promote_tuple promote_group promote_anyrow promote_anyrow!

#[[ promote_array ]]
#(A::Array{T,N},B::Array{T,N})=promote_array(A::Array{T1,N1}, B::Array{T2,N2})
#promote the dimention and element type, which based on promote_rule. The dimention only support 1~2 so far.
#See also: promote_tuple, promote_group, promote, promote_rule, promote_type, promote_anyrow!
#Xiong Jieyi, 9 Sep 2014 > 15 Jun 2021

function promote_array(A::AbstractVector{AT},B::AbstractVector{BT}) where {AT,BT}
    if AT==BT
        return (A,B)
    end
    T=promote_type(AT,BT)
    if T==AT
        B=convert(AbstractVector{T},B)
    else
        A=convert(AbstractVector{T},A)
    end
    return (A,B)
end
function promote_array(A::AbstractMatrix{AT},B::AbstractMatrix{BT}) where {AT,BT}
    @assert(size(A,2)==size(B,2),"Two matrix groups have different column number.")
    if AT==BT
        return (A,B)
    end
    T=promote_type(AT,BT)
    if T==AT
        B=convert(AbstractMatrix{T},B)
    else
        A=convert(AbstractMatrix{T},A)
    end
    return (A,B)
end
promote_array(A::AbstractVector,B::AbstractMatrix)=promote_array(repeat(A,1,1),B)
promote_array(A::AbstractMatrix,B::AbstractVector)=promote_array(A,repeat(B,1,1))

#[[ promote_tuple ]]
#((A1,A2,...), (B1,B2,...))=promote_tuple((A1,A2,...), (B1,B2,...))
#Promote each type of element in the tuple. If array, function will use promote_array instead.
#See also: promote_array, promote_group
#Xiong Jieyi, 9 Sep 2014

function promote_tuple(A::Tuple,B::Tuple)
    if typeof(A)==typeof(B)
        return (A,B)
    end
    len=length(A)
    @assert(len==length(B),"A and B haven't the same lengths.");
    nA=Array{Any}(undef, len)
    nB=Array{Any}(undef, len)
    for i=1:len
        if isa(A[i],AbstractArray)
            (nA[i],nB[i])=promote_array(A[i],B[i])
        else
            (nA[i],nB[i])=promote(A[i],B[i])
        end
    end
    return (tuple(nA...),tuple(nB...))
end

#[[ promote_group ]]
#(A,B)=promote_group(A,B)
#Promote each type of element in the group. Function will choose promote_array or promote_tuple to run depend on the types. Not support Dict group.
#See also: promote_array, promote_tuple
#Xiong Jieyi, 9 Sep 2014
promote_group(A::AbstractArray,B::AbstractArray)=promote_array(A,B)
promote_group(A::Tuple,B::Tuple)=promote_tuple(A,B)

#[[ promote_anyrow promote_anyrow! ]]
# (A',B') = promote_anyrow!(A::Dict,B::Dict)
# (A',B') = promote_anyrow(A::Dict,B::Dict)
# promote_anyrow[!](A,B)=promote_group(A,B) #with or without ! is identical.
#
#When input is Dict:
#For promote_anyrow!, the input A and B may change if promote is needed.
#For promote_anyrow, the input A and B will be deepcopied if promote is needed.
#
#See also: promote_array, promote_tuple, promote_group
#Xiong Jieyi, May 23, 2015 > 15 Jun 2021

function promote_anyrow!(A::Dict,B::Dict)
    for (k,v) in A
        # A[k],B[k]=promote_group(v,B[k])
        A[k],B[k]=promote_anyrow!(v,B[k]) #14 Jun 2021: To support nested dict in dtshift().
    end
    (A,B)
end

#14 Jun 2021: Changed this function to supported nest Dict, although so far no other function using promote_anyrow(). Just for my OCD....
# function promote_anyrow(A::Dict,B::Dict)
#     intact=true
#     for (k,v) in A
#         vB=B[k]
#         if intact && typeof(v)!=typeof(vB)
#             intact=false
#             A=deepcopy(A)
#             B=deepcopy(B)
#         end
#         A[k],B[k]=promote_group(v,vB)
#     end
#     (A,B)
# end
function promote_anyrow(A::Dict, B::Dict)
    intactA=true
    intactB=true
    for (k, vA) in A
        vB=B[k]
        pA, pB=if intactA || intactB
            promote_anyrow(vA, vB)
        else
            promote_anyrow!(vA, vB)
        end
        if pA!==vA
            if intactA
                intactA=false
                A=deepcopy(fd_i(A, k))
            end
            A[k]=pA
        end
        if pB!==vB
            if intactB
                intactB=false
                B=deepcopy(fd_i(B, k))
            end
            B[k]=pB
        end
    end
    (A, B)
end
promote_anyrow!(A,B)=promote_group(A,B)
promote_anyrow(A,B)=promote_group(A,B)

export promote_array, promote_tuple, promote_group, promote_anyrow, promote_anyrow!

#}}

#{{ splitdim splitdimv [abolished outside]
# These functions is abolished. Using eachrow|eachcol() instead.
# Tuple = splitdim(AbstractMatrix, Dim=1|;dims=1)
# Vector = splitdimv(...)
#Cut the matrix into slices or vector and return the tuple or vector of slices. Each element of output container is always a vector. splitdimv has faster speed than splitdim.
#See also: dtshift, rowfun
#Xiong Jieyi, 12 Sep 2014>October 24, 2014>October 31, 2014>Aug 31, 2015>21 Nov 2018

function splitdimv(X::AbstractMatrix, D::Integer=1; dims::Integer=D)
    if dims==1
        V=[X[:,i] for i=1:size(X,2)]
    elseif dims==2
        V=[X[i,:] for i=1:size(X,1)]
    else
        error("dims only support 1 or 2 so far.")
    end
    return V
end
splitdim(args...;wargs...)=tuple(splitdimv(args...;wargs...)...)
# export splitdim, splitdimv
#}}

#{{ op [ obsoleted ]
#=
#...=op(fun, A, B)
#Expand the row or column of A or B, in order to make then have the same size, then run fun. Input could be scalar, vector or matrix. The vector input will be regarded as column matrix.
#The difference between bradcast: the result type will be promotered if nessary, rather than output a error.
#See also: splitdim, broadcast
#Xiong Jieyi, 14 Sep 2014>December 8, 2014

function op(fun::Function,A::AbstractMatrix,B::AbstractMatrix)
    Ar,Ac=size(A)
    Br,Bc=size(B)
    if Ar==Br && Ac==Bc
        #do nothing
    elseif Ar==Br
        if Ac==1
            A=repmat(A,1,Bc)
        elseif Bc==1
            B=repmat(B,1,Ac)
        else
            error("Invalid shape. A has $Ac columns but B has $Bc columns.")
        end
    elseif Ac==Bc
        if Ar==1
            A=repmat(A,Br,1)
        elseif Br==1
            B=repmat(B,Ar,1)
        else
            error("Invalid shape. A has $Ar rows but B has $Br rows.")
        end
    else
        error("Invalid shape. A: $Ar x $Ac but B: $Br x $Bc.")
    end
    return fun(A,B)
end
op(fun::Function,A::AbstractVector,B::AbstractMatrix)=op(fun,A[:,:],B)
op(fun::Function,A::AbstractMatrix,B::AbstractVector)=op(fun,A,B[:,:])
op(::Function,::AbstractArray,::AbstractArray)=error("Invalid shape.")
function op(fun::Function,A::AbstractArray,B)
    @assert(!isa(B,AbstractArray),"Invalid shape.")
    fun(A,fill(B,size(A)))
end
function op(fun::Function,A,B::AbstractArray)
    @assert(!isa(A,AbstractArray),"Invalid shape.")
    fun(fill(A,size(B)),B)
end
# export op
=#
#}}

#{{ onlyone findonly
#[[ onlyone findonly ]]
#X[1]=onlyone(X)
# findonly(...) equal to onlyone(find(...))
#If X is not an one-element vector or tuple, a error will occur.
#See also: truesbut, falsesbut, findfirst, findlast
#Xiong Jieyi, Jul 22, 2015>Aug 6, 2015>Nov 23, 2015>21 Feb 2019

export onlyone, findonly
function onlyone(X)
    out=iterate(X)
    out==nothing && error("Input is empty.")
    X1, stat=out
    iterate(X, stat)==nothing || error("Input has more than one element.")
    X1
end
function findonly(args...)
    I=findfirst(args...)
    if VERSION>=v"0.7.0"
        (I==nothing || I==0) && error("No any TURE in vector.")
        t=findnext(args..., I+1)
        (t==nothing || t==0) || error("More than one TRUE in vector.")
    else
        I!=0 || error("No any TURE in vector.")
        findnext(args...,I+1)==0 || error("More than one TRUE in vector.")
    end
    I
end
#}}

#{{ inp
#[[ inp ]]
#(args,kwargs)=inp(var1, var2, ..., key1=var3, key2=var4, ...)
#Returne the argument list of input.
#See also:
#Xiong Jieyi, December 9, 2014

inp(args...;wargs...)=(args,wargs)
export inp
#}}

#{{ retype retype_key

#[[ retype retype_key ]]
#Array = retype(Array[, :check|:share][, Type = typeof(Array[1])])
#Dict = retype_key(Dict[, :share | :check ][, Type])
#Convert Array-element or Dict-key to given type.
# :share Convert to the shared ancester type, or using given type if it is higher than all element types.
# :check Error if the type in array is not consistent, or if the given type is the highest. Note that only in retype but retype_key, AbstractFloat is also include Integer.
# Missing 2nd arg: Convert to the type of first elements without checking others.
#See also: demissing, promote, typein, tidy, tidytype
#Xiong Jieyi, January 7, 2015

export retype,retype_key
retype(X::Array,T::Type)=convert(Array{T},X)
function retype(X::Array)
    if isempty(X)
        return Union{}[]
    end
    typ=typeof(X[1])
    retype(X,typ)
end
function retype(X::Array,P::Symbol,T::Type=Union{})
    if isempty(X)
        return T[]
    end
    if P==:share
        typ=T
        for cX in X
            typ=promote_type(typ,typeof(cX))
        end
        retype(X,typ)
    elseif P==:check
        if T==Union{}
            typ=typeof(X[1])
            for i=2:length(X)
                @assert(typ==typeof(X[i]),"Array element types are inconsistent. $typ-$(typeof(X[i]))")
            end
            retype(X,typ)
        else
            typ=T<:AbstractFloat ? Union{T,Integer} : T
            for cX in X
                @assert(typeof(cX)<:typ,"Array element types are inconsistent. $T-$(typeof(cX))")
            end
            retype(X,T)
        end
    else
        error("Unknown parameter $P (Only support :share and :check).")
    end
end
# function retype_key{Tk,Tv}(X::Dict{Tk,Tv},T::Type)
#     if Tk==T
#         return X
#     end
#     D=Dict{T,Tv}()
#     for (k,v) in X
#         D[k]=v
#     end
#     D
# end
retype_key(X::Dict{Tk,Tv},T::Type) where {Tk,Tv}=convert(Dict{T,Tv},X)
function retype_key(X::Dict{Tk,Tv}) where {Tk,Tv}
    if isempty(X)
        return Dict{Union{},Tv}()
    end
    T=typeof(first(keys(X)))
    retype_key(X,T)
end
function retype_key(X::Dict{Tk,Tv},P::Symbol,T::Type=None) where {Tk,Tv}
    if isempty(X)
        return Dict{T,Tv}()
    end
    if P==:share
        typ=T
        for ck in keys(X)
            typ=promote_type(typ,typeof(ck))
        end
        retype_key(X,typ)
    elseif P==:check
        if T==Union{}
            typ=typeof(first(keys(X)))
            for ck in keys(X)
                @assert(typ==typeof(ck),"Dict key types are inconsistent. $typ-$(typeof(ck))")
            end
            retype_key(X,typ)
        else
            for ck in keys(X)
                @assert(typeof(ck)<:T,"Dict key types are inconsistent. $T-$(typeof(ck))")
            end
            retype_key(X,T)
        end
    else
        error("Unknown parameter $P (Only support :share and :check).")
    end
end

#deany(X::Array)=convert(Array{typeof(X[1])},X)
#deany_key(Dict)|deany_value(Dict) #Only deany key or value.
#Convert any or union-type array to specific element type according to its first element.
#See also: promote, typein, tidy, tidytype
#Xiong Jieyi,11 Oct 2014>November 23, 2014>January 7, 2015

# function deany(X::Array)
#     warn("Will be obsoleted. Instead by retype.")
#     typ=typeof(X[1])
#     if typ<:AbstractFloat
#         float(X)
#     elseif typ<:Integer
#         int(X)
#     else
#         convert(Array{typ},X)
#     end
# end
# function deany_key(X::Dict)
#     warn("Will be obsoleted. Instead by retype.")
#     keytyp,valtyp=typein_dict(X)
#     keytyp=typeof(first(keys(X)))
#     T=Dict{keytyp,valtyp}()
#     for (k,v) in X
#         T[k]=v
#     end
#     T
# end
# function deany_value(X::Dict)
#     warn("Will be obsoleted. Instead by retype.")
#     keytyp,valtyp=typein_dict(X)
#     valtyp=typeof(first(values(X)))
#     T=Dict{keytyp,valtyp}()
#     for (k,v) in X
#         T[k]=v
#     end
#     T
# end
#export deany, deany_key, deany_value

#}}

#{{ disent
#[[ disent ]]
# A, B, ... = disent( [(a1, b1, ...), (a2, b2, ...); ...] )
# A, B, ... = disent( [[a1, b1, ...], [a2, b2, ...]; ...] )
# A, B, ... = disent( [a1=>b1, a2=>b2, ...]; ...] )
# Disentangle tuple array into seperated arrays.
# The input array elements could be array, tuple or pair with the same lengths. Function only use the length of the first tuple, and the lengths of other tuples will not be checked. In the case any other tuples are longer than the first tuple, function will ignore the remains without any warning or error. Input an empty array is not allowed.
# See also: mxexpand, pmapwd
# Xiong Jieyi, 8 Feb 2019

export disent
function disent(X::Union{AbstractArray})
    t=map(1:length(X[1])) do i
        map((x::Union{AbstractVector, Tuple, Pair})->x[i], X)
    end
    tuple(t...)
end
#}}

#{{ mxexpand
#[[ mxexpand ]]
#Y1, Y2, ... = mxexpand( X1, X2, ... )
#Expand arrays in a similar rule as broadcast().
#See also: disent, coo2mx, eachcombr, mxasvecdo
#Xiong Jieyi, 1 Apr 2019

export mxexpand
mxexpand(X...)=disent(broadcast((x...)->x, X...))
#}}

#{{ demissing
#[[ demissing ]]
# Array = demissing(MissingArray; rp=missing_replacement|justtry=false)
# Convert a missing array to normal array. If the array contains any missing elements, it will throw an error unless 1) rp=... is assigned or 2) justtry=true (return itself without copy).
# See also: withmissing, retype, retype_key, dictfun, dictfun!
# Xiong Jieyi, 2 Mar 2019

export demissing
function demissing(X::AbstractArray{Union{T, Missing}}; rp::Union{T, UndefInitializer}=undef, justtry::Bool=false) where T
    if rp!=undef
        X=copy(X)
        X[ismissing.(X)].=rp
    end
    if justtry && any(ismissing, X)
        X
    else
        convert(AbstractArray{T}, X)
    end
end
#}}

#{{ withmissing
#[[ withmissing ]]
# Array{Union{Missing, ...}} = withmissing(Array)
# Convert input array to the missing-tolerate array. Note that if input array contains missing, the function just output this array itself without copy.
# See also: demissing, grpfunexp
# Xiong Jieyi, 2 Mar 2019

export withmissing
withmissing(X::AbstractArray{T}) where {T>:Missing} = X
withmissing(X::AbstractArray{T}) where {T} = convert(AbstractArray{Union{Missing, T}}, X)
#}}

#{{ pyfun @pyf_str
#[[ pyfun @pyf_str ]]
# ... = pyfun("python_fun"|:python_fun, args...;wargs...)
#  Call buildin python function.
#
# ... = pyf"python_cmd_with_$1_$2..."(...)
# PyObject = pyf"..."o(...)
#  @pyf_str will create a lambda function. Only support one line. Require PyCall being pre-loaded.
#
#See also: @pyhelp
#Xiong Jieyi, 16 Sep 2014>December 8, 2014>4 Mar 2019>26 Jun 2019>29 Nov 2019

export pyfun, @pyf_str
function pyfun(f::AbstractString,args...;wargs...)
    @pkgfun(PyCall, pycall, pybuiltin)
    pycall(pybuiltin(f), Main.PyCall.PyAny, args...; wargs...)
end
pyfun(f::Symbol,args...;wargs...)=pyfun(string(f),args...;wargs...)

macro pyf_str(F::String, M::String="j")
    if occursin(r"(?<!\\)\$(\d+)", F) || occursin(r"(?<!\\)\$\((\d+)\)", F)
        if occursin('\n', F)
            error("Python command in multiple lines is not supported so far.")
        end
        F=replace(F, r"(?<!\\)\$(\d+)" => s"(insidE_valuE_dO_noT_usE[\1-1])")
        F=replace(F, r"(?<!\\)\$\((\d+)\)" => s"(insidE_valuE_dO_noT_usE[\1-1])")
        F="lambda *insidE_valuE_dO_noT_usE : $F"
    end
    if M=="o"
        :((X...)->Main.PyCall.pycall(Main.PyCall.@py_str($F, "o"), Main.PyCall.PyObject, X...))
    else
        @eval(@macroexpand(Main.PyCall.@py_str($F)))
    end
end
#}}

#{{ randint

#[[ randint ]]
#Syntax: randint([Type=Int, ]N[, row_num, col_num, ...])
#Return matrix of integer random number 1~N.
#randint(10,5,3) return a 5 x 3 matrix with the number of 1~10
#Xiong Jieyi, Mar 12, 2014>December 12, 2014 >25 Jun 2018

randint(x::Real,arg...)=ceil.(Int, rand(arg...)*x)
randint(T::Type,x::Real,arg...)=ceil.(T, rand(arg...)*x)
export randint

#}}

#{{ rvorder [ = invperm(), abolished outside ]
#Syntax: rvidx=rvorder(idx)
#When A[idx]==B, B[rvidx]==A
#Xiong Jieyi, Mar 12, 2014>4 Oct 2014

const rvorder = invperm
# function rvorder(x::AbstractVector{T}) where {T<:Real}
#     t=zeros(T,size(x))
#     t[x]=1:length(x)
#     return t
# end
# export rvorder
#}}

#{{ drawer drawergrp
#[[ drawer drawergrp ]]
#Vector = drawer(N, Group_Num)
#  Devide N into given group as even as possible, return number of each group.
#GroupNO = drawergrp(N, Group_Num; shuffle=false)
#  Return a ordered or randomly shuffled group NO vector, composited by 1:Group_Num,  with the length of N.
# For both functions, the higher groups may less one element.
# e.g. drawer(11,3)->[4,4,3]
# See also: segpileup
# Xiong Jieyi, October 19, 2014 > 13 May 2020

export drawer, drawergrp
function drawer(N::Integer, G::Integer)
    D,R=divrem(N,G)
    O=fill(D,G)
    O[1:R]+=1
    O
end
function drawergrp(N::Integer, G::Integer; shuffle::Bool=false)
    Ls=drawer(N, G)
    C=map(enumerate(Ls)) do (i, L)
        fill(i, L)
    end
    O=vcat(C...)
    if shuffle
        O[randperm(length(O))]
    else
        O
    end
end
#}}

#{{ segpileup
#[[ segpileup ]]
# Layer=segpileup(pos)
#Pile segments up with the priority of lower layer. Return layer number. Input is a n x 2 matrix.
#See also: drawer
#Xiong Jieyi, December 9, 2014

#<v0.6# function segpileup{T<:Real}(pos::Matrix{T})
function segpileup(pos::Matrix{T}) where {T<:Real}
    (spos,idx)=sortr(pos)
    pileNO=ones(Int,size(spos,1))
    pileend=[spos[1,2]]
    for i=2:size(spos,1)
        flag=true
        for j=1:length(pileend)
            if spos[i,1]>pileend[j]
                pileNO[i]=j
                pileend[j]=spos[i,2]
                flag=false
                break
            end
        end
        if flag
            pileNO[i]=length(pileend)+1
            push!(pileend,spos[i,2])
        end
    end
    return pileNO[rvorder(idx)]
end
export segpileup
#}}

#{{ filefun filefunwith fileloop fileiofun
#[[ filefun filefunwith fileloop fileiofun ]]
# (Array1,Array2,...)=filefun(function, filename|Cmd|IO; input_rowno|no=false, mrows=false, skipline=0, ignore::Char='\0', ischomp=true, each=0, maxrow=0, nowarn=false)
# rownum=fileloop(function, filename|Cmd|IO; input_rowno|no=false, skipline=0, ignore::Char='\0', ischomp=true, each=0, maxrow=0)
# (infile_rownum, outfile_rownum) = fileiofun(function, in_filename|Cmd|IO, out_filename|cmd|IO="/dev/tty"; input_rowno|no=false, skipline=0, ignore::Char='\0', ischomp=true, each=0, maxrow=0, append=false) #'append' assign the open behave of output file.
# Run function(arg1_rowi,arg2_rowi,...) for each row of file and pile-up the results.
# If mrows is true, outputs of given function will be collected by addrows!() instead of addrow!(), which means multiple row output is allowed.
# If the input_rowno is true, the first argument for given function is row number, or row_num / each when each is given.
# If ischomp is true, the newliner at the end of each line will be chomped.
# The first #skipline line(s) and the later lines start with ignore charactor will be ignored. The ignored lines will not be counted as row_NO.
# If maxrow>0, function will stop after get the max rowno.
# In fileiofun, when the "each" args >0, function will input a N length vector rather than a string. It triggle the function for every N lines. Note that each=0 is different from each=1.
#In both filefun and filefunwith function, if the given function return a nothing, this record will be ignored.
#row loop just run the function but do not record any result.
#In fileiofun, function can only output three types: AbstractString, AbstractString vector( writed as multi-lines) and nothing(ignore writting). A new-line will be added automatically.
#
#Abolishing this function in the future:
# ((Array1,Array2,...), row_NO)=filefunwith(function, filename|Cmd|IO; input_rowno|no=false, multi_output=false, rowexp=false, skipline=0, ignore::Char='\0', ischomp=true, each=0, maxrow=0, nowarn=false)
# Consider using filefun(...; mrows=true) instead.
# (The differences between rowfun and rowfunwith: rowfun require the output of fun should be one row, while in rowfunwith fun can output any number of lines. Note that the vector output in rowfun is regarded as a row, while in rowfunwith, be regarded as a column. rowfunwith is also supported empty output. rowexp is only useful while multi_output is true.)
#The enter at the end of each line will be removed.
#
#See also: getfiles, chomp, dictfun, grpfun, grpfunwith, grpfunexp, rowfun rowfunwith rowloop, readheadmx, read2headmx, flow
#Xiong Jieyi, 5 Sep 2014>1 Oct 2014>February 23, 2015>May 16, 2015>Feb 23, 2016>Apr 16, 2016>28 Mar 2024

mutable struct _filefuncore
    buf::Vector{AbstractString}
    rowno::Int
    eachp::Int
    function _filefuncore()
        new(AbstractString[],0,1)
    end
end
function _runfilefun(obj::_filefuncore, fun::Function, cline::AbstractString, islastline::Bool; no::Bool=false, input_rowno::Bool=no, each::Int=0, ischomp::Bool=true)
    if ischomp
        cline=chomp(cline)
    end
    crlt=nothing
    if each==0
        obj.rowno+=1
        try
            if input_rowno
                crlt=fun(obj.rowno,cline)
            else
                crlt=fun(cline)
            end
        catch err
            println("Error occured in file line $(obj.rowno):")
            sk(cline)
            rethrow(err)
        end
        (crlt, obj.rowno, true)
    else
        push!(obj.buf,cline)
        if obj.eachp<each
            obj.eachp+=1
            (nothing, obj.rowno, false)
        else
            obj.rowno+=1
            try
                if input_rowno
                    crlt=fun(obj.rowno,obj.buf)
                else
                    crlt=fun(obj.buf)
                end
            catch err
                println("Error occured in file line-group $(obj.rowno):")
                sk(obj.buf)
                rethrow(err)
            end
            obj.buf=AbstractString[]
            obj.eachp=1
            (crlt, obj.rowno, true)
        end
    end
end

function filefun(fun::Function,fid::IO; multi_output::Bool=false, skipline::Integer=0, ignore::Char='\0', maxrow::Integer=0, nowarn::Bool=false, mrows::Bool=false, wargs...)
    for i=1:skipline
        readline(fid)
    end
    rlt=rowpile()
    rowno=0
    stt=_filefuncore()
    for cline in eachline(fid)
        if ignore!='\0' && cline[1]==ignore
            continue
        end
        crlt,rowno,t=_runfilefun(stt,fun,cline,eof(fid);wargs...)
        if !t
            continue
        end
        try
            if !isnothing(crlt) #!isa(crlt,Void)
                if multi_output
                    # This case is abolished. Kept here only for compatible reason.
                    @assert(!mrows)
                    addrow!(rlt,crlt...)
                else
                    if mrows
                        addrows!(rlt,crlt)
                    else
                        addrow!(rlt,crlt)
                    end
                end
            end
        catch err
            println("Error occured in file line $rowno")
            sk(cline)
            println("With function returns:")
            sk(crlt)
            rethrow(err)
        end
        if maxrow>0 && rowno>=maxrow
            break
        end
    end
    if rownum(rlt)==0
        nowarn || @warn("No output collected, return nothing instead.")
        return nothing
    else
        return value(rlt, keep_tuple=multi_output)
    end
end
filefun(fun::Function,fn;args...)=open(x->filefun(fun,x;args...),fn)

function filefunwith(fun::Function,fid::IO;multi_output::Bool=false,skipline::Integer=0,ignore::Char='\0', maxrow::Integer=0, nowarn::Bool=false, wargs...)
    #This function has been abolished since 28 Mar 2024. Kept here only for compatible reason.
    for i=1:skipline
        readline(fid)
    end
    rlt=rowpile()
    rowno=0
    stt=_filefuncore()
    for cline in eachline(fid)
        if ignore!='\0' && cline[1]==ignore
            continue
        end
        crlt,rowno,t=_runfilefun(stt,fun,cline,eof(fid);wargs...)
        if !t
            continue
        end
        try
            if !isa(crlt,Void)
                rprowno=fill(rowno,rownum(crlt))
                if multi_output
                    addrows!(rlt,rprowno,crlt...)
                else
                    addrows!(rlt,rprowno,crlt)
                end
            end
        catch err
            println("Error occured in file line $rowno")
            sk(cline)
            println("With function returns:")
            sk(crlt)
            rethrow(err)
        end
        if maxrow>0 && rowno>=maxrow
            break
        end
    end
    if isvirgin(rlt)
        nowarn || @warn("No any output collected, return nothing instead.")
        return (nothing,Int[])
    else
        t=value(rlt)
        if multi_output
            return (t[2:end],t[1])
        else
            return (t[2],t[1])
        end
    end
end
filefunwith(fun::Function,fn;args...)=open(x->filefunwith(fun,x;args...),fn)

function fileloop(fun::Function,fid::IO;skipline::Integer=0,ignore::Char='\0',maxrow::Integer=0,wargs...)
    for i=1:skipline
        readline(fid)
    end
    rowno=0; trownum=0
    stt=_filefuncore()
    for cline in eachline(fid)
        if ignore!='\0' && cline[1]==ignore
            continue
        end
        _,rowno,t=_runfilefun(stt,fun,cline,eof(fid);wargs...)
        if t && maxrow>0 && rowno>=maxrow
            break
        end
    end
    return rowno
end
fileloop(fun::Function,fn;args...)=open(x->fileloop(fun,x;args...),fn)

function fileiofun(fun::Function,infn::IO,outfn::IO="/dev/tty";
                   skipline::Integer=0,ignore::Char='\0',maxrow::Integer=0,wargs...)
    for i=1:skipline
        readline(infn)
    end
    inrowno=0
    outrowno=0
    stt=_filefuncore()
    for cline in eachline(infn)
        if ignore!='\0' && cline[1]==ignore
            continue
        end
        crlt,inrowno,t=_runfilefun(stt,fun,cline,eof(infn);wargs...)
        if !t
            continue
        end
        if isa(crlt,AbstractString)
            # write(outfn,crlt*"\n")
            println(outfn,crlt)
            outrowno+=1
        elseif isa(crlt,Vector) && eltype(crlt)<:AbstractString && !isempty(crlt)
            # write(outfn,join(map(x->x*"\n",crlt)))
            foreach(x->println(outfn,x),crlt)
            outrowno+=length(crlt)
        else
            @assert(isa(crlt,Void) || isempty(crlt),
                    "Invalid function output type $(typeof(crlt)). Only AbstractString, Vector of string and nothing are supported.")
        end
        if maxrow>0 && inrowno>=maxrow
            break
        end
    end
    (inrowno, outrowno)
end
fileiofun(fun::Function,infn,args...;wargs...)=open(x->fileiofun(fun,x,args...;wargs...), infn)
fileiofun(fun::Function,infn::IO,outfn,args...;append::Bool=false,wargs...)=open(x->fileiofun(fun,infn,x,args...;wargs...),outfn,append ? "a" : "w")
export filefun, filefunwith, fileloop, fileiofun

#}}

#{{ readheadmx read2headmx
#[[ readheadmx read2headmx ]]
# (Matrix, RowHead)=readheadmx(filename, fun::Function_or_Type=float; dm='\t', skipline=0, ignore='')
# (Matrix, RowHead, ColumnHead)=read2headmx(...; nowarn=false, ...)
# ... = readheadmx|read2headmx(Function, filename; ...) #is also supported.
# Read matrix tsv data with row head (readheadmx) or with both head (read2headmx).
#  * The first line after skipline and ignore-lines will be regarded as column head.
# For read2headmx:
#  * Both RowHead and ColumnHead are vectors.
#  * The top-left cell will be ignored if its length just be matrix width+1. If this cell is not empty, a warning will show unless nowarn=true.
#See also: filefun, readtb.
#Xiong Jieyi, Nov 27, 2015 > 6 Feb 2019 > 30 Jul 2019

function readheadmx(filename::AbstractString,fun::Function=float;
                    dm::Union{Char,AbstractString}='\t',skipline::Integer=0,ignore::Char='\0',headfun::Function=ascii)
    filefun(filename;multi_output=true,skipline=skipline,ignore=ignore) do li
        C=split(li,dm)
        (map(fun,C[2:end]), headfun(C[1]))
    end
end
function read2headmx(filename::AbstractString,fun::Function=float;
                     dm::Union{Char,AbstractString}="\t",skipline::Integer=0,ignore::Char='\0',headfun::Function=ascii,colheadfun::Function=ascii, nowarn::Bool=false)
    flag=true
    colhead=nothing
    rlt=filefun(filename;multi_output=true,skipline=skipline,ignore=ignore) do li
        C=split(li,dm)
        if flag
            flag=false
            colhead=C
            nothing
        else
            (map(fun,C[2:end]), headfun(C[1]))
        end
    end
    if length(colhead)==size(rlt[1], 2)+1
        if !nowarn && !isempty(colhead[1])
            @warn(f"Droped non-empty top-left corner \"$1\"."(colhead[1]))
        end
        deleteat!(colhead, 1)
    end
    (rlt..., map(colheadfun, colhead))
end
readheadmx(filename::AbstractString, typ::Type; wargs...)=readheadmx(filename, x->parse(typ, x); wargs...)
read2headmx(filename::AbstractString, typ::Type; wargs...)=read2headmx(filename, x->parse(typ, x); wargs...)
readheadmx(fun::Function, filename::AbstractString; wargs...)=readheadmx(filename, fun; wargs...)
read2headmx(fun::Function, filename::AbstractString; wargs...)=read2headmx(filename, fun; wargs...)
export readheadmx, read2headmx
#}}

#{{ writeall
#[[ writeall ]]
#writeall(filename|CMD, data; append=false)
#Write all the data to file. If data is a string, number, symbol or bool, it will be "print"-ed to file directly. If data is vector of string, number, symbol or bool, each elements will be "println"-ed as one line. Not support other types or matrix.
#See also: writedlm
#Xiong Jieyi, 9 Oct 2014 >December 23, 2014>9 Jan 2019

writeall(fn, Data::Union{AbstractString, Symbol, Number, Bool}; append::Bool=false)=open(x->print(x, Data), fn, append ? "a" : "w")
# writeall(fn,Data::Real; append::Bool=false)=writeall(fn,num2str(Data);append=append)
# writeall(fn,Data::Vector{T};append::Bool=false) where {T<:AbstractString} =open(x->write(x,join(Data,"\n")),fn,append ? "a" : "w")
# writeall(fn,Data::AbstractVector{T};append::Bool=false) where {T<:Real} =writeall(fn,num2str(Data);append=append)
function writeall(fn, Data::AbstractVector{T}; append::Bool=false) where {T<:Union{AbstractString, Number, Symbol, Bool}}
    open(fn, append ? "a" : "w") do fid
        for x in Data
            println(fid, x)
        end
    end
end
# writeall(fn,Data::AbstractArray; append::Bool=false)=error("writeall only support scalar or vector. For matrix, try writedlm instead.")
export writeall
#}}

#{{ getfiles
#[[ getfiles readdir_recur ]]
#(file_name_vec, out2_vec)=getfiles(path=".", nothing|r"Regex"|"ContainingStr"|Function|(cond1, cond2, ...); fullpath=true, recur=false, max_depth=Inf, typ=:all(default)|:file|:dir, out2=:name|:nameext|:folder|:path, dirfailed=:error|:warn|:ignore)
#Or ... = getfiles(Function, path="."; ...)
#Return a vector of full file-filenames match the given filter, as well as a vector of pure_name (no path and no ext).
#When 2nd input is a tuple for multiple conditions, funciton only output the files/dirs match all the conditions.
#When recur is true, function will search subdirectory recurly.
#When typ is :file, these is no directory in output. When typ is :dir, there is no file in output. When typ is :all, the directories in the output will ended with '/', and pure_name_vec will be empty "" for the directories items.
# typ could be function(fullname::String)::Bool, used to filter the output files.
#out2 designs what is the second output, the file name without ext (:name, default), the file name with ext (:nameext, default when typ=:dir), the last folder name (:folder) or the relative path (:path).
#
#file_vec = readdir_recur(path="."; max_depth::Real=Inf, flt=nothing, typ=:all)
#Similar as readdir, but do it recurly. flt is the filter input similar to the 2nd input of getfiles(), so does the typ.
#dirfailed= assigns the behave when reading a dir with system error, such as Permission Denied. It is not effect for the most-top directory.
#
#See also: matchfile readdir, joinpath, splitdir, findinfile, filefun, filefunwith, fileloop
#Xiong Jieyi, 16 Sep 2014>January 12, 2015>Jul 22, 2015>Apr 30, 2016>2 Jan 2019>29 Jan 2019

export getfiles, readdir_recur
function __runUntilSystemError__(fun::Function, X...; dirfailed::Symbol=:error, default=missing)
    try
        fun(X...)
    catch err
        if isa(err, SystemError) || isa(err, Base.IOError)
            if dirfailed==:warn
                @warn err
            elseif dirfailed!=:ignore
                rethrow(err)
            end
        else
            rethrow(err)
        end
        default
    end
end
function readdir_recur(path::AbstractString="."; max_depth::Real=Inf, flt::Union{Regex,AbstractString,Function,Void,Tuple}=nothing, typ::Symbol=:all, dirfailed::Symbol=:error)
    out=AbstractString[]
    function _recur(path, depth)
        fns=__runUntilSystemError__(readdir, path, dirfailed=dirfailed, default=String[])
        if path!="."
            fns=map(x->joinpath(path,x),fns)
        end
        for fn in fns
            ofn, flag=__check_file__(fn, flt, typ, dirfailed=dirfailed)
            if flag
                push!(out, ofn)
            end
            if __runUntilSystemError__(isdir, fn, dirfailed=dirfailed, default=false) && depth < max_depth
                _recur(fn, depth+1)
            end
        end
    end
    _recur(path, 1)
    out
end
function __check_file__(fn::AbstractString, flt::Union{Regex,AbstractString,Function,Void,Tuple}, typ::Symbol; dirfailed::Symbol=:error)
    E=if typ==:file
        # isfile(fn)
        __runUntilSystemError__(isfile, fn, dirfailed=dirfailed, default=false)
    elseif typ==:dir
        # isdir(fn)
        __runUntilSystemError__(isdir, fn, dirfailed=dirfailed, default=false)
    else
        typ==:all || error( "typ can only be :all(default)|:file|:dir.")
        t=__runUntilSystemError__(isdir, fn, dirfailed=dirfailed, default=missing)
        if ismissing(t)
            false
        else
            if t
                fn=fn*"/"
            end
            true
        end
    end
    if E && flt!=nothing
        E=if isa(flt, Tuple)
            flag=true
            for cflt in flt
                flag=if isa(cflt, Function)
                    cflt(fn)::Bool
                else
                    occursin(cflt, fn)
                end
                if !flag
                    break
                end
            end
            flag
        else
            if isa(flt, Function)
                flt(fn)::Bool
            else
                occursin(flt, fn)
            end
        end
    end
    (fn, E)
end
function getfiles(path::AbstractString=".", flt::Union{Regex,AbstractString,Function,Void,Tuple}=nothing; fullpath::Bool=true, recur::Bool=false, max_depth::Real=Inf, typ::Symbol=:all, out2::Symbol=typ==:dir ? :nameext : :name, dirfailed::Symbol=:error)
    isCurPath=path=="."
    if !isCurPath && !isdir(path)
        error("Path $path does not exists.") #Making sure the most top dir exists.
    end
    if recur
        fns=readdir_recur(path, max_depth=max_depth, flt=flt, typ=typ, dirfailed=dirfailed)
    else
        fns=AbstractString[]
        for fn in readdir(path)
            ofn, flag=__check_file__(isCurPath ? fn : joinpath(path, fn), flt, typ, dirfailed=dirfailed)
            if flag
                push!(fns, ofn)
            end
        end
    end

    if out2==:name
        purename=map(x->splitext(splitdir(x)[2])[1],fns)
    elseif out2==:nameext
        purename=map(x->splitdir(x)[2],fns)
    elseif out2==:folder
        purename=map(x->splitdir(splitdir(x)[1])[2],fns)
    elseif out2==:path
        purename=map(x->splitdir(x)[1],fns)
    else
        error("Invalid out2. out2 can only be :name, :name.ext, :folder, :path.")
    end
    
    if fullpath
        fns=map(abspath,fns)
    end
    (fns, purename)
end
getfiles(fun::Function, path::AbstractString="."; warg...)=getfiles(path, fun; warg...)
#}}

#{{ matchfile
#[[ matchfile ]]
# filename|"" = matchfile("filename_with_?|*"; onlyone=true(D), hasone=onlyone)
# filename_vec|String[] = matchfile("filename_with_?|*"; onlyone=false, hasone=onlyone)
#Match filename and return the completed one(s).
#In default (onlyone=true), function always return a string, and blank string for mismatch. When onlyone=false, function always return a string vector, and empty vector for mismatch.
#When hasone is true, funciton will throw an error when mismatch.
#See also: getfiles, isfile, isdir
#Xiong Jieyi, Jul 8, 2016 > 14 Jan 2019

export matchfile
function matchfile(nm; onlyone::Bool=true, hasone::Bool=onlyone)
    nm=strip(nm)
    if in('?',nm) || in('*',nm)
        S=rstrip(readchomp(sfc"printf '%s\0' $1"(nm)), '\0')
        if S==nm
            if hasone
                error("No file can match $nm")
            elseif onlyone
                ""
            else
                String[]
            end
        elseif onlyone
            if in('\0', S)
                error("More than one file has matched $nm.")
            else
                S
            end
        else
            string.(split(S, '\0'))
        end
    else
        error("No ? or * in the input #1.")
        # if isfile(nm) || isdir(nm)
        #     nm
        # elseif hasone
        #     error("No file can match $nm")
        # else
        #     ""
        # end
    end
end
#}}

#{{ capstdout
#[[ capstdout ]]
# out_string, fun_output = capstdout( f::Function )
#Capture function stand output to a string. The captured outputs will not be showed on screen.
#See also: KK, writeall, runstr, successstr, @repress_stdout, @repress_stderr
#Xiong Jieyi, Jun 17, 2015 >Jan 15, 2016>4 Mar 2022

export capstdout
function capstdout(fun::Function)
    oldout=stdout
    fout=nothing
    out=Base.PipeEndpoint()
    try
        out,=redirect_stdout()
        fout=fun()
    catch err
        rethrow(err)
    finally
        redirect_stdout(oldout)
    end
    String(readavailable(out)), fout
end
#}}

#{{ my.isobj my.getobj
#[[ my.isobj my.getobj ]]
# T|F = my.isobj( "A.[:B].C..." )
# Obj|undef = my.getobj( "A.[:B].C..." )
#Return if the named object is in Main.
#This function is not exported since 13 Jun 2019.
#See also: 
#Xiong Jieyi, March 9, 2015 > 13 Jun 2019

#=
function isobj(S::AbstractString)
    md=Main
    C=split(S, '.')
    N=length(C)
    for (i, t) in enumerate(C)
        ss=Symbol(t)
        if isa(md, Module)
            if !isdefined(md, ss)
                return false
            end
        elseif isdefined(Main, :PyObject) && isa(md, Main.PyObject)
            if !pyfun("hasattr", md, ss)
                return false
            end
        else
            if !(ss in propertynames(md))
                return false
            end
        end
        # md=Base.MainInclude.eval(Meta.parse(S[1:i-1]))
        if i<N
            md=getproperty(md, ss)
        end
    end
    return true
end
function getobj(S::AbstractString)
    if isempty(S)
        return nothing
    end
    md=Main
    C=split(S, '.')
    N=length(C)
    for (i, t) in enumerate(C)
        ss=Symbol(t)
        if isa(md, Module)
            if !isdefined(md, ss)
                return nothing
            end
        elseif isdefined(Main, :PyObject) && isa(md, Main.PyObject)
            if !pyfun("hasattr", md, ss)
                return false
            end
        else
            if !(ss in propertynames(md))
                return nothing
            end
        end
        md=getproperty(md, ss)
    end
    return md
end
=#
isobj(S::AbstractString)=!isa(getobj(S), UndefInitializer)
function getobj(S0::AbstractString)
    function _getobj(S::Symbol, get_variable::Bool=false)
        if get_variable
            if isdefined(Main, S)
                getproperty(Main, S)
            else
                undef
            end
        else
            S
        end
    end
    _getobj(S::QuoteNode, _::Bool=false)=S.value
    _getobj(S, _::Bool=false)=S
    function _getobj(S::Expr, _::Bool=false)
        md=_getobj(S.args[1], true)
        if isa(md, UndefInitializer)
            return undef
        end
        if S.head==:.
            if length(S.args)!=2
                return undef
            end
            ss=_getobj(S.args[2])
            if !(isa(ss, Symbol) || isa(ss, AbstractString))
                return undef
            end
            if isa(md, Module)
                if isdefined(md, ss)
                    return getproperty(md, ss)
                else
                    return undef
                end
            elseif isdefined(Main, :PyObject) && isa(md, Main.PyObject)
                if Main.PyCall.pycall(Main.PyCall.pybuiltin("hasattr"), Bool, md, ss)
                    return getindex(md, ss)
                else
                    return undef
                end
            else
                if ss in propertynames(md)
                    return getproperty(md, ss)
                else
                    return undef
                end
            end
        elseif S.head==:ref
            if length(S.args)<2
                return undef
            end
            ss=_getobj.(S.args[2:end])
            if isa(md, Union{AbstractArray, Tuple})
                if all(x->isa(x, Int64), ss)
                    try
                        return getindex(md, ss...)
                    catch err
                        if isa(err, BoundsError)
                            return undef
                        else
                            rethrow(err)
                        end
                    end
                else
                    return undef
                end
            elseif length(ss)==1 && (isa(ss[1], AbstractString) || isa(ss[1], Symbol))
                ss1=ss[1]
                if isdefined(Main, :PyObject) && isa(md, Main.PyObject)
                    if Main.PyCall.pycall(Main.PyCall.pybuiltin("hasattr"), Bool, md, ss1)
                        getindex(md, ss1)
                        # if isa(ss1, AbstractString)
                        #     return Main.PyCall.pycall(Main.PyCall.pybuiltin("getattr"), Main.PyObject, ss1)
                        # else
                        #     return getindex(md, ss1)
                        # end
                    else
                        return undef
                    end
                else
                    return get(md, ss1, undef)
                end
            else
                return undef
            end
        else
            error("Invalid object format: $S0")
        end
    end
    _getobj(Meta.parse(S0), true)
end
#}}

#{{ tryfun
#[[ tryfun ]]
# tryfun(function, times=2)
#If function has a error occur, try it again.
#See also: 
#Xiong Jieyi, Nov 18, 2015

export tryfun
function tryfun(fun::Function,N::Integer=2)
    for i=1:N-1
        try
            fun()
        catch err
            @warn("[ Tryfun ERROR ] "*err.msg)
            println("Try again $(i+1) / $N ...")
        end
    end
    fun()
end
#}}

#{{ isijulia
#[[ isijulia ]]
#T|F=isijulia()
#Showing if in the IJulia env.
#See also: importpy(:IPython!display)
#Xiong Jieyi, December 12, 2014

export isijulia
isijulia()=isdefined(Main, :IJulia) && Main.IJulia.inited
#}}

#{{ getmethods
#[[ getmethods ]]
#methods_vec = getmethods(Type|X)
#methods_vec = getmethods(Type|X, Function|"fun_name")
#Return a string vector of methods of given type, or the function args with given types. If the first input is not a Type, it equal to input typeof(X).
#When the 2nd input is string rather than function, all sub-types will be considered.
#See also: addhelp
#Xiong Jieyi, March 10, 2015

export getmethods
function getmethods(T::Type)
    ms=methodswith(T,true)
    fs=map(ms)do x
        s=match(r"^[^\(\{]+",string(x))
        isa(s,Void) ? string(x) : s.match
    end
    sort(unique(fs))
end
function getmethods(T::Type, S::AbstractString)
    ms=map(string,methodswith(T,true))
    filter!(x->startswith(x,S),ms)
end
function getmethods(T::Type,F::Function)
    ms=ASCIIString[ string(x) for x in methods(F)]
    filter!(Regex("\\:$T[\\,\\)\\}]"),ms)
end
getmethods(::T,args...) where {T}=getmethods(T,args...)
#}}

#{{ trsp
#[[ trsp ]]
#M' = trsp(M)
#Transpose matrix but not do it recurly.
#BTW, I also write Base.transpose(::AbstractString and Symbol) so "'" can also transpose string array.
#See also: h2v, mxasvecdo
#Xiong Jieyi, Jun 5, 2015 >Sep 7, 2015 >Nov 11, 2016>Feb 12, 2017

export trsp, transpose
import Base.transpose
transpose(x::AbstractString)=x
transpose(x::Symbol)=x

# trsp(X)=X'
trsp(X::AbstractMatrix)=permutedims(X, [2,1])
trsp(X::AbstractVector)=reshape(X,1,length(X))
function trsp(X::AbstractMatrix{Any})
    rN,cN=size(X)
    O=Array{Any}(undef, cN, rN)
    for ri=1:rN, ci=1:cN
        O[ci,ri]=X[ri,ci]
    end
    O
end
# trsp(X::AbstractVector{Any})=trsp(hcat(X))
#<v0.6# function trsp{T<:AbstractArray}(X::AbstractMatrix{T})
function trsp(X::AbstractMatrix{T}) where {T<:AbstractArray}
    rN,cN=size(X)
    O=fill(eltype(T)[],cN,rN)
    for ri=1:rN, ci=1:cN
        O[ci,ri]=X[ri,ci]
    end
    O
end
# trsp{T<:AbstractArray}(X::AbstractVector{T})=trsp(hcat(X))
#}}

#{{ h2v
#[[ h2v ]]
# vec = h2v( one-row-matrix )
# Convert a one-row-matrix to a vector. An error will occur if input contain more than one row.
# See also: trsp, col2vec, mxasvecdo
# 24 Mar 2022.

export h2v
function h2v(M::AbstractMatrix)
    if size(M, 1)!=1
        error("Input is not an one-row-matrix.")
    end
    vec(M)
end
#}}

#{{ totalsizeof
#[[ totalsizeof ]]
# sizebits=totalsizeof(X)
#Return the size of memory occupied by the X.
#Code is downloaded from https://gist.github.com/avitale/0070b629b89350b39c21#file-totalsize-jl
#and modified by XJY for compat to Julia 1.0.
#See also: saynum
#Xiong Jieyi, Jun 13, 2015 > 3 Feb 2019

# include("totalsize.jl")
# totalsizeof=Totalsizeof.totalsizeof
export totalsizeof

# Helper: Pointer cache is used to break circular references in objects and avoid double countings
function is_seen!(x, ptr_cache)
    if isimmutable(x) #Add in 3 Feb 2019 by XJY to avoid error in Julia 1.0
        false
    else
        in(pointer_from_objref(x), ptr_cache) || (push!(ptr_cache, pointer_from_objref(x)); false)
    end
end

# Helper: Catch types without size method
sizeof_catch(x) = try sizeof(x) catch; 0 end

# Composite Types
function totalsizeof(x, ptr_cache = Set())
  is_seen!(x, ptr_cache) && return 0
  result = sizeof_catch(x)

  this_names = nothing

  # Check if type has names method
  try
    this_names = fieldnames(x)
  catch
    this_names = []
  end

  el = nothing
  isempty(this_names) ? result :
    result + mapreduce(i -> 
      # Check if some elements are not reachable
      (try el = getfield(x,i) catch; el = nothing end; totalsizeof(el, ptr_cache))
      , +, this_names)
end

# Tuples
function totalsizeof(x::Tuple, ptr_cache = Set())
  is_seen!(x, ptr_cache) && return 0
  isempty(x) ? 0 : mapreduce(i -> totalsizeof(i, ptr_cache), +, x)
end

# Array
#<v0.6# function  totalsizeof{T<:Any}(x::Array{T}, ptr_cache = Set())
function  totalsizeof(x::Array{T}, ptr_cache = Set()) where {T<:Any}
    is_seen!(x, ptr_cache) && return 0
    result = sizeof(x)

    #Add by Xiong Jieyi, Jun 13, 2015
    if T<:Number || T<:AbstractString
        return result
    end
    #Add end
    
    el = nothing
    for i in 1:length(x)
        # Check if some elements are undefined
        try el = x[i] catch; el = nothing end
        result += totalsizeof(el, ptr_cache)
    end
    result
end

# Dict #Added by Jieyi Xiong for Julia 1.0, 3 Dec 2018
function totalsizeof(x::Dict, ptr_cache = Set())
  is_seen!(x, ptr_cache) && return 0
    isempty(x) ? 0 : mapreduce(k->totalsizeof(k, ptr_cache), +, keys(x))+mapreduce(v->totalsizeof(v, ptr_cache), +, values(x))
end

# AbstractString, Function, IntrinsicFunction
totalsizeof(x::Union{AbstractString, Function}, ptr_cache = Set()) = is_seen!(x, ptr_cache) ? 0 : sizeof(x)
#}}

#{{ parsestr ParseStr
#[[ parsestr parsestr1 ParseStr next next1 ]]
# Any[value1, value2, ...] = parsestr(txt, commands...; start=1|start_ref=[start])
# value1|vec_of_value1s = parsestr1(txt|txt_vec, commands...; default=...) #Only capture one output. Input can be a vector of string. When unmatched, default value will be returned if assigned.
# obj=ParseStr(txt; start=1)
# Any[value1, value2, ...] = next(obj, commands...)
# value1 = next1(obj, commands...; default=...)
#
# Find txt from the start position. The commands is found by order. Supported commands are: AbstractString, Regex, Int. Int means moving cursor or fetch next/pervious worlds.
# (:search, "substring"|r"pattern"|r"(token)"|offset) --Search pattern or move cursor. Error if mismatch.
# (:test|:try|:no, "substring"|r"pattern"|r"(token)") --Similar as :search. :try will ignore the mismatch error. :test will return true|false. :no will error if matched.
# (:rsearch(=:r)|:rtry, "substring"|offset) --Search from right side. Match should not left to the last match.
# (:get|, "substring"|r"pattern"|r"(token)"|offset[, function[, defval]]) --Fetch contents. In substring model, return the substring from last match to this match. If mismatch, an error will be thrown unless the "defval" was set. defval is not supported in the "offset" mode.
# (:getto(=:gt)|:getuntil(=:gu), "substring"|r"pattern"|r"(token)"[, function[, defval]]) --Fetch the substring from last match to the beginning(:getuntil) or end(:getto) this match. Note (:getuntil, "substring") is invalid, replace it to (:get, "substring").
# (:rget(=:rg), "substring"[, function[, defval]]) --Reversely search pattern and fetch substring between this match and last match.
# (:before|:back, commmand) Run one command only in the region between last match and the match before last match. :back will also move start point back but :before will not. Only one command is supported so far.
# start_ref can be used to both input and return the next start position before and after processed all commands.
#Function always return a Any-vector.
#See also: readall
#Xiong Jieyi, Sep 22, 2015 >Oct 7, 2015>Nov 17, 2015>Feb 5, 2016>Apr 15, 2016>May 12, 2016>May 20, 2016>Jan 18, 2017>13 Jun 2017>9 Jan 2019>11 May 2022>20 Aug 2024

export parsestr, parsestr1
macro __ferror__tmp(msg) #In v1.10, macro cannot be defined inside let-block.
    quote
        if __repress_error
            return nothing
        else
            error($msg)
        end
    end |>esc
end
let
    global parsestr
    # __repress_error=false
    function parsestr(txt::AbstractString, cmds...; start_ref::Vector{Int}=[1],start::Int=start_ref[1],prev_ref::Vector{Int}=[1],lastMatchRange_ref::Vector{Trg}=[0:-1], no_err::Bool=false) where {Trg<:Range}
        #19 Aug 2024: no_err is a hiden argument so far. When no_err=true, function return nothing if any error throwed.
        # prev_ref is a one-element vector of the left boundary of :rsearch and :rget.
        global __repress_error=no_err
        curp=start
        prvp=prev_ref[1]
        mrg::Range=lastMatchRange_ref[1]
        rlts=Any[]
        for cmd in cmds
            t=__parsestr_core(txt,curp,prvp,mrg,cmd)
            if isnothing(t) #no_err=true && an error was throwed.
                return nothing
            end
            rlt,curp,prvp,mrg=t
            rlt!=nothing && push!(rlts,rlt)
        end
        start_ref[1]=curp
        prev_ref[1]=prvp
        lastMatchRange_ref[1]=mrg
        rlts
    end
    function __parsestr_core(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,cmd)
        #curp: current point
        #prvp: pervious point
        #mrg: last matched range
        fun=Dict{Symbol,Function}(:test=>__parsestr_core_test,
                                  :check=>__parsestr_core_test, #Just for compatible.
                                  :rtest=>__parsestr_core_rtest,
                                  :try=>__parsestr_core_try,
                                  :rtry=>__parsestr_core_rtry,
                                  :search=>__parsestr_core_search,
                                  :rsearch=>__parsestr_core_rsearch,
                                  :r=>__parsestr_core_rsearch,
                                  :get=>__parsestr_core_fetch,
                                  :getto=>__parsestr_core_fetchto,
                                  :gt=>__parsestr_core_fetchto,
                                  :getuntil=>__parsestr_core_fetchuntil,
                                  :gu=>__parsestr_core_fetchuntil,
                                  :rg=>__parsestr_core_rfetch,
                                  :rget=>__parsestr_core_rfetch,
                                  :no=>__parsestr_core_no,
                                  :end=>__parsestr_core_end,
                                  :back=>__parsestr_core_back,
                                  :before=>__parsestr_core_before)
        if isa(cmd,Tuple)
            if isa(cmd[1],Symbol)
                fun[cmd[1]](txt,curp,prvp,mrg,cmd[2:end]...)
            else
                __parsestr_core_fetch(txt,curp,prvp,mrg,cmd...)
            end
        elseif cmd == Symbol("end")
            __parsestr_core_end(txt,curp,prvp,mrg)
        else
            __parsestr_core_search(txt,curp,prvp,mrg,cmd)
        end
    end
    
    function __parsestr_core_test(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString)
        R=search(txt,P,curp)
        if isempty(R)
            (false,curp,prvp,R)
        else
            (true,last(R)+1,curp,R)
        end
    end
    function __parsestr_core_rtest(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString)
        R=rsearch(txt[prvp:curp-1],P)
        if isempty(R)
            (false,curp,prvp,R)
        else
            (true,last(R)+1,curp,R)
        end
    end
    function __parsestr_core_test(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Regex)
        R=match(P,txt,curp)
        if R==nothing
            (false,curp,prvp,0:-1)
        elseif isempty(R.captures)
            (true,R.offset+length(R.match),curp,R.offset:length(R.match)-1)
        else
            (true,onlyone(R.offsets)+length(onlyone(R.captures)),curp,R.offsets[1]:R.offsets[1]+length(R.captures[1])-1)
        end
    end
    function __parsestr_core_no(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P)
        rlt,=__parsestr_core_test(txt,curp,prvp,mrg,P)
        rlt && @__ferror__tmp("Pattern $P is found at position $curp.")
        (nothing,curp,prvp,mrg)
    end
    function __parsestr_core_try(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P)
        rlt,ncurp,nprvp,nmrg=__parsestr_core_test(txt,curp,prvp,mrg,P)
        if rlt
            (nothing,ncurp,nprvp,nmrg)
        else
            (nothing,curp,prvp,mrg)
        end
    end
    function __parsestr_core_rtry(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString)
        rlt,ncurp,nprvp,nmrg=__parsestr_core_rtest(txt,curp,prvp,mrg,P)
        if rlt
            (nothing,ncurp,nprvp)
        else
            (nothing,curp,prvp)
        end
    end
    function __parsestr_core_search(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString)
        R=search(txt,P,curp)
        isempty(R) && @__ferror__tmp("Cannot find $P at position $curp.")
        (nothing,last(R)+1,curp,R)
    end
    function __parsestr_core_search(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Regex)
        R=match(P,txt,curp)
        if R==nothing
            @__ferror__tmp("Cannot match $P at position $curp.")
        elseif isempty(R.captures)
            (nothing,R.offset+length(R.match),curp,R.offset:length(R.match)-1)
        else
            (nothing,onlyone(R.offsets)+length(onlyone(R.captures)),curp,R.offsets[1]:R.offsets[1]+length(R.captures[1])-1)
        end
    end
    function __parsestr_core_search(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Int)
        (nothing,curp+P,prvp,curp+P:curp+P)
    end
    function __parsestr_core_rsearch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString)
        R=rsearch(txt,P,curp-1)
        (isempty(R) || first(R)<prvp)  && @__ferror__tmp("Cannot find $P between position $prvp to $curp.")
        (nothing,last(R)+1,curp,R)
    end
    function __parsestr_core_rsearch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Int)
        (nothing,curp-P,prvp,curp-P:curp-P)
    end
    function __parsestr_core_fetch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString,fun::Function=x->x,defval=undef)
        R=search(txt,P,curp)
        # isempty(R) && error("Cannot find $P at position $curp.")
        if isempty(R)
            if defval==undef
                @__ferror__tmp("Cannot find $P at position $curp.")
            else
                return (defval, curp, prvp, 0:-1)
            end
        end
        (fun(txt[curp:first(R)-1]),last(R)+1,curp,R)
    end
    __parsestr_core_fetchuntil(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString,fun::Function=x->x)=__parsestr_core_fetch(txt,curp,prvp,P)
    function __parsestr_core_fetchto(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString,fun::Function=x->x,defval=undef)
        R=search(txt,P,curp)
        if isempty(R)
            if defval==undef
                @__ferror__tmp("Cannot find $P at position $curp.")
            else
                return (defval, curp, prvp, 0:-1)
            end
        end
        (fun(txt[curp:last(R)]),last(R)+1,curp,R)
    end
    function __parsestr_core_fetch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Regex,fun::Function=x->x,defval=undef)
        R=match(P,txt,curp)
        if R==nothing
            if defval==undef
                @__ferror__tmp("Cannot match $P at position $curp.")
            else
                return (defval, curp, prvp, 0:-1)
            end
        elseif isempty(R.captures)
            (fun(R.match),R.offset+length(R.match),curp,R.offset:length(R.match)-1)
        else
            (fun(onlyone(R.captures)),R.offsets[1]+length(R.captures[1]),curp,R.offsets[1]:R.offsets[1]+length(R.captures[1])-1)
        end
    end
    function __parsestr_core_fetch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Int,fun::Function=x->x)
        (fun(txt[curp:curp+P-1]),curp+P,curp,curp+P:curp+P)
    end
    function __parsestr_core_fetchuntil(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Regex,fun::Function=x->x,defval=undef)
        R=match(P,txt,curp)
        if R==nothing
            if defval==undef
                @__ferror__tmp("Cannot match $P at position $curp.")
            else
                return (defval, curp, prvp, 0:-1)
            end
        elseif isempty(R.captures)
            (fun(txt[curp:R.offset-1]),R.offset+length(R.match),curp,R.offset:length(R.match)-1)
        else
            (fun(txt[curp:onlyone(R.offsets)-1]),R.offsets[1]+length(R.captures[1]),curp,R.offsets[1]:R.offsets[1]+length(R.captures[1])-1)
        end
    end
    function __parsestr_core_fetchto(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Regex,fun::Function=x->x,defval=undef)
        R=match(P,txt,curp)
        if R==nothing
            if defval==undef
                @__ferror__tmp("Cannot match $P at position $curp.")
            else
                return (defval, curp, prvp, 0:-1)
            end
        elseif isempty(R.captures)
            (fun(txt[curp:R.offset+length(R.match)-1]),R.offset+length(R.match),curp,R.offset:length(R.match)-1)
        else
            (fun(txt[curp:onlyone(R.offsets)+length(R.captures[1])-1]),R.offsets[1]+length(R.captures[1]),curp,R.offsets[1]:R.offsets[1]+length(R.captures[1])-1)
        end
    end
    function __parsestr_core_rfetch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::AbstractString,fun::Function=x->x,defval=undef)
        ncurp=first(mrg)-1
        R=rsearch(txt,P,first(mrg)-1)
        if isempty(R) || first(R)<prvp
            if defval==undef
                @__ferror__tmp("Cannot rsearch $P at position $ncurp.")
            else
                return (defval, curp, prvp, 0:-1)
            end
        end
        (fun(txt[last(R)+1:ncurp]),curp,prvp,R)
    end
    function __parsestr_core_rfetch(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,P::Int,fun::Function=x->x)
        R=first(mrg)|>x->x-P:x-1
        (fun(txt[R]),curp,prvp,R)
    end
    function __parsestr_core_back(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,cmd)
        rlt,curp,prvp=__parsestr_core(txt[1:curp-1],prvp,prvp,mrg,cmd)
        (rlt,curp,prvp,mrg)
    end
    function __parsestr_core_before(txt::AbstractString,curp::Int,prvp::Int,mrg::Range,cmd)
        rlt,_,_=__parsestr_core(txt[1:curp-1],prvp,prvp,mrg,cmd)
        (rlt,curp,prvp,mrg)
    end
    function __parsestr_core_end(txt::AbstractString,curp::Int,prvp::Int,mrg::Range)
        (nothing,length(txt),curp,length(txt)+1:length(txt)+1)
    end
end

function parsestr1(txt::AbstractString, arg...; default=undef, kw...)
    rlt=parsestr(txt, arg...; no_err=!isa(default, UndefInitializer), kw...)
    if isnothing(rlt)
        default
    else
        only(rlt)
    end
end
parsestr1(txtvec::AbstractVector{<:AbstractString}, arg...; kw...)=map(x->parsestr1(x, arg...; kw...), txtvec)

export ParseStr, next, next1

struct ParseStr
    str::AbstractString
    start_ref::Vector{Int}
    prev_ref::Vector{Int}
    lastMatchRange_ref::Vector{Range}
end

ParseStr(txt::AbstractString; start::Integer=1)=ParseStr(txt,[start],[1],[0:-1])
# import Base.next
next(obj::ParseStr, args...; kw...)=parsestr(obj.str, args...; start_ref=obj.start_ref, prev_ref=obj.prev_ref, lastMatchRange_ref=obj.lastMatchRange_ref, kw...)
function next1(obj::ParseStr, args...; default=undef, kw...)
    rlt=next(obj, args...; no_err=!isa(default, UndefInitializer), kw...)
    if isnothing(rlt)
        default
    else
        only(rlt)
    end
end
#}}

#{{ RefAny
# R = RefAny([Init_value = nothing])
# Create a reference. R.x is the value.
# The difference between Base.Ref is RefAny is not for C-interfact purpose. So R.x can also change its type in program.
# See also: Ref
# Jieyi, Xiong, Dec 30, 2015

mutable struct  RefAny
    x
end
RefAny()=RefAny(nothing)
export RefAny
#}}

#{{ @repress_stdout @repress_stderr
#[[ @repress_stdout @repress_stderr ]]
# @repress_stdout|err any_command
#Repress stdout/stderr printing on screen.
#See also: capstdout, @timetip
#Xiong Jieyi, 17 Jun 2019

export @repress_stdout, @repress_stderr
macro repress_stdout(cmd)
    quote
        __oLDsTDoUT__=stdout
        try
            redirect_stdout()
            $cmd
        catch err
            redirect_stdout(__oLDsTDoUT__)
            rethrow(err)
        finally
            redirect_stdout(__oLDsTDoUT__)
        end
    end |> esc
end
macro repress_stderr(cmd)
    quote
        __oLDsTDeRR__=stderr
        try
            redirect_stderr()
            $cmd
        catch err
            redirect_stderr(__oLDsTDeRR__)
            rethrow(err)
        finally
            redirect_stderr(__oLDsTDeRR__)
        end
    end |> esc
end
#}}

#{{ @timetip
#[[ @timetip ]]
# @timetip[ Seconds=300] command
# Run command if time lapsed 300 sec. It could be used inside loop to print process tip.
# See also: @nowarn, capstdout
# Xiong Jieyi, 12 Jun 2019

__lAST_tIME_tIP__=0.0
function _reset_last_time_tip()
    global __lAST_tIME_tIP__
    __lAST_tIME_tIP__=time()
end
export @timetip
macro timetip(laps::Real, X)
    # @eval(Main, __lAST_tIME_tIP__=time())
    laps=Float64(laps)
    quote
        # global __lAST_tIME_tIP__
        if time()-my.__lAST_tIME_tIP__::Float64 > $laps
            my._reset_last_time_tip()
            # __lAST_tIME_tIP__=time()
            print(my.Dates.format(my.Dates.now(), "HH:MM")*"> ")
            $X
        end
    end |> esc
end
macro timetip(X)
    # @eval(Main, __lAST_tIME_tIP__=time())
    quote
        # global __lAST_tIME_tIP__
        if time()-my.__lAST_tIME_tIP__::Float64 > 300.0
            my._reset_last_time_tip()
            # __lAST_tIME_tIP__=time()
            print(my.Dates.format(my.Dates.now(), "HH:MM")*"> ")
            $X
        end
    end |> esc
end
#}}

