#{{ tb ds tbexp
#[[ tb tb! ds ds! tbexp tbexp! ]]
#Dict{String,Any}(...)=tb|ds([Dict, ]"key1"=>val1, "key2"=>val2; key3=val3, key4=val4, ...)
# ... = tb|ds([Dict, ]c"key1, key2, ..." => (val1, val2, ...), ...; ...)
# ... = tb|ds([Dict, ]("key1"|:, nothing, ("key2", "key3")) => (value1|Dict, ignored, (value2, value3)), ...; ...)
# ... = tb|ds(A, B; ...) = tb|ds(A => B; ...) #Only support 2-3 input.
# ... = tb!|ds!(Dict, ...; ...) #The inputed Dict will be modified.
# ... = (value1, value2, ...) |> tb|ds(c"key1, key2|_, ...") #tb|ds output a function for pipeline convenience. Recur supported also.
# ... = theOnlyValue |> tb|ds("theOnlyKey")
# ... = tbexp|tbexp!( multi_rows_table, one_row_tb_params...;...) # Row-repeat the all-but-first inputed table and merge with the first inputs. At least two inputs. Conflicted keys will be covered by the rightest one.
# ... = tbexp( one_row_table, multi_rows_tb_params...;...) # [ Abolished with warning! ] Row-repeat the first inputed table and merge with other inputs. At least two inputs. Conflicted keys will be covered by the rightest one.
# ... = tbexp(("field"=>value, ...); ...) | tbexp((field=value, ...); ...) #Also acceptable.
#Create a table-style Dict (Dict{String,Any}()). In tb, the rownum will be checked but not in ds.
#In the `=>' mode, if any key is a Colon (:) rather than a string, the value should be a Dict, and it will be merged to the table by dictconc; if any key is nothing rather than a string, the paired value will be ignored, and its row number will not be checked.
#In tbexp(), if any field names are exists in both repeating-data and non-repeating-data, only the non-repeating-data will be used, and no error or warning will occured. However, using this feature is NOT recommended as all the data will be repeated firstly without checking the fieldnames.
#For tb|tbexp(T1::Dict, key=T2), the row consistency for both T1 and T2 will be checked. However, for tb!|tbexp!(T1::Dict, key=T2), the row consistency of T1 will NOT be checked. You can use tb!(T, key=value) as a safe way of T["key"]=value.
#See also: dictfun, reprow, tb2DataFrame, pandasDF, rec, fd*
#Xiong Jieyi, February 9, 2015 >May 27, 2015>Jun 3, 2015>Nov 4, 2015>Nov 4, 2015>23 Jul 2017>18 Dec 2018>18 Apr 2020>13 Jul 2020>25 Apr 2022>26 Oct 2023

export tb, ds, tb!, ds!, tbexp, tbexp!

function __deepzip__(P::Pair)
    #7 May 2022: to solve bug tb(split("a,b,c", ','), 1:3)
    # if isa(P[1], Union{AbstractString, Colon})
    #     return [convert(Pair{Union{String, Colon}, Any}, P)]
    # end
    if isa(P[1], AbstractString)
        if isa(P[1], String)
            return [convert(Pair{Union{String, Colon}, Any}, P)]
        else
            return [convert(Pair{Union{String, Colon}, Any}, string(P[1])=>P[2])]
        end
    elseif isa(P[1], Colon)
        return [convert(Pair{Union{String, Colon}, Any}, P)]
    end
    isa(P[1], Union{AbstractVector, Tuple}) || error("Key of a Table should always be String.")
    length(P[1])==length(P[2]) || error("Two inputs cannot be paired due to different lengths or strcutures.")
    O=Pair{Union{String, Colon}, Any}[]
    if isa(P[1], Tuple)
        for (a, b) in zip(P[1], P[2])
            if !isnothing(a) #Added in 25 Oct 2023
                append!(O, __deepzip__(a=>b))
            end
        end
    else
        for (a, b) in zip(P[1], P[2])
            # if a!="_" #It may not be a good idea to support "_" as an ignore marker while all other strings not.
            append!(O, __deepzip__(a=>b))
            # end
        end
    end
    O
end

function tb!(D::Dict{T, Any}, args::Pair...; wargs...) where T<:AbstractString
    if isempty(D)
        len=-1
    else
        # len=rnum(D)
        len=rownum(D) #Modified in 25 Apr 2022
    end
    for carg in args
        for (k, v) in __deepzip__(carg)
            if len==-1
                len=rnum(v)
            else
                len==rnum(v) || error("Inputs $k doesn't have the same row number $len.")
            end
            if isa(k, Colon)
                isa(v, Dict) || error("For tb!(..., : => V, ...), v should be a Dict.")
                dictconc!(D, v)
            else
                D[k]=v
            end
        end
    end
    for (k, v) in wargs
        if len==-1
            len=rnum(v)
        else
            len==rnum(v) || error("Inputs $k doesn't have the same row number $len.")
        end
        D[string(k)]=v
    end
    D
end

tb(fds::Union{AbstractVector,Tuple}, X::Union{Tuple, AbstractVector}; warg...)=tb!(Dict{String, Any}(), fds => X; warg...)
tb(K::Pair, args::Pair...; wargs...)=tb!(Dict{String, Any}(), K, args...; wargs...)
tb(; wargs...)=tb!(Dict{String, Any}(); wargs...)
function tb(D::Dict, args...; wargs...)
    rnum(D) #To check the row consistency.
    tb!(deepcopy(D), args...; wargs...)
end
tb(fds::Union{AbstractVector,Tuple})=x->tb(fds=>x)
tb(fd::AbstractString)=x->tb(fd=>x)

function tbexp!(R::Dict{Tr, Any}, args...; wargs...) where {Tr<:AbstractString}
    # N=rnum(R)
    N=rownum(R)
    merge!(R, reprow(ds(args...; wargs...), N))
end
function tbexp(R::Dict{Tr, Any}, args...; wargs...) where {Tr<:AbstractString}
    T=ds(args...; wargs...)
    N=rownum(T)
    if N>1
        @warn "tbexp( one_row_table, multi_rows_tb_params...;...) will not be supported in the future. Using tbexp(multi_rows_table, one_row_tb_params...;...) instead."
        merge!(reprow(R, N), T) #Shall I support it?
    else
        N=rnum(R)
        merge(R, reprow(T, N))
    end
end
tbexp(R::NamedTuple, arg...; kw...)=tbexp(ds(;R...), arg...; kw...)
tbexp(R::Tuple, arg...; kw...)=tbexp(ds(R...), arg...; kw...)
# function tbexp(R::Dict{T, Any}, TB1::Dict{T2, Any}, args::Pair...; wargs...) where {T<:AbstractString, T2<:AbstractString}
#     N=rownum(TB1)
#     if N>1
#         tb!(merge!(reprow(R, N), TB1), args...; wargs...)
#     else
#         N=rnum(R)
#         # tbexp!(deepcopy(R), TB1, args...; wargs...)
#         merge(R, reprow(ds(TB1, args...; wargs...), N))
#     end
# end
# function tbexp(R::Dict{T, Any}, TB1::Pair, args::Pair...; wargs...) where T<:AbstractString
#     N=rownum(TB1.second)
#     if N>1
#         tb!(reprow(R, N), TB1, args...; wargs...)
#     else
#         N=rnum(R)
#         # tbexp!(deepcopy(R), TB1, args...; wargs...)
#         merge(R, reprow(ds(TB1, args...; wargs...), N))
#     end
# end
# function tbexp(R::Dict{T, Any}; wargs...) where T<:AbstractString
#     N=rownum(wargs[1])
#     if N>1
#         tb!(reprow(R, N); wargs...)
#     else
#         N=rnum(R)
#         # tbexp!(deepcopy(R), ; wargs...)
#         merge(R, reprow(ds(wargs...), N))
#     end
# end

function ds!(D::Dict{T, Any}, args::Pair...; wargs...) where T<:AbstractString
    for carg in args
        for (k, v) in __deepzip__(carg)
            if isa(k, Colon)
                isa(v, Dict) || error("For ds!(..., : => V, ...), v should be a Dict.")
                dictconc!(D, v)
            else
                D[k]=v
            end
        end
    end
    for (k, v) in wargs
        D[string(k)]=v
    end
    D
end

ds(K::Pair, args...; wargs...)=ds!(Dict{String, Any}(), K, args...; wargs...)
ds(; wargs...)=ds!(Dict{String, Any}(); wargs...)
ds(D::Dict,args...;wargs...)=ds!(deepcopy(D),args...;wargs...)
ds(fds::Union{AbstractVector,Tuple}, X::Union{Tuple, AbstractVector}; warg...)=ds!(Dict{String, Any}(), fds => X; warg...)
ds(fds::Union{AbstractVector,Tuple})=x->ds(fds=>x)
ds(fd::AbstractString)=x->ds(fd=>x)
#}}

#{{ dictfun dictfun!
#[[ dictfun dictfun! ]]
#outD =  dictfun(func, D::Dict...; key=keys(D1)|:shared, copy=false, inkey=false)
#  D1 = dictfun!(func, D::Dict...; key=keys(D1)|:shared, replace=true, inkey=false)
#  DF = dataframefun(func, DF::DataFrame...; inkey=false)
#Run function on each value of Dict under the given keys. Note that both dictfun and dictfun! will change the value of input Dict if the given function do that.
#The only difference of these two function is dictfun collect the output of function as a new dict and return it, while dictfun! just simplely return the first argument.
#If copy=true, the missing D-fields in key will be deepcopy to outD.
#If replace=true, D1[key] will be assigned to function output.
#If key=:shared, only the shared keys among all inputed Dict are calculated. Note that Dict number should be at least two.
#If inkey=true, the first input of called function is the key.
#Change in 16 May 2019: changed the default value of dectfun!( replace=false) to true.
#Xiong Jieyi, 13 May 2014>11 Sep 2014>October 14, 2014>Nov 25, 2015>Jul 14, 2016>28 Jan 2019>16 May 2019>29 Mar 2021

function dictfun(f::Function,D::Dict{Tk,}; key=keys(D), copy::Bool=false, inkey::Bool=false) where {Tk}
    if isa(key, Symbol)
        if key==:shared
            error("key=:shared is not supported when only one Dict inputted.")
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    O=Dict{Tk,Any}()
    if inkey
        for k in key
            O[k]=f(k, D[k])
        end
    else
        for k in key
            O[k]=f(D[k])
        end
    end
    if copy
        for k in keys(D)
            if !in(k,key)
                O[k]=deepcopy(D[k])
            end
        end
    end
    return O
end
function dictfun(f::Function,args::Dict...;key=keys(args[1]), copy::Bool=false, inkey::Bool=false)
    if isa(key, Symbol)
        if key==:shared
            key=intersect(keys.(args)...)
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    O=Dict{keytype(args[1]),Any}()
    if inkey
        for k in key
            O[k]=f(k, map(x->x[k],args)...)
        end
    else
        for k in key
            O[k]=f(map(x->x[k],args)...)
        end
    end
    if copy
        for k in keys(args[1])
            if !in(k,key)
                O[k]=deepcopy(args[1][k])
            end
        end
    end
    return O
end
function dataframefun(f::Function, args::AbstractDataFrame...; inkey::Bool=false)
    O=DataFrame()
    nms=names(args[1])
    if inkey
        for k in nms
            O[k]=f(k, map(x->x[k],args)...)
        end
    else
        for k in nms
            O[k]=f(map(x->x[k],args)...)
        end
    end
    return O
end

function dictfun!(f::Function, D::Dict; key=keys(D), replace::Bool=true, inkey::Bool=false)
    if isa(key, Symbol)
        if key==:shared
            error("key=:shared is not supported when only one Dict inputted.")
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    if replace
        if inkey
            for k in key
                D[k]=f(k, D[k])
            end
        else
            for k in key
                D[k]=f(D[k])
            end
        end
    else
        if inkey
            for k in key
                f(k, D[k])
            end
        else
            for k in key
                f(D[k])
            end
        end
    end
    return D
end
function dictfun!(f::Function, args::Dict...; key=keys(args[1]), replace::Bool=true, inkey::Bool=false)
    if isa(key, Symbol)
        if key==:shared
            key=intersect(keys.(args)...)
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    if replace
        if inkey
            for k in key
                args[1][k]=f(k, map(x->x[k],args)...)
            end
        else
            for k in key
                args[1][k]=f(map(x->x[k],args)...)
            end
        end
    else
        if inkey
            for k in key
                f(k, map(x->x[k],args)...)
            end
        else
            for k in key
                f(map(x->x[k],args)...)
            end
        end
    end
    return args[1]
end
export dictfun, dictfun!, dataframefun
#}}

#{{ istable [private]
# istable(Dict)
# ** This function is nolonger expoerted. Use isa(x, Dict) && isgroup(x) instead. * 29 Sep 2020 **
#If the Dict is a table, i.e., have the same row number for each elements, return true. Otherwise, return false.
#See also: vcatr, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dictconc, dictconc!, tbcomp
#Xiong Jieyi, 13 May 2014 > 29 Sep 2020

function istable(D::Dict)
    return iseqrownum(D)
    # rn=-1;
    # for (k,v) in D
    #     if rn==-1
    #         rn=rownum(v)
    #     elseif rn!=rownum(v)
    #         return false
    #     end
    # end
    # return true
end
# export istable
#}}

#{{ fd fd_i fd! fd_i!
#[[ fd fd_i fd! fd_i! ]]
#sub_Dict = fd|fd!(Dict, fields::AnyIterateable|Regex|AbstractString)
#Or ...= fd_i|fd_i!(Dict, fields::AnyIterateable|Regex|AbstractString)
#Or fd|fd_i|fd_i!(Dict, X1, X2...) #Union or substract all given field|Regex. [*]
#Get structure with subfields. fd_i only reture the Dict with the feilds not included or matched the Regex. Both funcitons do NOT deepcopy the fields.
#The ! function will also change the input Dict.
#The difference between delete!(D::Dict,f::AbstractString) and fd_i!(D::Dict,f::AbstractString): fd_i! will check if given field exists in the Dict, but delete! will not. fd_i(D, f::AbstractString) will also check the field existence. However, when the 2nd input is not a string, the field existence will not be checked.
#[*] Try to avoid fd(D,c"a",c"b",...) but use fd(D, c"a,b,...") instead for speed reason. You can use fd(D, c"fields, ...", r"Regex1", r"Regex2"...). fd! doesn't support this syntax.
#See also: sk, rename, dictfun
#Xiong Jieyi, 3 Sep 2014 > February 10, 2015 >Oct 9, 2015 >Aug 12, 2016 >4 Mar 2020

import Base.fd #fd has already been used in Base but for other propose.
export fd, fd_i, fd!, fd_i!
function fd(D::Dict{Tk, Tv}, F) where {Tk<:AbstractString, Tv}
    O=Dict{Tk, Tv}()
    for f in F
        O[f]=D[f]
    end
    return O
end
fd(D::Dict{Tk,}, f::AbstractString) where {Tk<:AbstractString} =fd(D, [f])
function fd(D::Dict{Tk,}, R::Regex) where {Tk<:AbstractString}
    filter(((k, _),)->occursin(R, k), D)
end
function fd_i(D::Dict{Tk,}, F) where {Tk<:AbstractString}
    K=collect(keys(D))
    fd(D,setdiff(K,F))
end
function fd_i(D::Dict{Tk,}, f::AbstractString) where {Tk<:AbstractString}
    haskey(D, f) || error("Field $f is not exist.")
    fd_i(D, [f])
end
function fd_i(D::Dict{Tk,}, R::Regex) where {Tk<:AbstractString}
    filter(((k, _),)->!occursin(R, k), D)
end
function fd_i!(D::Dict{Tk,}, F) where {Tk<:AbstractString}
    for f in F
        haskey(D, f) && delete!(D, f)
    end
    D
end
function fd!(D::Dict{Tk,}, F) where {Tk<:AbstractString}
    sF=Set(F)
    for f in keys(D)
        if !in(f, sF)
            delete!(D, f)
        end
    end
    length(D)==length(sF) || error("Not all given fields are in table.")
    D
end
fd!(D::Dict{Tk,}, f::AbstractString) where {Tk<:AbstractString} = fd!(D, [f])
function fd!(D::Dict{Tk,}, R::Regex) where {Tk<:AbstractString}
    filter!(((k, _),)->occursin(R, k), D)
end
function fd_i!(D::Dict{Tk,}, R::Regex) where {Tk<:AbstractString}
    filter!(((k, _),)->!occursin(R, k), D)
end
function fd_i!(D::Dict{Tk,}, f::AbstractString) where {Tk<:AbstractString}
    haskey(D,f) || error("Field $f is not exist.")
    delete!(D,f)
end
fd_i(D::Dict{Tk,}, X, P...) where {Tk<:AbstractString} = fd_i(fd_i(D, X), P...)
fd_i!(D::Dict{Tk,}, X, P...) where {Tk<:AbstractString} = fd_i!(fd_i!(D, X), P...)
fd(D::Dict{Tk,}, P...) where {Tk<:AbstractString} = merge!(map(x->fd(D, x), P)...)
#}}

#{{ @d_str

#[[ @d_str  ]]
#d"field1[, field2, ...]"dict or d"""..."""dict
#return (dict["field1"], dict["field2"],...)
#field will be strip.
#d"field"dict equal to dict["field"]
#For the simplicity to form a tuple with the table fields.
#See also: @c_str, @dd_str, dict2tuple, @with
#Xiong Jieyi, 9 Sep 2014>January 6, 2015>Oct 19, 2015>May 13, 2016

macro d_str(S,P...)
    fields=map(strip,split(S,','))
    # for i=1:length(fields)
    #     fields[i]=replace(fields[i],r"^[\ \t]*\n","")
    # end
    cmd=if isempty(P)
        "D->("*join(["""D["$f"]""" for f in fields],',')*")"
    else
        dct=onlyone(P)
        join(["""$dct["$f"]""" for f in fields],',')
    end
    esc(Meta.parse(cmd))
end
export @d_str

# if VERSION<v"0.4"
#     macro d_mstr(fields,dct)
#         fields=map(lstrip,split(fields,','))
#         esc(parse(join(["""$dct["$f"]""" for f in fields],',')))
#     end
#     export @d_mstr
# end
#}}

#{{ @dd_str @dd_i_str
#[[ @dd_str @dd_i_str ]]
# SubDict = dd" field1[ => rename1], field2[ => prefix_# ],..."Dict # fd and rename! a Dict.
# SubDict = dd_i" remove_field, ..., rename_field => ..."Dict       # fd_i and rename! a Dict.
# Or ...  = dd|dd_i"..."(Dict) # Used as a function.
# fd and rename! Dict.
# If # is in new field name but not in the old one, it will be replaced by the old field name. All the blanks at the flanks of names will be stripped.
#See also: @d_str, fd*, rename!, @with
#Xiong Jieyi, May 13, 2016 > 5 Oct 2020

export @dd_str, @dd_i_str, @DD_str
macro dd_str(S::String, P...)
    fds1=AbstractString[]
    fds2=AbstractString[]
    fds3=AbstractString[]
    for C in split(S,',')
        C=strip(C)
        tC=split(C,"=>")
        if length(tC)>1
            tC[1]=strip(tC[1])
            tC[2]=strip(tC[2])
            push!(fds1,tC[1])
            push!(fds2,tC[1])
            if !in('#',tC[1])
                tC[2]=replace(tC[2], "#" => tC[1])
            end
            push!(fds3,tC[2])
        else
            push!(fds1,C)
        end
    end
    cmd=if isempty(P)
        if isempty(fds2)
            "D->fd(D,$fds1)"
        else
            "D->rename!(fd(D,$fds1),$fds2,$fds3)"
        end
    else
        D=onlyone(P)
        if isempty(fds2)
            "fd($D,$fds1)"
        else
            "rename!(fd($D,$fds1),$fds2,$fds3)"
        end
    end
    esc(Meta.parse(cmd))
end
macro dd_i_str(S::String, P...)
    fds1=AbstractString[]
    fds2=AbstractString[]
    fds3=AbstractString[]
    for C in split(S,',')
        C=strip(C)
        tC=split(C,"=>")
        if length(tC)>1
            tC[1]=strip(tC[1])
            tC[2]=strip(tC[2])
            # push!(fds1,tC[1])
            push!(fds2,tC[1])
            if !in('#',tC[1])
                tC[2]=replace(tC[2], "#" => tC[1])
            end
            push!(fds3,tC[2])
        else
            push!(fds1,C)
        end
    end
    cmd=if isempty(P)
        if isempty(fds2)
            "D->fd_i(D,$fds1)"
        else
            "D->rename!(fd_i(D,$fds1),$fds2,$fds3)"
        end
    else
        D=onlyone(P)
        if isempty(fds2)
            "fd_i($D,$fds1)"
        else
            "rename!(fd_i($D,$fds1),$fds2,$fds3)"
        end
    end
    esc(Meta.parse(cmd))
end
macro DD_str(S, P...)
    error("@DD_str is renamed to @dd_i_str.")
end
#}}

#{{ @with @withlet @withrow @filterrow
#[[ @with @withlet @withrow @filterrow ]]
#... = @with|@withlet dict any-cmd-with-$field-or-$"field"-or-$_
#  A lasy way to use dict["field"]. In the second parameter, any $field or $"field" will be replaced as dict["field"]. $_ will be replaced as dict itself. `dict' must be in string-like keys. In the nested @with, $$abc will be parsed to $abc.
#  @withlet: the 1st input can also be an expression, and the commands will be run in a let...end block.
#  @with: the 1st input can only be a variable name.
#anonymous_fun = @with( (..., $_, ...) -> any-cmd-with-$field-or-$"field"-or-$_ )
#generate = @with( any-cmd-with-$field-or-$"field"-or-$_ for $_|(..., $_, ...) =|in ...)
#vec = [ @with(above arguments...)... ]
#... = @with fun(...) do ..., $_, ... end [ |> ... ] #followed pipe is not effected by @with.
#Or @with for $_|(..., $_, ...) =|in eachr(...)
#      any-cmd-with-$field-or-$"field"-or-$_
#  end
#... = @withrow table any-cmd-with-$field-or-$"field"-or-$_-or-$^
#  The same as ... = rowfun(T) do X
#                  @with(X, cmdWith($field, $_, $^))
#              end
#  If $_ is not used, function will do fd(T, used_keys) first to speed up the rowfun().
#  $^ is the current row number.
#... = @filterrow table any-cmd-with-$field-or-$"field"-or-$_-or-$^
#  Filter table by row, like rec(table, @withrow(table, ...)). When table is zero-row, it always output deepcopy(table).
# For @withrow and @filterrow, the 1st input can either be a variable name or be an expression.
# For all macros, the first input could be either a variable name or an expression.
#See also: @d_str, @dd_str, @rec
#Xiong Jieyi, 27 Oct 2017 > 7 Jan 2022 > 15 Feb 2022 > 5 Jan 2023 >17 Jul 2024

macro with(T::Symbol, X::Expr)
    function recur(A::Vector, I::Int)
        C=A[I]
        isa(C, Expr) || return
        if isa(C, Symbol)
            return
        elseif C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            if C.args[1]==:_
                A[I]=T
            else
                C.head=:ref
                C.args=Any[T, string(C.args[1])]
            end
        elseif C.head==:$ && length(C.args)==1 && isa(C.args[1], Expr) && C.args[1].head==:$ && length(C.args[1].args)==1 && (isa(C.args[1].args[1], Symbol) || isa(C.args[1].args[1], String))
            # In nested @with, $$xxx->$xxx. 17 Jul 2024
            C.args=C.args[1].args
        else
            for i=1:length(C.args)
                recur(C.args, i)
            end
        end
    end
    for i=1:length(X.args)
        recur(X.args, i)
    end
    esc(X)
end
macro with(X::Expr)
    internalT=:_witH_insidE_variablE_dO_noT_usE_
    flag=false
    function recur(A::Vector, I::Int, selfTypeCheck::Bool, childTypeCheck::Bool)
        C=A[I]
        isa(C, Expr) || return
        if isa(C, Symbol)
            return
        elseif C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            if C.args[1]==:_
                flag=true
                if selfTypeCheck
                    C.head=:(::)
                    C.args=Any[internalT, :(Dict{<:AbstractString,})]
                else
                    A[I]=internalT
                end
            elseif C.head==:$ && length(C.args)==1 && isa(C.args[1], Expr) && C.args[1].head==:$ && length(C.args[1].args)==1 && (isa(C.args[1].args[1], Symbol) || isa(C.args[1].args[1], String))
            # In nested @with, $$xxx->$xxx. 17 Jul 2024
                C.args=C.args[1].args
            else
                C.head=:ref
                C.args=Any[internalT, string(C.args[1])]
            end
        else
            for i=1:length(C.args)
                recur(C.args, i, childTypeCheck, childTypeCheck)
            end
        end
    end
    if X.head==:->
        if isa(X.args[1], Expr) && X.args[1].head==:tuple
            for i=1:length(X.args[1].args)
                recur(X.args[1].args, i, true, false)
            end
        else
            recur(X.args, 1, true, false)
        end
        flag || error("No `\$_' was found in the argument list.")
        recur(X.args, 2, false, false)
    elseif X.head==:do && X.args[2].head==:->
        if isa(X.args[2].args[1], Expr) && X.args[2].args[1].head==:tuple
            for i=1:length(X.args[2].args[1].args)
                recur(X.args[2].args[1].args, i, true, false)
            end
        else
            recur(X.args[2].args, 1, true, false)
        end
        flag || error("No `\$_' was found in the argument list.")
        recur(X.args[2].args, 2, false, false)
    elseif X.head==:call && X.args[1]==:|> &&
        X.args[2].head==:do && X.args[2].args[2].head==:->
        #For the style @with fun(...) do $_ ... end |> ...
        if isa(X.args[2].args[2].args[1], Expr) && X.args[2].args[2].args[1].head==:tuple
            for i=1:length(X.args[2].args[2].args[1].args)
                recur(X.args[2].args[2].args[1].args, i, true, false)
            end
        else
            recur(X.args[2].args[2].args, 1, true, false)
        end
        flag || error("No `\$_' was found in the argument list.")
        recur(X.args[2].args[2].args, 2, false, false)
    elseif X.head==:for && length(X.args)==2 &&
        X.args[1].head==:(=) && X.args[2].head==:block
        recur(X.args[1].args, 1, true, false)
        flag || error("No `for \$_ =...' was found.")
        for i=1:length(X.args[2].args)
            recur(X.args[2].args, i, false, false)
        end
    elseif X.head==:generator && length(X.args)==2 &&
        X.args[2].head==:(=) && length(X.args[2].args)==2
        #For style [@with($a for $_ in eachr(...))...]
        recur(X.args[2].args, 1, true, false)
        flag || error("No `for \$_ =...' was found.")
        recur(X.args, 1, false, false)
    else
        error("No `->`, `for' or `do' syntax was found.")
    end
    esc(X)
end
macro withlet(T::Union{Expr, Symbol}, X::Expr)
    internalT=:_withroW_insidE_variablE_dO_noT_usE_
    quote
        let $internalT=$T
            @with($internalT, $X)
        end
    end |> esc
end
macro withrow(T::Union{Symbol, Expr}, X::Expr)
    internalT=:_withroW_insidE_variablE_dO_noT_usE_
    internalN=:_withroW_insidE_variablE_dO_noT_usE_rownuM_
    used_keys=String[]
    used_whole_table=false
    used_row_num=false
    function recur(A::Vector, I::Int)
        C=A[I]
        isa(C, Expr) || return
        if C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            if C.args[1]==:_
                used_whole_table=true
                # C.head=:(::)
                # C.args=Any[:_withroW_insidE_variablE_dO_noT_usE_, :(Dict{<:AbstractString,})]
                A[I]=internalT
            elseif C.args[1]==:^
                used_row_num=true
                C.head=:(::)
                C.args=Any[internalN, :Int]
            else
                push!(used_keys, string(C.args[1]))
                C.head=:ref
                C.args=Any[internalT, string(C.args[1])]
            end
        elseif C.head==:$ && length(C.args)==1 && isa(C.args[1], Expr) && C.args[1].head==:$ && length(C.args[1].args)==1 && (isa(C.args[1].args[1], Symbol) || isa(C.args[1].args[1], String))
            # In nested @with, $$xxx->$xxx. 17 Jul 2024
            C.args=C.args[1].args
        else
            for i=1:length(C.args)
                recur(C.args, i)
            end
        end
    end
    for i=1:length(X.args)
        recur(X.args, i)
    end
    tbcmd=used_whole_table ? T : :(fd($T, $used_keys))
    if used_row_num
        quote
            rowfun($tbcmd, no=true) do $internalN, $internalT
                $X
            end
        end
    else
        quote
            rowfun($tbcmd) do $internalT
                $X
            end
        end
    end |> esc
end
macro filterrow(T::Symbol, X::Expr)
    quote
        if rnum($T::Dict{<:AbstractString,})==0
            deepcopy($T)
        else
            getrow($T, @withrow($T, ($X)::Bool)::AbstractVector{Bool})
        end
    end |> esc
end
macro filterrow(T::Expr, X::Expr)
    quote
        let _filterroW_insidE_variablE_dO_noT_usE_=$T
            @filterrow(_filterroW_insidE_variablE_dO_noT_usE_, $X)
        end
    end |> esc
end
export @with, @withlet, @withrow, @filterrow
#}}

#{{ rec rec_i
#[[ rec rec_i @rec @rec_i ]]
#Table = rec(Table, index)
#Table = rec_i(Table, index)
# ...  = @rec[_i](Table, $field.>0)
#Get the rows of each element of Dict (rec), or remove the rows of Dict (rec_i)
#The difference between rec and getrow: rec will check if the input is a table firstly, but getrow will not.
#17 Jul 2024: For @rec[_i] nested with @withlet, $$abc will be parsed to $abc similar as @with[let].
#See also: getrow, fd, fd_i, ft, sk, dictconc, tbx, dict2tuple
#Xiong Jieyi, 8 Jul 2014>28 Sep 2014>October 16, 2014>14 Nov 2017

function rec(D::Dict,I::Union{Int,AbstractVector})
    istable(D) || error("Input is not a table.")
    getrow(D,I)
end
rec_i(D::Dict,I::BitVector)=rec(D, .!I)
rec_i(D::Dict,I::AbstractVector{Bool})=rec(D, .!I)
rec_i(D::Dict,I::AbstractVector)=rec(D,fill(true,rownum(D))|>x->(x[I].=false;x))
rec_i(D::Dict,i::Integer)=rec(D,fill(true,rownum(D))|>x->(x[i]=false;x))

macro rec(T, C)
    function recur(C)
        if C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            C.head=:ref
            C.args=Any[T,string(C.args[1])]
        elseif C.head==:$ && length(C.args)==1 && isa(C.args[1], Expr) && C.args[1].head==:$ && length(C.args[1].args)==1 && (isa(C.args[1].args[1], Symbol) || isa(C.args[1].args[1], String))
            # In nested @with, $$xxx->$xxx. 17 Jul 2024
            C.args=C.args[1].args
        else
            for (i, v) in enumerate(C.args)
                if isa(v, Expr)
                    recur(v)
                end
            end
        end
    end
    isa(C, Expr) && recur(C)
    esc(:(rec($T, $C)))
end

macro rec_i(T, C)
    function recur(C)
        if C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            C.head=:ref
            C.args=Any[T,string(C.args[1])]
        elseif C.head==:$ && length(C.args)==1 && isa(C.args[1], Expr) && C.args[1].head==:$ && length(C.args[1].args)==1 && (isa(C.args[1].args[1], Symbol) || isa(C.args[1].args[1], String))
            # In nested @with, $$xxx->$xxx. 17 Jul 2024
            C.args=C.args[1].args
        else
            for (i, v) in enumerate(C.args)
                if isa(v, Expr)
                    recur(v)
                end
            end
        end
    end
    isa(C, Expr) && recur(C)
    esc(:(rec_i($T, $C)))
end

export rec, rec_i, @rec, @rec_i
#}}

#{{ dictconc dictconc!
#[[ dictconc dictconc! ]]
# CombindDict = dictconc(A::Dict, B::Dict; delconflict=false | aprefix="Pre-"|function, bprefix=... | anest="akey", bnest="bkey")
# ... = dictconc(A::Dict, B::Dict, C::Dict, ...; delconflict=false) #Other keywords are not supported.
#   A = dictconc!(A::Dict, ...; ...) #Directly merge other dicts to A.
#Combind all the fields of two dictionaries. If two field's values with the same fieldnames are inconsistent, a error will occur in default, or this field will be deleted (delconflict=true) or renamed (a/bprefix=...). aprefix/bprefix can also be a function, similar as rename!(function, ...). aprefix/bprefix can only be used when input Dict number = 2. If anest or bnest was assigned, the conflict keys will be nested keys with the given names.
#In dictconc!, A will be added B's fields and be outputed, but B is ALWAYS intact (even when bprefix=... is assigned).
#See also: vcatr, vcatr_with, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dict2tuple, tbcomp, rename!, dictshare
#Xiong Jieyi, 16 May 2014 >May 22, 2015 >Oct 9, 2015>Apr 27, 2016>13 Jan 2018>3 May 2019>29 Jan 2021

function dictconc!(O::Dict,B::Dict; delconflict::Bool=false,
                   aprefix::Union{AbstractString, Function}="",
                   bprefix::Union{AbstractString, Function}="",
                   anest::AbstractString="",
                   bnest::AbstractString="")
    confkys=collect(filter(ky->!isequal(O[ky], B[ky]), intersect(keys(O), keys(B))))
    if !isempty(confkys)
        if delconflict
            for ky in confkys
                delete!(O, ky)
            end
        else
            changed_flag=false
            if !isempty(aprefix)
                isempty(anest) || error("aprefix=... and anest=... cannot be used simultaneously.")
                changed_flag=true
                for ky in confkys
                    nky=isa(aprefix, AbstractString) ? aprefix*ky : aprefix(ky)
                    (haskey(O, nky) || haskey(B, nky)) && error("Conflicted key: \"$nky\"")
                    rename!(O, ky, nky)
                end
            elseif !isempty(anest)
                changed_flag=true
                t=fd(O, confkys)
                for ky in confkys
                    delete!(O, ky)
                end
                (haskey(O, anest) || haskey(B, anest)) && error("Conflicted key: \"$nky\"")
                O[anest]=t
            end
            if !isempty(bprefix)
                isempty(bnest) || error("bprefix=... and bnest=... cannot be used simultaneously.")
                changed_flag=true
                for ky in confkys
                    nky=isa(bprefix, AbstractString) ? bprefix*ky : bprefix(ky)
                    (haskey(O, nky) || haskey(B, nky)) && error("Conflicted key: \"$nky\"")
                    O[nky]=B[ky]
                end
            elseif !isempty(bnest)
                changed_flag=true
                (haskey(O, bnest) || haskey(B, bnest)) && error("Conflicted key: \"$nky\"")
                O[bnest]=fd(B, confkys)
            end
            changed_flag || error(f"The coexisted field \"$1\" is heterogeneous."(first(confkys)))
        end
    end
    merge!(O, fd_i(B, confkys))
    # if !(isempty(aprefix) && isempty(bprefix))
    #     delconflict && error("rename_conflict=true conflicts with field renaming.")
    #     confkys=collect(filter(ky->!isequal(O[ky], B[ky]), intersect(keys(O), keys(B))))
    #     if !isempty(aprefix)
    #         O=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, O, confkys)
    #     end
    #     # if !isempty(bprefix)
    #     #     B=rename!(isa(bprefix, AbstractString) ? x->bprefix*x : bprefix, deepcopy(B), confkys)
    #     # end
    # end
    
    # for (k, v) in B
    #     if haskey(O, k)
    #         if !isequal(O[k], v)
    #             if delconflict
    #                 delete!(O, k)
    #             elseif !isempty(bprefix)
    #                 nk=isa(bprefix, AbstractString) ? bprefix*k : bprefix(k)
    #                 haskey(O, nk) && error("New field $k => $nk is still conflicted.")
    #                 O[nk]=B[k]
    #             else
    #                 error("The coexisted field \"$k\" is heterogeneous.")
    #             end
    #         end
    #     else
    #         O[k]=B[k]
    #     end
    # end
    # return O
end
dictconc!(O::Dict,B::Dict,P::Dict...; delconflict::Bool=false)=dictconc!(dictconc!(O::Dict, B::Dict; delconflict=delconflict), P...; delconflict=delconflict)
dictconc(A::Dict,P::Dict...; warg...)=dictconc!(deepcopy(A),P...; warg...)
export dictconc,dictconc!
#}}

#{{ dictshare
#[[ dictshare ]]
#shared, A_spec, B_spec = dictshare(A::Dict, B::Dict)
#Seperate the shared sub-dictionary and each specific sub-dictionary from two input tables. The values with the same key in A and B will be compared by isequal(), i.e., NaN==NaN. A_spec and B_spec will contain the specific keys in each input dictionary, or the same key but with different values in A and B.
#See also: tbx, dictconc
#Xiong Jieyi, 16 Jun 2020

export dictshare
function dictshare(A::Dict, B::Dict)
    sharedkeys=filter(x->isequal(A[x], B[x]), intersect(keys(A), keys(B)))
    (fd(A, sharedkeys), fd_i(A, sharedkeys), fd_i(B, sharedkeys))
end
#}}

#{{ tbuniq, tbuniq_i
#[[ tbuniq tbuniq_i ]]
# UniqTb=tbuniq(Tb)
# UniqTb=tbuniq(Tb, fields|nothing*[, "field"|fields_vec|Regex => function, ...]; trim|trim_recur=false, check=false, stable=false)
#  ...  =tbuniq(function, Tb, fields, field(_vec)=>:, ...; ...) # ":" will be replaced by the first funciton input.
# UniqTb=tbuniq_i(Tb, ignore_fields)
# Get the unique table by the given fields. If the fields is omitted in tbuniq, function using all the fields.
# If trim, trim_recur or check is true, function will check each of non-key-fields to make sure it can be unique following the key-fields. If failed, function will remove this field (trim=true) or throw an error (check=true).
# trim_recur=true will recursly trim the nested table, while trim=true regards the nested table as a whole.
# When "field"|Vector|Regex=>function is given, these fields will be firstly excluded in the unique step, and then be added as To["field"]=grpfun(fun, ..., Ti["field"]). If any Regex is used, function require either trim(_recur)=true or check=true should also be set just for fool-proofing.
# * Arg#2=nothing means all fields except the ones specificly handled in 3+ args. Only allowed with 3+ args.
#See also: coo2mx, tbcomp, grpfun
#Xiong Jieyi, 30 Aug 2014>7 Oct 2014>Sep 8, 2017>9 May 2018>26 Aug 2019>22 Apr 2020>30 Jul 2021>8 Jul 2024

tbuniq(T::Dict{Tk,}; stable::Bool=false) where {Tk<:AbstractString} =
    rec(T, uniqr(Tuple(values(T)); stable=stable)[1]) 

function tbuniq(T::Dict{Tk,}, fld, Pr::Pair...; stable::Bool=false, trim::Bool=false, trim_recur::Bool=false, check::Bool=false) where {Tk<:AbstractString}
    fg=nothing
    if isempty(Pr)
        isnothing(fld) && error("Arg#2 = nothing is not allowed when no \"field\"=>function behand. Just use tbuniq(T) instead of tbuniq(T, nothing)")
    else
        T0=T
        #30 Jul 2021: Support Regex.
        Pr_fds=map(Pr) do (x, _)
            if isa(x, Regex)
                (trim || check) || error("When Regex=>fun is used, either trim=true or check=true should be set.")
                collect(filter(x::Regex, keys(T0)))
            else
                x::Union{AbstractVector, Tuple, AbstractString}
            end
        end
        T=fd_i(T, vcat(Pr_fds...,))
        fg=isnothing(fld) ? fastgrp(dict2tuple(T); stable=stable) : fastgrp(dict2tuple(T, fld); stable=stable)
    end
    T=if trim || check || trim_recur
        isnothing(fld) && error("Arg#2 = nothing is not allowed when trim=true or check=true.")
        otherfd=setdiff(collect(keys(T)), fld)
        if isnothing(fg)
            fg=fastgrp(dict2tuple(T, fld); stable=stable)
        end
        if trim_recur #Added in 8 Jul 2024
            check && error("check=true conflicts with trim_recur=true.")
            function recurfun(cT::Dict{<:AbstractString,})
                nT=typeof(cT)()
                for (ky, val) in cT
                    if isa(val, Dict{<:AbstractString,}) && !isempty((t=recurfun(val);))
                        nT[ky]=t
                    elseif all(isrowsame, eachgrp(fg, val))
                        nT[ky]=getrow(val, fg.grpfirsti)
                    end
                end
                nT
            end
            merge!(getrow(fd(T, fld), fg.grpfirsti), recurfun(fd(T, otherfd)))
        else
            isuni=map(otherfd) do ofd
                all(isrowsame, eachgrp(fg, T[ofd]))
            end
            if check
                if !all(isuni)
                    t=otherfd[.!isuni]
                    if length(t)>50
                        t=[t[1:50];"..."]
                    end
                    error("Below fields are not intra-group unique: "*join(t, ", "))
                end
                getrow(T, fg.grpfirsti)
            else
                getrow(fd(T, [fld; otherfd[isuni]]), fg.grpfirsti)
            end
        end
    else
        if isnothing(fg)
            rec(T, uniqr(dict2tuple(T, fld); stable=stable)[1])
        else
            rec(T, fg.grpfirsti)
        end
    end
    #30 Jul 2021: Avoid repeatly calculate repeated keys.
    if !isempty(Pr)
        for ((_, fun), ky) in zip(Pr, Pr_fds)
            if isa(ky, AbstractString)
                haskey(T, ky) && error("Field $ky has been assigned to more than one function.")
                T[ky]=grpfun(fun, fg, T0[ky])
            else
                for cky in ky::Union{AbstractVector, Tuple}
                    haskey(T, ky) && error("Field $ky has been assigned to more than one function.")
                    T[cky]=grpfun(fun, fg, T0[cky])
                end
            end
        end
    end
    T
end
tbuniq(T::Dict{Tk,}, fld::AbstractString, Pr::Pair...; kw...) where {Tk<:AbstractString} = tbuniq(T, [fld], Pr...; kw...)
tbuniq(fun::Function, T::Dict{<:AbstractString,}, fld, Pr1::Pair{<:Union{<:AbstractString, <:AbstractVector{<:AbstractString}}, Colon}, Pr::Pair...; kw...)=tbuniq(T, fld, Pr1.first=>fun, Pr...; kw...)

tbuniq_i(T::Dict{Tk,},fld) where {Tk<:AbstractString} =
    rec(T,uniqr(tuple(values(fd_i(T,fld))...))[1])
tbuniq(T::Dict{Tk,}, fld::AbstractString; warg...) where {Tk<:AbstractString} =tbuniq(T, [fld]; warg...)
tbuniq_i(T::Dict{Tk,}, fld::AbstractString; warg...) where {Tk<:AbstractString} =tbuniq_i(T, [fld]; warg...)
export tbuniq, tbuniq_i
#}}

#{{ coo2mx
#[[ coo2mx ]]
# row_id, col_id, Matrix1, ... = coo2mx(row_id::Group, col_id::Group, Vector1, ...;
#                     default=...|nothing, defaults=(d1, d2, ...), row_order=..., col_order=...)
# Change 3-column-table data to a matrix. When default is given, all the inputs will using the same default value. When the default value is nothing, the emptyval() will be used. When default value is ignored, unassinged matrix element will trigger an error.
# See also: pysparse2coo, tbuniq, mxexpand, disent, fullfactors
# Xiong Jieyi, 26 Sep 2019 > 20 Mar 2023

export coo2mx
function coo2mx(R::Group, C::Group, Vs::AbstractVector...; default=undef, defaults::Tuple=(), row_order::Union{Group, Nothing}=nothing, col_order::Union{Group, Nothing}=nothing)
    # siR, iR=uniqr(R)
    # sR=getrow(R, siR)
    # nrow=length(siR)
    # siC, iC=uniqr(C)
    # sC=getrow(C, siC)
    # ncol=length(siC)
    
    if !isnothing(row_order)
        s_row_order, siRO=sortr(row_order)
        # viRO=invperm(siRO)
        sR=unival(s_row_order, sorted=true, sortcheck=false)
        nrow=rownum(sR)
        iR=memberr(R, sR, Bsorted=true, sortcheck=false)
        iR, R, C, Vs=getrow((iR, R, C, Vs), iR.>0)
    else
        iR=nothing
    end
    if !isnothing(col_order)
        s_col_order, siCO=sortr(col_order)
        # viCO=invperm(siCO)
        sC=unival(s_col_order, sorted=true, sortcheck=false)
        ncol=rownum(sC)
        iC=memberr(C, sC, Bsorted=true, sortcheck=false)
        l=iC.>0
        iC, R, C, Vs=getrow((iC, R, C, Vs), l)
        if !isnothing(iR)
            iR=iR[l]
        end
    end
    if isnothing(row_order)
        siR, iR=uniqr(R)
        sR=getrow(R, siR)
        nrow=length(siR)
    end
    if isnothing(col_order)
        siC, iC=uniqr(C)
        sC=getrow(C, siC)
        ncol=length(siC)
    end

    Ms=map(1:length(Vs)) do i
        cdef=isempty(defaults) ? default : defaults[i]
        if isnothing(cdef)
            cdef=emptyval(eltype(Vs[i]))
        end
        if cdef===undef
            Matrix{eltype(Vs[i])}(undef, nrow, ncol)
        else
            fill(cdef, nrow, ncol)
        end
    end
    ML=falses(nrow, ncol)
    for (i, (ir, ic)) in enumerate(zip(iR, iC))
        # if ir==0 || ic==0
        #     continue
        # end
        if ML[ir, ic]
            for mi=1:length(Vs)
                isequal(Ms[mi][ir, ic], Vs[mi][i]) || error("Duplicated row-column with inconsistent values was detected.")
            end
        else
            ML[ir, ic]=true
            for mi=1:length(Vs)
                Ms[mi][ir, ic]=Vs[mi][i]
            end
        end
    end
    if default===undef && isempty(defaults) && !all(ML)
        error(f"$1 out of $2 matrix elements are unassigned. Try set default value."(sum(.!ML), length(ML)))
    end
    if !isnothing(row_order)
        nrow=length(siRO)
        OMs=map(1:length(Vs)) do i
            cdef=isempty(defaults) ? default : defaults[i]
            if isnothing(cdef)
                cdef=emptyval(eltype(Vs[i]))
            end
            if cdef===undef
                Matrix{eltype(Vs[i])}(undef, nrow, ncol)
            else
                fill(cdef, nrow, ncol)
            end
        end
        Rfg=fastgrp(sR, s_row_order, Asorted=true, Bsorted=true, sortcheck=false)
        for gi=1:Rfg.grpnum
            if !isempty((rl=want(Rfg, gi);))
                kri=siRO[gi]
                for mi=1:length(Vs)
                    OMs[mi][kri, :]=Ms[mi][rl, :]
                end
            end
        end
        Ms=OMs
        sR=row_order
    end
    if !isnothing(col_order)
        # ncol=rnum(col_order)
        ncol=length(siCO)
        OMs=map(1:length(Vs)) do i
            cdef=isempty(defaults) ? default : defaults[i]
            if isnothing(cdef)
                cdef=emptyval(eltype(Vs[i]))
            end
            if cdef===undef
                Matrix{eltype(Vs[i])}(undef, nrow, ncol)
            else
                fill(cdef, nrow, ncol)
            end
        end
        Cfg=fastgrp(sC, s_col_order, Asorted=true, Bsorted=true, sortcheck=false)
        for gi=1:Cfg.grpnum
            if !isempty((cl=want(Cfg, gi);))
                kci=siCO[gi]
                for mi=1:length(Vs)
                    OMs[mi][:, kci]=Ms[mi][:, cl]
                end
            end
        end
        Ms=OMs
        sC=col_order
    end
    (sR, sC, Ms...)
end
#}}

#{{ tbx
#[[ tbx ]]
#merged = tbx( A::Dict, B::Dict, key;
#             auniq=false, buniq=false, afull=false, bfull=false, akeep=false, bkeep=false,
#             exp=:AB(Default)|:A|:none, order=:none(Default)|:A|:B,
#             aprefix="Pre-"|function, bprefix=..., rename_conflict=true, must_rename=Set(fields)|fields_vec|"1field", delconflict=false,
#             fillempty=:none(def)|:A|:B|:AB, emptyA=ds(...)|missing, emptyB=... )
#outA, outB = tbx(...; nomerge=true, rmdup=:none(def)|:A|:B, ...)
#key could be "fieldname", ("fieldname_A", "fieldname_B"), (group1, group2), (group1, "fieldname_B"). If the key is omitted, the only shared key will be used.
#Or ... = tbx (A, B;< key=c"key1, key2, ..."|"the_only_key" >|< akey=[c]"...", bkey=[c]"..." >, ...)
#Merge two tables together. When key=(keyA, keyB), the keyB will be not in the output.
#The aprefix and bprefix will be added to the names of conflict fields (rename_conflict=true) or all non-key-fields (rename_conflict=false) in output table. In the rename_conflict mode, fields include must_rename will always be renamed without equality checks.
#delconflict: For the field with conflict values, throw an error (false) or delete this filed (true).
#aprefix/bprefix can also be a function, similar as rename!(function, ...). 
#emptyA and emptyB is Dict for the empty values in A and B. If any fields in emptyA/B are missing, emptyval(...) will be used. Note the field names in emptyA/B should be as the names in the input A/B, rather than the modified names in the output.
#When nomerge=true, tbx() will output two tables with marched rows. The key field(s) in B will be kept, and `delconflict` should only be false.
#rmdup=:A or :B: remove duplicated fields (identical values by isequal()) in either one of two outputs in nomerge=true model.
#See also: vcatr, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dictconc, dictconc!, tbcomp, dictshare
#Xiong Jieyi, 16 May 2014>11 Sep 2014>December 7, 2014>March 26, 2015>Jun 2, 2015>14 Oct 2017>13 Jan 2018>13 Mar 2018>7 Feb 2019>17 Sep 2019>8 Jul 2020>8 Jul 2024

export tbx
function __tbx__(A::Dict, B::Dict, key::Tuple{Group,Group},
             keynameA::AbstractVector{T1}, keynameB::AbstractVector{T2};
             exp::Symbol=:AB, fillempty::Symbol=:none,
             auniq::Bool=false, buniq::Bool=false, afull::Bool=false, bfull::Bool=false, akeep::Bool=false, bkeep::Bool=false, order::Symbol=:none, aprefix::Union{AbstractString, Function}="", bprefix::Union{AbstractString, Function}="", delconflict::Bool=false, rename_conflict::Bool=true, must_rename=Set{String}(), emptyA::Union{Dict, Missing}=ds(), emptyB::Union{Dict, Missing}=ds(), nomerge::Bool=false, rmdup::Symbol=:none) where {T1<:AbstractString,T2<:AbstractString}

    # Input check
    if !isa(must_rename, Set)
        must_rename=if isa(must_rename, AbstractString)
            Set([must_rename])
        else
            Set(must_rename)
        end
    end
    in(fillempty, [:A,:B,:AB,:none]) || error("'fillempty' should only be :A, :B, :AB or :none.")
    in(order, [:A,:B,:none]) || error("'order' should only be :A, :B, :AB or :none.")
    if akeep
        afull=buniq=true
    end
    if bkeep
        bfull=auniq=true
    end
    if fillempty!=:none
        # delconflict && error("delconflict=true is not supported in fillempty model so far.")
        afull && (fillempty==:A || fillempty==:AB) && error("afull/akeep=true conflicts with fillempty=:A/:AB.")
        bfull && (fillempty==:B || fillempty==:AB) && error("bfull/bkeep=true conflicts with fillempty=:B/:AB.")
    end

    # Index match calculation
    if exp==:AB
        (uA,uB)=setxri(key...)
    elseif exp==:A
        (fillempty==:B || fillempty==:AB) && error("When exp=:A, fillempty must be :none or :A.")
        (uA,uB)=setxri(key...;uniqB=true)
    elseif exp==:none
        fillempty==:none || error("When exp=:none, fillempty must be :none.")
        (uA,uB)=intersectri(key...)
    else
        error("exp should only be :AB,:A or :none.")
    end

    auniq && (isuniqr(uB) || error("The A key is not row-unique"))
    buniq && (isuniqr(uA) || error("The B key is not row-unique"))

    afull && (length(uniqr(uA)[1])==rownum(A) || error("A is not full."))
    bfull && (length(uniqr(uB)[1])==rownum(B) || error("B is not full."))

    if order==:A && (fillempty==:none || fillempty==:A)
        (uA,idx)=sortr(uA)
        uB=uB[idx]
    elseif order==:B && (fillempty==:none || fillempty==:B)
        (uB,idx)=sortr(uB)
        uA=uA[idx]
    end
    oA=rec(A,uA)
    oB=rec(B,uB)
    
    #Fill empty
    if fillempty!=:none
        sharedfields=if isempty(aprefix) && isempty(bprefix)
            # setdiff(intersect(keys(A), keys(B)), keyname)
            intersect(setdiff(keys(A), keynameA), setdiff(keys(B), keynameB))
        elseif rename_conflict
            filter(x->isequal(oA[x], oB[x]), intersect(setdiff(keys(A), keynameA), setdiff(keys(B), keynameB)))
        else
            String[]
        end
        if fillempty==:B || fillempty==:AB
            cAi=setdiff(1:rownum(key[1]), uA)
            cA=getrow(A, cAi)
            oA=vcatr(oA, cA)
            
            # cB=emptyrows(B, length(cAi))
            cB=tb()
            for (k, v) in B
                cB[k]=if ismissing(emptyB)
                    missing
                elseif haskey(emptyB, k)
                    reprow(emptyB[k], length(cAi))
                else
                    emptyrows(v, length(cAi))
                end
            end
            
            for x in sharedfields
                cB[x]=cA[x]
            end
            oB=vcatr(oB, cB)

            if order==:A
                idx=sortri([uA; cAi])
                oA=getrow(oA, idx)
                oB=getrow(oB, idx)
            elseif order==:B
                uB=[uB; rownum(key[2])+(1:length(cAi))]
            end
        end
        if fillempty==:A || fillempty==:AB
            cBi=setdiff(1:rownum(key[2]), uB)
            cB=getrow(B, cBi)
            oB=vcatr(oB, cB)
            
            # cA=emptyrows(A, length(cBi))
            cA=tb()
            for (k, v) in A
                cA[k]=if ismissing(emptyA)
                    missing
                elseif haskey(emptyA, k)
                    reprow(emptyA[k], length(cBi))
                else
                    emptyrows(v, length(cBi))
                end
            end
            
            if !isempty(keynameA)
                for (cka, ckb) in zip(keynameA, keynameB)
                    cA[cka]=cB[ckb]
                end
                # cA[keyname[1]]=cB[keyname[2]]
            end
            for x in sharedfields
                cA[x]=cB[x]
            end
            oA=vcatr(oA, cA)

            if order==:B
                idx=sortri([uB; cBi])
                oA=getrow(oA, idx)
                oB=getrow(oB, idx)
            end
        end
    end
    
    #Note that I cannot operate delete! or rename! on A or B directly, since A and B are not copied.
    #Remove B key field
    if !nomerge && !isempty(keynameB)
        for ckb in keynameB
            delete!(oB, ckb)
        end
    end
    # Field rename
    if !(isempty(aprefix) && isempty(bprefix))
        if rename_conflict
            delconflict && error("rename_conflict=true conflicts with delconflict=true")
            confkys=collect(filter(ky->in(ky, must_rename::Set) || !isequal(oA[ky], oB[ky]), intersect(keys(oA), keys(oB))))
            if !isempty(aprefix)
                oA=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, oA, confkys)
            end
            if !isempty(bprefix)
                oB=rename!(isa(bprefix, AbstractString) ? x->bprefix*x : bprefix, oB, confkys)
            end
        else
            if !isempty(aprefix)
                if isempty(keynameA)
                    oA=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, oA)
                else
                    oA=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, oA, exclude=keynameA)
                end
            end
            if !isempty(bprefix)
                oB=rename!(isa(bprefix, AbstractString) ? x->bprefix*x : bprefix, oB)
            end    
        end
    end
    if nomerge
        delconflict && error("nomerge=true conflicts with delconflict=true.")
        if rmdup!=:none
            if rmdup==:B
                filter!(oB) do (k, v)
                    !(in(k, keynameB) || (haskey(oA, k) && isequal(v, oA[k])))
                end
            elseif rmdup==:A
                filter!(oA) do (k, v)
                    !(in(k, keynameA) || (haskey(oB, k) && isequal(v, oB[k])))
                end
            else
                error("Invalid rmdup=$rmdup. Only :A, :B or :none is allowed.")
            end
        end
        (oA, oB)
    else
        rmdup==:none || error("rmdup=$rmdup conflicts with nomerge=true.")
        dictconc!(oA,oB; delconflict=delconflict)
    end
end
tbx(A::Dict, B::Dict, key::Tuple{T1,T2}; P...) where {T1<:AbstractString,T2<:AbstractString} =
    __tbx__(A, B, (A[key[1]], B[key[2]]), [key[1]], [key[2]]; P...)
tbx(A::Dict, B::Dict, key::Tuple{T,Group}; P...) where {T<:AbstractString} =
    __tbx__(A, B, (A[key[1]], key[2]), [key[1]], String[]; P...)
tbx(A::Dict, B::Dict, key::Tuple{Group, T}; P...) where {T<:AbstractString} =
    __tbx__(A, B, (key[1], B[key[2]]), String[], [key[2]]; P...)
tbx(A::Dict, B::Dict, key::Tuple{Group, Group}; P...)=__tbx__(A, B, key, String[], String[]; P...)
tbx(A::Dict, B::Dict, key::AbstractString; P...)=tbx(A, B, (key, key); P...)
function tbx(A::Dict,B::Dict; key=String[], akey::Union{AbstractVector{T1}, T1}=key, bkey::Union{AbstractVector{T2}, T2}=key, wargs...) where {T1<:AbstractString,T2<:AbstractString}
    if isempty(akey)
        sharekey=intersect(keys(A),keys(B))
        length(sharekey)==1 || error("More than one key is shared, or no key is shared. Please assign the linking key.")
        tbx(A,B,first(sharekey); wargs...)
    else
        if isa(akey, AbstractString)
            akey=[akey]
        end
        if isa(bkey, AbstractString)
            bkey=[bkey]
        end
        length(akey)==length(bkey) || error("akey=... and bkey=... should be in the same lengths.")
        __tbx__(A, B, (dict2tuple(A, akey), dict2tuple(B, bkey)), akey, bkey; wargs...)
    end
end
#}}

#{{ tbcomp
#[[ tbcomp ]]
#T|F=tbcomp(Tb_Small, Tb_Big[, c"fields"]; eq=false)
#Test if the combinations of given fields in first table is exist in second table (eq=fasle), or if their uniqued combinations are the same (eq=true). The repeated combinations will be uniqued.
#See also: vcatr, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dictconc, dictconc!
#Xiong Jieyi, January 8, 2015 >Apr 19, 2016

export tbcomp
function tbcomp(A::Dict, B::Dict, fds::AbstractVector; eq::Bool=false)
    @assert(istable(A),"A is not a table.")
    @assert(istable(B),"B is not a table.")
    GA=dict2tuple(A,fds)
    GB=dict2tuple(B,fds)
    @assert(typeof(GA)==typeof(GB),"Type of two tables is inconsistent. A: $(typeof(GA)), B: $(typeof(GB)).")
    if eq
        unival(GA)==unival(GB)
    else
        all(ismemberr(GA,GB))
    end
end
tbcomp(A::Dict, B::Dict; wargs...)=tbcomp(A::Dict, B::Dict, collect(intersect(keys(A), keys(B))); wargs...)
#}}

#{{ readtb writetb autotypetb
#[[ readtb writetb autotypetb! ]]
# table = readtb(filename[.gz]; dm='\t', head=c"data1, N::data2, ...", skipline=0, autotype=false, quotes=false, corner=e.g."N::NO", NA2NaN=false, headincom=false)
# autotypetb!(Table; keepstr=c"field1, ...", NA2NaN=false)
# writetb(filename[.gz], Table, "key1[=>rename1], key2[=>rename1], ..."; dm='\t', with_type=true)
# writetb(filename[.gz], Table, c"key1, key2, ..."; dm='\t', with_type=true) #not parse '=>'.
#Read tsv file, using the first line as table field name (or use the head=c"..." instead). The #-started lines will be regarded as comments and be skipped.
#The elements of first line could be started with below symbol. Otherwise, data will be readed as string:
#N:: Int number *
#n:: Int number, empty will be converted to 0.
#F:: Float number, only NaN will be converted to NaN. *
#f:: Float number, empty and NaN, NA, na, Na, n.a., "NaN", "NA", "na", "Na", "n.a.", "NaN" will be converted to NaN.
#B:: Bool number (lowercase of values should be any of "t,true,f,false,y,yes,n,no,1,0") *
#S:: Symbol *
#T:: String, keeping quote mark (") also. In default, the "..." will be striped.
#C:: Char (Only support one charactor) *
#X:: Ignore this column
#The head could also be < , which means this column will be concaterate with its left column and output a matrix. *
# --Items with * are supported by writetb.
#In the case header is in a shorter length than column number, set the last header as "..." will ignore the right columns, set the last header as <<< will merge all the right columns (equal to '<' for all right columns), Otherwise an error will be throwed.
#autotype: transfer all string fields to BitVector, Int or Float64 once possible ("NaN" is accepted as Float64 but "nan", et al. is not accepted; BitVector only accept TRUE, True, true, FALSE, False, false).
#autotypetb!(): Just seperated the autotype function from readtb, so user can read multiple tables and merge them by vcatr() before autotypetb!().
#NA2NaN: In autotype mode, transfer NA to NaN also. It is useful when handling R's data.
#quotes: When quotes=true, all "..." will be read as a whole even if there is dm between. The quote marks will be stripped.
#corner: When the field name of the first column is missed in the file, corner= can add this field name. Only work when head= is unsigned.
#headincom: Whether the head line is started with #. Note that when headincom=true, the head line must be the last '#...' line just before data line, and it cannot be quoted even under quotes=true.
#with_type: weather write with the type mark, e.g. N::xxx. If with_type is false, 'field < ...' will be replaced to 'field.1 field.2 ...'.
#Support gzip compression. Just use `xxxx.gz' as filename.
#See also: filefun
#Xiong Jieyi, February 22, 2015 >Sep 4, 2015 >Aug 17, 2016 >15 May 2017>31 May 2017>22 Nov 2018>17 Dec 2018>12 Apr 2021>22 Oct 2021>28 Mar 2022

export readtb, writetb, autotypetb!

function __parse_tb_head__(head::Vector{T}) where {T<:AbstractString}
    
    #In the case field name is quited by `"':
    for i=1:length(head)
        if length(head[i])>=2 && head[i][1]=='"' && head[i][end]=='"'
            head[i]=head[i][2:end-1]
        end
    end
    
    if head[end]=="..."
        rightCols=-1
        head=head[1:end-1]
    elseif head[end]=="<<<"
        rightCols=1
        head=head[1:end-1]
    else
        rightCols=0
    end
    dtfun=Function[]
    fdnm=String[]
    colidx=Array{Int}[]
    for (i,cf) in enumerate(head)
        if startswith(cf,"X::")
            continue
        elseif cf=="<"
            push!(colidx[end],i)
            continue
        else
            push!(colidx,[i])
        end
        if startswith(cf,"N::")
            push!(dtfun,int)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"n::")
            push!(dtfun,x::AbstractString->isempty(x) ? 0 : int(x))
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"F::")
            push!(dtfun,float)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"f::")
            push!(dtfun,x::AbstractString->(isempty(x)||in(x,["NA","na","Na","n.a.","\"NA\"","\"na\"","\"Na\"","\"n.a.\"","\"NaN\""])) ? NaN : float(x))
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"B::")
            push!(dtfun,x::AbstractString->begin
                  if lowercase(x) in c"t,true,1,y,yes"
                  true
                  else
                  @assert(lowercase(x) in c"f,false,0,n,no",
                          "Invalid logical symbol $x")
                  false
                  end
                  end)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"C::")
            push!(dtfun,x::AbstractString->begin
                  @assert(length(x)==1,
                          "Char value should only be one charactors.")
                  x[1]
                  end)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"S::")
            push!(dtfun,x::AbstractString->Symbol(x))
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"T::")
            push!(dtfun,x::AbstractString->String(x))
            push!(fdnm,cf[4:end])
        else
            push!(dtfun,x::AbstractString->(length(x)>=2 && x[1]=='"' && x[end]=='"') ? String(x[2:end-1]) : String(x))
            push!(fdnm,cf)
        end
    end
    (dtfun,fdnm,colidx,rightCols)
end

function readtb(fn; dm='\t', head::Vector{Ts}=String[], skipline=0, quotes::Bool=false, autotype::Bool=false, corner::AbstractString="", NA2NaN::Bool=false, headincom::Bool=false) where {Ts<:AbstractString}
    column=Ref(-1)
    rightCols=Ref(0)
    if isempty(head)
        dtfun=Function[]
        fdnm=String[]
        colidx=Array{Int}[]
    else
        isempty(corner) || error("corner=... cannot be used when head=... is assigned.")
        (dtfun,fdnm,colidx,rightCols.x)=__parse_tb_head__(head)
        column.x=length(head)-(rightCols.x!=0)
    end
    strForSure=Set([])
    if isa(fn, AbstractString) && occursin(r"[^\\]\.gz$", fn)
        fn=`zcat $fn`
    end
    lastline=""
    isheadline=isempty(head)
    T=filefun(fn; input_rowno=true, skipline=skipline)do no, ln
        if (length(ln)>=1) && ln[1]=='#'
            # m=match(r"^#\:\:\s+",ln)
            # if m==nothing
            lastline=ln
            return nothing
            # else
            #     ln=ln[length(m.match)+1:end]
            # end
        end
        C=split(ln, dm)
        if quotes
            isdel=falses(length(C))
            isqt=fill(false, length(C)) #The collect() is used to avoid a Julia bug in v1.4.1
            i=1
            while !isnothing((i=findnext(x->length(x)>0 && x[1]=='"', C, i);))
                j=findnext(x->length(x)>0 && x[end]=='"', C, i)
                if isnothing(j)
                    break
                else
                    C[i]=join(C[i:j], dm)[2:end-1]
                    isdel[i+1:j].=true
                    isqt[i]=true
                    i=j+1
                end
            end
            deleteat!(C, isdel)
            deleteat!(isqt, isdel) #Here, if isqt is a BitArray, it will trigger a bug in Julia v1.4.1
            
            for x in findall(isqt)
                push!(strForSure, x)
            end
        end
        if isheadline
            isheadline=false
            if headincom && !isempty(lastline)
                let C=split(lstrip(lstrip(lastline, '#')), dm)
                    if !isempty(corner)
                        C=[corner, C...]
                    end
                    (dtfun,fdnm,colidx,rightCols.x)=__parse_tb_head__(C)
                    column.x=length(C)-(rightCols.x!=0)
                end
            else
                if !isempty(corner)
                    C=[corner, C...]
                end
                (dtfun,fdnm,colidx,rightCols.x)=__parse_tb_head__(C)
                column.x=length(C)-(rightCols.x!=0)
                empty!(strForSure)
                return nothing
            end
        end
        if rightCols.x==0 ? length(C)!=column.x : length(C)<column.x
            error("Invalid column number $(length(C)) (expected $(column.x)) in row $no.")
        end
        T=ds()
        for i=1:length(colidx)
            cidx=colidx[i]
            if i==length(colidx) && rightCols.x==1
                T[fdnm[i]]=map(dtfun[i],C[cidx[1]:end])
            elseif length(cidx)==1
                T[fdnm[i]]=dtfun[i](C[cidx[1]])
            else
                T[fdnm[i]]=map(dtfun[i],C[cidx])
            end
        end
        T
    end
    if autotype
        strForSureD=Set(fdnm[collect(strForSure)])
        for (k, v) in T
            in(k, strForSureD) && continue
            if eltype(v)<:AbstractString
                v=strip.(v)
                if all(x->in(x, c"TRUE, True, true, FALSE, False, false"), v)
                    T[k]=BitVector(parse.(Bool, lowercase.(v)))
                elseif all(x->occursin(r"^\-?\d+$", x), v)
                    T[k]=parse.(Int, v)
                elseif NA2NaN && all(x->occursin(r"^[eE\d\-\+\.]+$|^NaN$|^NA$", x), v)
                    v[v.=="NA"].="NaN"
                    try
                        T[k]=parse.(Float64, v)
                    catch
                    end
                elseif all(x->occursin(r"^[eE\d\-\+\.]+$|^NaN$", x), v)
                    try
                        T[k]=parse.(Float64, v)
                    catch
                    end
                end
            end
        end
    elseif NA2NaN
        error("NA2NaN=true can only be used in autotype mode.")
    end
    T
end

function autotypetb!(T::Dict{<:AbstractString,}; keepstr::Union{<:AbstractVector{<:AbstractString}, Set{<:AbstractString}, Nothing}=nothing, NA2NaN::Bool=false)
    if !isnothing(keepstr)
        strForSureD=if isa(keepstr, Set)
            keepstr
        else
            Set(keepstr)
        end
    end
    for (k, v) in T
        !isnothing(keepstr) && in(k, strForSureD) && continue
        if eltype(v)<:AbstractString
            v=strip.(v)
            if all(x->in(x, c"TRUE, True, true, FALSE, False, false"), v)
                T[k]=BitVector(parse.(Bool, lowercase.(v)))
            elseif all(x->occursin(r"^\-?\d+$", x), v)
                T[k]=parse.(Int, v)
            elseif NA2NaN && all(x->occursin(r"^[eE\d\-\+\.]+$|^NaN$|^NA$", x), v)
                v[v.=="NA"].="NaN"
                try
                    T[k]=parse.(Float64, v)
                catch
                end
            elseif all(x->occursin(r"^[eE\d\-\+\.]+$|^NaN$", x), v)
                try
                    T[k]=parse.(Float64, v)
                catch
                end
            end
        end
    end
    T
end

function writetb(filename::AbstractString, T::Dict{Ts1, Any}, keyls::Union{AbstractString, AbstractVector{Ts2}}=sort(collect(keys(T))); dm='\t', with_type::Bool=true) where {Ts1<:AbstractString, Ts2<:AbstractString}
    istable(T) || error("Input is not a table.")
    if isa(keyls, AbstractString)
        keyls=strip.(split(keyls, ','))
        parseflag=true
    else
        parseflag=false
    end
    cols=map(enumerate(keyls)) do (i, ky)
        if parseflag
            t=split(ky, "=>")
            if length(t)==2
                ky, wky=strip.(t)
                keyls[i]=ky
            else
                wky=ky
            end
        else
            wky=ky
        end
        occursin(dm, wky) && error("Delim char (dm=$dm) is not allowed in key $wky.")
        vv=T[ky]
        if with_type
            vt=eltype(vv)
            pf=if vt<:Bool
                "B::"
            elseif vt<:Integer
                "N::"
            elseif vt<:AbstractFloat
                "F::"
            elseif vt<:Symbol
                "S::"
            elseif vt<:Char
                "C::"
            elseif vt<:AbstractString
                ""
            else
                error("Unsupported type $vt in field $ky.")
            end
            ["$pf$wky", fill("<", size(vv, 2)-1)...]
        else
            if size(vv, 2)>1
                f"$wky.$1".(1:size(vv, 2))
            else
                [wky]
            end
        end
    end
    isgz=occursin(r"[^\\]\.gz$", filename)
    if isgz
        filename=replace(filename, r"\.gz$"=>"")
        if isfile(filename)
            error("File $filename exists, which is conflicted with the temp file for gzip.")
        end
    end
    writedlm(filename, asmx(Pair.(cols, dict2tuple(T, keyls))...), dm)
    if isgz
        run(`gzip -f $filename`)
    end
end

#}}

#{{ rename!
#[[ rename! ]]
#Dict=rename!(Dict, Key1=>NewKey1, Key2=>NewKey2, ...) or rename!(Dict, Key_vec, NewKey_vec)
#... =rename!(Dict, Regex=>Replace_string[, c"field1,field2,..."|Regex];
#            exclude=Vector|Tuple|Set|AbstractString|Regex) #use replace() for fields.
#... =rename!(function, Dict[, key_vec|Regex]; exclude=...) #replace key to function(key).
#Change the key in Dict to the new name, and return the changed Dict.
#e.g. rename!(D, regex1=>s"...", regex2; exclude=regex3), the keys of D matches regex2 and not match regex3 will be renamed as replace(key, regex1 => s"...").
#For compat reason, rename!(Dict, "key"|Regex, newkey, ...) is still supported but it is not recommended.
#See also: fd, fd_i, delete!
#Xiong Jieyi, 10 Sep 2014>December 12, 2014>March 25, 2015>Jun 2, 2015>Nov 12, 2015>29 Apr 2019

export rename!
function rename!(D::Dict{Tk,},A,B) where {Tk<:AbstractString}
    if A!=B
        D[B]=D[A]
        delete!(D,A)
    else
        D
    end
end
function rename!(D::Dict{Tk,},A::T,B::T) where {Tk<:AbstractString, T<:Union{Vector,Tuple}}
    for i=1:length(A)
        rename!(D,A[i],B[i])
    end
    D
end
function rename!(D::Dict{Tk,},A::Regex,B::AbstractString,fds::Union{AbstractVector,Tuple}=collect(keys(D)); exclude::Union{AbstractVector,AbstractSet,Tuple,AbstractString,Regex}=()) where {Tk<:AbstractString}
    if !isempty(exclude)
        if isa(exclude,AbstractString)
            fds=fds[fds.!=exclude]
        elseif isa(exclude, Regex)
            fds=filter(x->!occursin(exclude, x), fds)
        else
            fds=setdiff(fds,exclude)
        end
    end
    nfds=map(x->replace(x, A => B),fds)
    l=fds.!=nfds
    fds=fds[l]
    nfds=nfds[l]
    isuniqr(nfds) || error("Duplicated fieldname occurred after replacement.")
    rename!(D,fds,nfds)
end
function rename!(F::Function, D::Dict{Tk,}, fds::Union{AbstractVector,Tuple}=collect(keys(D)); exclude::Union{AbstractVector,AbstractSet,Tuple,AbstractString,Regex}=())  where {Tk<:AbstractString}
    if !isempty(exclude)
        if isa(exclude,AbstractString)
            fds=fds[fds.!=exclude]
        elseif isa(exclude, Regex)
            fds=filter(x->!occursin(exclude, x), fds)
        else
            fds=setdiff(fds,exclude)
        end
    end
    nfds=map(F, fds)
    l=fds.!=nfds
    fds=fds[l]
    nfds=nfds[l]
    isuniqr(nfds) || error("Duplicated fieldname occurred after replacement.")
    rename!(D,fds,nfds)
end
rename!(F::Function, D::Dict{Tk,}, ft::Regex; kw...) where {Tk<:AbstractString}=rename!(F, D, collect(filter(ft, keys(D))); kw...)
rename!(D::Dict{Tk,}, P::Pair, ft::Regex; kw...) where {Tk<:AbstractString} = rename!(D, P, collect(filter(ft, keys(D))); kw...)
rename!(D::Dict{Tk,}, P::Pair{Regex, T}, arg...; kw...) where {Tk<:AbstractString, T<:Union{AbstractString, Function}} = rename!(D, P.first, P.second, arg...; kw...)
function rename!(D::Dict{Tk,}, pr1::Pair, pr::Pair...) where {Tk<:AbstractString}
    for (a,b) in (pr1, pr...)
        rename!(D,a,b)
    end
    D
end
#}}

#{{ tb2DataFrame, DataFrame2tb
#[[ tb2DataFrame ]]
# dataframe = tb2DataFrame( Table[, fields ])
#Convert table to DataFrame
#If field "F" is a matrix, its column will be splitted into :F_1, :F_2, ....
#See also: readtb, sk, sk_s, pandasDF, DataFrame2tb
#Xiong Jieyi, Jul 10, 2015 > 2 May 2017

export tb2DataFrame
#<v0.6# function tb2DataFrame{Tk<:AbstractString}(T::Dict{Tk,}, fields=keys(T))
function tb2DataFrame(T::Dict{Tk,}, fields=keys(T)) where {Tk<:AbstractString}
    D=importpkg(:DataFrames, preloaded=true).DataFrame()
    function foo(k,v::AbstractMatrix)
        for i=1:size(v,2)
            D[Symbol("$(k)_$i")]=v[:,i]
        end
    end
    function foo(k,v)
        D[Symbol(k)]=v
    end
    for k in fields
        foo(k,T[k])
    end
    D
end

#[[ DataFrame2tb ]]
# dataframe = DataFrame2tb( DataFrame )
#Convert DataFrame to table.
#NA in data will be converted to empty value, and a warning will occur.
#See also: readtb, sk, sk_s, pandasDF, tb2DataFrame
#Xiong Jieyi, Oct 6, 2015

export DataFrame2tb
function DataFrame2tb(D::DataFrame)
    T=tb()
    # isna(x)=invokelatest(importpkg(:DataFrames).isna, x)
    # @pkgfun(DataFrames,DataFrame,isna,getindex,names,collect)
    # DF=importpkg(:DataFrames)
    for k in names(D)
        F=D[k]
        l=ismissing.(F)
        if any(l)
            emp=emptyval(F[findfirst(.!l)])
            V=fill(emp,length(F))
            V[.!l]=F[.!l]
            @warn(f"$1/$2 are NA in field $3. Replaced by empty value `$4'($5)."(sum(l),length(l),k,emp,typeof(emp)))
        else
            # V=[F...]
            V=collect(F)
        end
        T[string(k)]=V
    end
    T
end
#}}

#{{ pandasDF pandas2tb
#[[ pandasDF pandas2tb ]]
# dataframe = pandasDF( Table )
# dataframe = pandasDF(tb_style_params...) # =pandasDF(tb(...))
# dataframe = pandasDF(data; index=..., colomns=...) # =pandas.DaraFrame(...)
# table = pandas2tb( pandas_dataframe; index_name="..."|nothing, key=c"...") # Requiring PyCall.
#Convert table to python pandas DataFrame or vise versa.
#For pandasDF, if field "F" is a matrix, its column will be splitted into "F:1", "F:2", ....
#For pandas2tb, index in pandas will be convert to index_name, or be discarded if index_name=nothing. When key=stringVec is given, only the assigned keys and the index will be converted.
#See also: readtb, sk, sk_s, tb, ds, tb2DataFrame
#Xiong Jieyi, Aug 16, 2015 > Jan 12, 2016 >26 Jun 2019>16 Dec 2020>16 Dec 2021

export pandasDF, pandas2tb
#<v0.6# function pandasDF{Tk<:AbstractString}(T::Dict{Tk,})
function pandasDF(T::Dict{Tk,}) where {Tk<:AbstractString}
    D=tb()
    function foo(k,v::AbstractMatrix)
        for i=1:size(v,2)
            D["$k:$i"]=v[:,i]
        end
    end
    function foo(k,v)
        D[k]=v
    end
    for (k,v) in T
        foo(k,v)
    end
    importpy(:pandas).DataFrame(D)
end
pandasDF(args...;wargs...)=pandasDF(tb(args...;wargs...))
pandasDF(D::AbstractArray;wargs...)=importpy(:pandas).DataFrame(D;wargs...)

function pandas2tb(pd; index_name::Union{AbstractString, Nothing}=something(pd.index.name, "index"), key::Union{AbstractVector{<:AbstractString}, Nothing}=nothing)
    T=tb()
    # pyclass=pd.__module__
    hasattr=importpkg(:PyCall, preloaded=true).pybuiltin("hasattr")
    if hasattr(pd, "__getitem__") && hasattr(pd, "index")
        for ky in (isnothing(key) ? pd : key)
            T[ky]=pd.__getitem__(ky).to_list()
        end
        if !isnothing(index_name)
            T[index_name]=pd.index.to_list()
        end
        T
    elseif pd.__module__=="pandas.core.series"
        pandas2tb(pd.to_frame(); index_name=index_name, key=key)
    else
        @error(f"Unsupported python class $1"(pd.__module__))
    end
end
#}}

#{{ asmx
#[[ asmx ]]
# Matrix = asmx(Matrix[, ...]; row=row_header_vec, col=col_header_vec, rowexp=N)
# Or ... = asmx(c"colname1, colname2, ..."=>Matrix|Tuple_of_columns, "header3"=>Vector,
#               "colname"=>Matrix, #The heads of 2:end colnums are '<'. Not support Tuple.
#               c"key1, key2=>colname1;colname2;..., '...', ..."=>Table, #You can literally use three dots '...' to represent all rest fields. Can be only used once.
#               Table; #Column is in alphabetic order. Tabel cannot be nested or contains Tuple.
#               key_dm="=>", col_dm=';', otherkeys_mk="...") #You can also change these markers.
#Composite a matrix, and ddd the row header and column header to a matrix.
#In the pairs input mode, the length of header and column number will be checked in each pair. It could avoid programmer's mistake.
#If rowexp>0, the single-row input will be expand to N, and other input should has N row.
#See also: writedlm, tocb, coo2mx
#Xiong Jieyi, Jul 6, 2017 > 28 Oct 2017 > 21 Oct 2019 > 18 Jun 2020 > 16 Feb 2023

export asmx
function asmx(M::Union{AbstractVector, AbstractMatrix}; row::Union{AbstractVector,Void}=nothing, col::Union{AbstractVector,Void}=nothing,corner="")
    if row==nothing && col==nothing
        M
    elseif row==nothing
        size(M,2)==length(col) || error(f"Column header vector length should be $1 but it is $2."(size(M,2), length(col)))
        [trsp(col); M]
    elseif col==nothing
        size(M,1)==length(row) || error(f"Row header vector length should be $1 but it is $2."(size(M,1), length(row)))
        [row M]
    else
        size(M,1)==length(row) || error(f"Row header vector length should be $1 but it is $2."(size(M,1), length(row)))
        size(M,2)==length(col) || error(f"Column header vector length should be $1 but it is $2."(size(M,2), length(col)))
        [trsp(Any[corner; col]); Any[row M]]
    end
end

function asmx(M1, Ms...; rowexp::Integer=0, wargs...)
    C=Any[M1, Ms...]
    if rowexp>0
        for (i, cM) in enumerate(C)
            rn=rownum(cM)
            if rn==1
                C[i]=reprow(cM, rowexp)
            else
                rn==rowexp || error("The $(i+1)-th input ($rn rows) is neither single row nor match the row number $rowexp.")
            end
        end
    end
    asmx(hcat(Array{Any}.(C)...); wargs...)
end

# function asmx(Ps::Pair... ; wargs...)
#     for (p1, p2) in Ps
#         if isa(p1, AbstractArray) && !isa(p2, Dict)
#             length(p1)==(isa(p2, Tuple) ? sum(size.(p2, 2)) : size(p2,2)) || error(f"Header \"$1\"...\"$2\" has $3 elements, but matrix has $4 columns."(p1[1],p1[end],length(p1),size(p2,2)))
#         else
#             isa(p2, AbstractVector) || isa(p2, Dict) || (size(p2,2)==1) || error("Single header \"$p1\" does not match a vector column.")
#         end
#     end
#     M=hcat(map(x->isa(x.second, Tuple) ? hcat(Array{Any}.(x.second)...) :
#                   isa(x.second, Dict) ? hcat(Array{Any}.(dict2tuple(x.second, x.first))...) : x.second, Ps)...)
#     header=Any[map(x->x.first, Ps)...;]
#     asmx(M; col=header, wargs...)
# end

# function asmx(Ts::Union{Pair, Dict}...; wargs...)
#     Ps=map(Ts) do T
#         if isa(T, Dict)
#             ky=sort(collect(keys(T)))
#             ky=>dict2tuple(T, ky)
#         else
#             T
#         end::Pair
#     end
#     asmx(Ps...; wargs...)
# end

function asmx(Ts::Union{Pair, Dict}...; key_dm="=>", col_dm=';', otherkeys_mk::AbstractString="...", wargs...)
    function foo(K::AbstractString, X::AbstractVector)
        ([K], Any[X])
    end
    function foo(K::AbstractString, X::AbstractMatrix)
        coln=size(X, 2)
        ([K; fill("<", coln-1)], Any[X])
    end
    function foo(K::AbstractVector{<:AbstractString}, X::Union{AbstractMatrix, AbstractVector})
        @assert(length(K)==size(X, 2), "Column name length is not equal to matrix width.")
        (K, Any[X])
    end
    function foo(K::AbstractVector{<:AbstractString}, X::Tuple)
        @assert(length(K)==length(X), "Column name length is not equal to tuple length. Note c\"key=>col1;col2;...\"=>table is only allowed for a table.")
        (K, Any[X...])
    end
    function foo(K::AbstractVector{<:AbstractString}, X::Dict{<:AbstractString,})
        oK=AbstractString[]
        oX=Any[]
        for ck in K
            t=if ck==otherkeys_mk
                l=K.==otherkeys_mk
                @assert(count(l)==1, "Other-keys-marker '$otherkeys_mk' was used more than once.")
                rK=K[.!l]
                foo(sort(collect(setdiff(keys(X), map(x->strip(split(x, key_dm, limit=2)[1]), rK)))), X)
            elseif contains(ck, key_dm)
                (kk, fk)=strip.(split(ck, key_dm, limit=2))
                if (isa(X[kk], AbstractMatrix) || isa(X[kk], Tuple)) && contains(fk, col_dm)
                    foo(strip.(split(fk, col_dm)), X[kk])
                else
                    foo(fk, X[kk])
                end
            else
                foo(ck, X[ck])
            end
            append!(oK, t[1])
            append!(oX, t[2])
        end
        (oK, oX)
    end
    function foo(X::Dict{<:AbstractString,})
        foo(sort(collect(keys(X))), X)
    end
    function foo(X::Pair)
        foo(X.first, X.second)
    end
    
    colnm=AbstractString[]
    dat=Any[]
    for cT in Ts
        t=foo(cT)
        append!(colnm, t[1])
        append!(dat, t[2])
    end
    asmx(dat...; col=colnm, wargs...)
end
#}}

#{{ sortri_as
#[[ sortri_as ]]
# Index = sortri_as(Group => :i|:d|:s|order_list, ...)
# Sorting each factor in given order. The order could be a list or
# :i -- increase order
# :d -- decrease order
# :s -- stable order
#Note that the speed of this function is slower than sortri((G1, G2, ...)).
#See also: sortri
#Xiong Jieyi, 28 Jan 2020

export sortri_as
function sortri_as(Ps::Pair...)
    C=map(Ps::Tuple) do (X, G)
        if G==:i
            uniqr(X)[2]
        elseif G==:d
            -uniqr(X)[2]
        elseif G==:s
            uniqr(X, stable=true)[2]
        else
            t=memberr(X, G::Group)
            if any(isequal(0), t)
                error("Some factors are missed in order list.")
            end
            t
        end
    end
    sortri(C)
end
#}}
