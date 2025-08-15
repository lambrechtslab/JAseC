#{{ quantilels
#[[ quantilels ]]
#     quantilels(Data1[, Data2, ...]; q=0:0.1:1)
# Or: quantilels(Data; grp=Group|fastgrp, q=...)
#List the quantile and value on screen.
#So far this function output a Float or Any Matrix, but it is not guarantied in the future version.
#See also: ls, showprop, twofactormx
#Xiong Jieyi, 5 Sep 2014 > 25 May 2018 > 1 Nov 2019

function quantilels(Xs::AbstractVector...;q::AbstractVector{T}=Float64[], grp::Union{Group, fastgrp, Nothing}=nothing) where {T<:AbstractFloat}
    if !isnothing(grp)
        vX, vgrp=if isa(grp, fastgrp)
            (grpvec(grp, onlyone(Xs)), grp.grpid)
        else
            grpvec(grp, onlyone(Xs))
        end
        out=quantilels(vX...; q=q)
        return ["q" vgrp'; out]
    end
    if isempty(q)
        N=maximum(map(rownum,Xs))
        if N>199
            q=[0.0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,1.0]
        elseif N>49
            q=0.0:0.1:1
        else
            q=[0.0,0.25,0.5,0.75,1.0]
        end
    end
    O=map(Xs,1:length(Xs)) do X, i
        l=isnan.(X)
        if any(l)
            @warn(f"$1 ($2%) NaNs in data #$i are ignored."(sum(l),mean(l)*100))
            X=X[.!l]
        end
        X
    end
    [q map(x->quantile(vec(x),q),O)...]
end
quantilels(X::AbstractMatrix;args...)=quantilels(splitdim(X,1)...;args...)
quantilels(Xs::AbstractVector{Any}...;wargs...)=quantilels(map(x->retype(x,Real),Xs)...;wargs...) #For the convience of julia -e 'quantilels(SIC[:,1])'.
export quantilels
#}}

#{{ twofactormx
#[[ twofactormx ]]
# AnyMx = twofactormx(factorA[ => orderA ], factorB[ => orderB ]; showsum=false, vb=0|1)
#Or FreqMx, row, col = twofactormx(...; out1mx=false, ...)
#Calculated frequenes for each factor combinations, and output the result as a matrix. Input can be groups, fastgrp objects or pairs. factorA and factorB will be as row and columns respectively.
#If showsum=true, the last row and column in the outputed matrix is the sum values. Note that when out1mx=false, the outputed row and col will not change no matter showsum=true or false.
#If vb=1, function will print the Fisher test result on screen, but not change the outputs.
#See also: freq
#Xiong Jieyi, 8 Nov 2017 > 9 Jun 2021

export twofactormx
function twofactormx(A::fastgrp, B::fastgrp; out1mx::Bool=true, showsum::Bool=false, vb::Int=0)
    A.itemnum == B.itemnum || error("The two factors are not in the same lengths.")
    M=fill(0, A.grpnum, B.grpnum)
    for ia=1:A.grpnum
        la=want(A,ia)
        for ib=1:B.grpnum
            lb=want(B,ib)
            M[ia, ib]=length(intersect(la,lb))
        end
    end
    if vb>0
        if size(M)==(2, 2)
            fishertest(M, vb=vb)
        else
            @warn("vb=1 only work for 2x2 table.")
        end
    end
    if out1mx
        _stringrow(R::Tuple)=join(_stringrow.(R), "|")
        _stringrow(R::Matrix)=join(_stringrow.(R), ",")
        _stringrow(R)=string(R)
        Aid=isa(A.grpid, Union{Tuple, AbstractMatrix}) ? map(_stringrow, eachr(A.grpid)) : A.grpid
        Bid=isa(B.grpid, Union{Tuple, AbstractMatrix}) ? map(_stringrow, eachr(B.grpid)) : B.grpid
    end
    if showsum
        M=[M sum(M,2);sum(M,1) sum(M)]
        if out1mx
            asmx(M, row=[Aid; "SUM"], col=[Bid; "SUM"])
        else
            (M, A.grpid, B.grpid)
        end
    else
        if out1mx
            asmx(M, row=Aid, col=Bid)
        else
            (M, A.grpid, B.grpid)
        end
    end
end
twofactormx(A, args...; wargs...)=twofactormx(fastgrp(A), args...; wargs...)
twofactormx(A::fastgrp, B; wargs...)=twofactormx(A, fastgrp(B); wargs...)
twofactormx(A::Pair, args...; wargs...)=twofactormx(fastgrp(A.first, A.second), args...; wargs...)
twofactormx(A::fastgrp, B::Pair; wargs...)=twofactormx(A, fastgrp(B.first, B.second); wargs...)
#}}

#{{ rankgrp

#[[ rankgrp ]]
# GroupID=rankgrp(V,N)
#Group by the rank of V based on drawer rule. Return a int vector of group NO.
#Note that considering the repeat elements in the vector, different groups could be have same values.
#See also: drawer, balancemask
#Xiong Jieyi, October 22, 2014

#<v0.6# function rankgrp{T<:Real}(V::AbstractVector{T},N::Real)
function rankgrp(V::AbstractVector{T},N::Real) where {T<:Real}
    gn=drawer(length(V),N)
    gl=vcat([fill(i,gn[i]) for i=1:N]...)
    return gl[rvorder(sortri(V))]
end
export rankgrp

#}}

#{{ kmean, kmeancost
#[[ kmean kmeancost ]]
# using Clustering
# ClusterNO=kmean(X, k; kmeans_param...)
# totalcost=kmeancost(X, k|k_vec; kmeans_param...)
#Calculate cluster of X' using kmeans, and return the assignments(group NO) or totalcost.
#Note that X are clusted by row, which is not like the Clustering.kmeans.
#See also: kmeancostplot, grpboxplot, fastcluster, fcluster, dendrogram.
#Xiong Jieyi, 3 Oct 2014 > 13 May 2019

# function needClustering()
#     isdefined(:Clustering) || require("Clustering")
#     nothing
# end
kmean(X::Matrix,k::Real;wargs...)=invokelatest(importpkg(:Clustering).kmeans, X', k; wargs...).assignments
kmeancost(X::Matrix,k::Real;wargs...)=invokelatest(importpkg(:Clustering).kmeans, X', k; wargs...).totalcost
function kmeancost(X::Matrix,ks::AbstractVector{T};wargs...) where {T<:Real} 
    @pkgfun(Clustering, kmeans)
    map(k->kmeans(X', k; wargs...).totalcost, ks)
end
export kmean, kmeancost
#}}

#{{ fastcluster,fcluster,dendrogram,leaves
#[[ fastcluster linkage fcluster dendrogram leaves ]]
# Z=fastcluster(Matrix; method="ward", metric="euclidean") #Make tree by data
# Z=linkage(DistancesMatrix|DistancesVector; method="ward") #Make tree by distance matrix.
# cluster_NO=fcluster(Z, cutoff; criterion="distance(D)|maxclust", ...)   #Cut tree
# dendrogram(Z; labels=None, color_threshold=height_value, ...)     #Draw tree(Only used in Qt so far.)
# leaves_index = leaves(Z)  #The return corresponds to the observation vector index as it appears in the tree from left to right. e.g. leaves_list(Z).+1 in Python.
#Use scipy.cluster.hierarchy and fastcluster to create a hierarchy tree.
# method could be single, complete, average, mcquitty, centroid, median, and ward.
# metric could be euclidean, maximum, manhattan, canberra, binary, minkowski, correlation, and spearman. The correlation distance is actually calculated by scipy.spatial.distance, and the spearman distance calculated by 1-corspearman(Matrix'). These two metrics are not originally supported by fastcluster.
# criterion in fcluster could be inconsistent, distance(default), maxclust, moncrit and maxclust_monocrit. For the detail, see: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
# fastcluster is a stable python package (https://pypi.org/project/fastcluster/). It needs to be preinstalled in PYTHONPATH, like  PYTHONPATH=/home/jxiong/programs/fastcluster/fastcluster-1.1.23/lib/python2.7/site-packages:$PATYONPATH .
#See also: kmean, kmeancost, distmx2vec
#Xiong Jieyi, 3 Oct 2014>October 31, 2014>Aug 14, 2015>11 Jan 2019

function fastcluster(X::AbstractMatrix;
                     method::Union{String,Symbol}="ward",
                     metric::Union{String,Symbol}="euclidean", wargs...)
    @pypkgfun(fastcluster, linkage)
    if metric=="correlation"
        @pkgfun(Distances,pairwise,CorrDist)
        D=distmx2vec(pairwise(CorrDist(),X'))
        linkage(D; method=method, wargs...)
    elseif metric=="spearman"
        D=max(distmx2vec(1-corspearman(X')),0)
        linkage(D; method=method, wargs...)
    else
        linkage(X; method=method, metric=metric, wargs...)
    end
end
linkage(X::AbstractMatrix; method::Union{String,Symbol}="ward", wargs...)=@il importpy("fastcluster").linkage(distmx2vec(X); method=method, wargs...)
linkage(X::AbstractVector; method::Union{String,Symbol}="ward", wargs...)=@il importpy("fastcluster").linkage(X; method=method, wargs...)
fcluster(Z::Matrix,t::Real; criterion::Union{String,Symbol}="distance",wargs...)=@il importpy("scipy.cluster.hierarchy").fcluster(Z, t; criterion=criterion, wargs...)
dendrogram(args...;wargs...)=@il importpy("scipy.cluster.hierarchy").dendrogram(args...; wargs...)
leaves(Z::Matrix)=@il(importpy("scipy.cluster.hierarchy").leaves_list(Z)).+1
export fastcluster, linkage, fcluster, dendrogram, leaves
#}}

#{{ distmx2vec
#[[ distmx2vec distvec2mx ]]
# distance_vector = dismx2vec(distance_matrix)
# distance_matrix = disvec2mx(distance_vector)
#Convert between distance matrix and consense distance vector.
#The distances are arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1)).
#See also: fastcluster, distvec
#Xiong Jieyi, 3 Oct 2014 > 3 Jun 2019

export distmx2vec, distvec2mx
function distmx2vec(M::AbstractMatrix{T}) where {T<:Real}
    n=size(M, 1)
    size(M, 2)==n || error("Inputed matrix is not a squire.")
    V=T[]
    for i=1:n
        for j=i+1:n
            M[i, j]==M[j, i] || error("Inputed matrix is not symmatric.")
            push!(V, M[j, i])
        end
    end
    V
end
function distvec2mx(V::AbstractVector{T}) where {T<:Real}
    t=sqrt(8*length(V)+1)/2+0.5
    if !isinteger(t)
        error("Invalid distance vector length.")
    end
    n=Int(t)
    
    M=fill(zero(T), n, n)
    p=1
    for i=1:n
        for j=i+1:n
            M[j, i]=V[p]
            M[i, j]=V[p]
            p+=1
        end
    end
    M
end
#}}

#{{ distancevec
#[[ distancevec ]]
# vector = distancevec(fun(row_i, col_i)::Float64, n::Int)
# Calculate consense distance vector. Note that the input of given function is a pair of scaler indice.
# The distances are arranged in the order (2,1), (3,1), ..., (n,1), (3,2), ..., (n,2), ..., (n,n–1)).
# See also: fastcluster, distmx2vec, distvec2mx
# Xiong Jieyi, 3 Jun 2019

export distancevec
function distancevec(fun::Function, N::Int)
    V=Float64[]
    for i=1:N
        for j=i+1:N
            push!(V, fun(j,i))
        end
    end
    V
end
#}}

#{{ louvain_cluster
#[[ louvain_cluster ]]
# clusterNO = louvain_cluster(adjecencyVec|Mx)
# Use Louvain community detection method to cluster. Note input is adjecencies (weights of edges) rather than distances.
# It used python package https://github.com/taynaud/python-louvain
# Outputed clusterNO is start from 1, which is different from the original python function.
# See also: distancevec, distmx2vec, fcluster
# Jieyi Xiong, 11 Dec 2020

export louvain_cluster
function louvain_cluster(adjmx::AbstractMatrix{T}; _no_input_check_::Bool=false) where {T<:Real}
    (_no_input_check_ || isapprox(adjmx, adjmx')) || error("Adjecency matrix is not symatric.")
    nx=importpy("networkx")
    @il G=nx.Graph()
    for i=1:size(adjmx, 1), j=i+1:size(adjmx, 1)
        if adjmx[i, j]>0
            @il G.add_edge(i, j, weight=adjmx[i, j])
        end
    end
    lov=importpy("community")
    @il pt=lov.best_partition(G)
    clusterNO=fill(0, size(adjmx, 1))
    for (k, v) in pt
        clusterNO[k]=v+1
    end
    @assert all(clusterNO.>0)
    clusterNO
end
louvain_cluster(adjvec::AbstractVector{T}) where {T<:Real} = louvain_cluster(distvec2mx(adjvec); _no_input_check_=false)
#}}

#{{ corspearman, nrank
#[[ corspearman ]]
# pho=corspearman(A[, B])
#Calculate Spearman correlation. Input could be vector or matrix. If input is matrix, function calculated column-wisely.
#See also: cor
#StatsBase Pkg, 3 Oct 2014

@pkgfun(StatsBase, corspearman)
export corspearman

#[[ nrank ]]
#rand=nrank(x[, dim]; rev=false, norm=false, ignore_nan=false)
#Compute tied ranking (also known as fractional ranking or 1 2.5 2.5 4 ranking) for x, using StatsBase Pkg.
# rev=true : output the reversed rank (the maximum is 1 and the minimum is N).
# ignore_nan=true : Ignore NaNs in calculation but refill NaNs in the rank.
# norm=true : normalize ranks to (0, 1), which was calculated by (rank-1)/(length(non_nan_value)-1).
#See also: rankgrp
#Xiong Jieyi, November 2, 2014 >Mar 3, 2016>May 9, 2016

function nrank(X::AbstractVector; rev::Bool=false, ignore_nan::Bool=false, norm::Bool=false)
    # sb=importpkg(:StatsBase)
    @pkgfun(StatsBase, tiedrank)
    l=isnan.(X)
    if !ignore_nan && any(l)
        error(f"$1 of $2 data are NaNs. They are also NaNs in the ranks."(sum(l),length(l)))
    end
    oR=tiedrank(X[.!l])
    R=fill(NaN, size(X))
    R[.!l]=if norm
        pR=(oR.-1) ./ (length(oR)-1)
        rev ? (1.0 .- pR) : pR
    else
        if rev
            length(oR).-oR.+1
        else
            oR
        end
    end
    R
end
function nrank(X::AbstractMatrix,dim::Real;wargs...)
    if dim==1
        trsp(rowfun(x->nrank(x;wargs...),X',as_vec=true))
    elseif dim==2
        rowfun(x->nrank(x;wargs...),X,as_vec=true)
    else
        error("dim only support 1 and 2 so far.")
    end
end
export nrank
#}}

#{{ HypothesisTest [ obsoleted ]
#pval=pvalue(fisherexacttest(a,b,c,d))
#pval=pvalue(fisherexacttest([a,b,c,d]))
#pval=pvalue(fisherexacttest(BitVector1,BitVector2))
#Calculate two-tailed Fisher exact test.
#See also: pvalue, fishertest
#Xiong Jieyi, 3 Oct 2014

# @pkgfun(HypothesisTests,FisherExactTest=>fisherexacttest,pvalue)
# function fisherexacttest(M::Matrix)
#     @assert(size(M)==(2,2),"FisherExactTest should be 2x2 matrix.")
#     fisherexacttest(M[1,1],M[1,2],M[2,1],M[2,2])
# end
# fisherexacttest(A::BitVector,B::BitVector)=fisherexacttest(sum(A),sum(B),sum(.!A),sum(.!B))
# export fisherexacttest, pvalue
#}}

#{{ glm12test

#[[ glm12test ]]
# pval=glm12test(Y,X;family="",test="F|Chisq",outall=false)
# Test if data fit linear or quadrotic model, using lm (when family is omitted) or glm. Return the minimal F-test p-value of these two models, or the both two p-values and the p-value between two models if outall is ture.
# family could be:
# binomial(link = "logit")
# gaussian(link = "identity")
# Gamma(link = "inverse")
# inverse.gaussian(link = "1/mu^2")
# poisson(link = "log")
# quasi(link = "identity", variance = "constant")
# quasibinomial(link = "logit")
# quasipoisson(link = "log")
# nb --Need run eR("library(MASS)") in advance.
#See also: callr, rdata, relem, rshow, callr("p.adjust",pval,method="BH")
#Xiong Jieyi, October 15, 2014>October 23, 2014

function glm12test(Y::Array,X::Vector;family::AbstractString="",test::AbstractString="F",outall::Bool=false)
    @assert(size(Y,1)==length(X),"Y and X have different length.")
    # Rif.isinitialized()||Rif.initr()
    # import needR
    @pkgfun(myR, callr, eR)
    if isa(Y,Matrix) && size(Y,2)==2 && size(Y,1)>1
        df=callr("data.frame",Y1=Y[:,1],Y2=Y[:,2],X=X)
        sY="cbind(Y1,Y2)"
    else
        df=callr("data.frame",Y=Y,X=X)
        sY="Y"
    end
    if isempty(family)
        md0=callr("lm",eR("$(sY)~1"),data=df)
        # md1=callr("lm",R("$(sY)~X"),data=df)
        # md2=callr("lm",R("$(sY)~I(X^2)+X+1"),data=df)
    elseif family=="nb"
        md0=callr("glm.nb",eR("$(sY)~1"),data=df)
        # md1=callr("glm.nb",R("$(sY)~X"),data=df)
        # md2=callr("glm.nb",R("$(sY)~I(X^2)+X+1"),data=df)
    else
        md0=callr("glm",eR("$(sY)~1"),data=df,family=eR(family))
        # md1=callr("glm",R("$(sY)~X"),data=df,family=R(family))
        # md2=callr("glm",R("$(sY)~I(X^2)+X+1"),data=df,family=R(family))
    end
    md1=callr("update",md0,eR(".~.+X"))
    md2=callr("update",md1,eR(".~.+I(X^2)"))

    pval1=callr("anova",md1,md0,test=test)|>x->relem(x,length(x))|>x->[x...]|>x->x[~isnan(x)]|>x->(@assert(length(x)==1);x[1])
    pval2=callr("anova",md2,md0,test=test)|>x->relem(x,length(x))|>x->[x...]|>x->x[~isnan(x)]|>x->(@assert(length(x)==1);x[1])
    if outall
        pval12=callr("anova",md2,md1,test=test)|>x->relem(x,length(x))|>x->[x...]|>x->x[~isnan(x)]|>x->(@assert(length(x)==1);x[1])
        return pval1,pval2,pval12
    else
        return min(pval1,pval2)
    end
end
export glm12test

#}}

#{{ nthmax nthmin nthmaxh nthminh
#[[ nthmin nthmax nthmaxh nthminh ]]
# N-th-max|min = nthmax|nthmin(N, Vector)
# N-th-max|min = nthmax|nthmin(N, Matrix, dim)
# Vector = nthmaxh|nthminh(N, Matrix) # equal to vec(nthmax|nthmin(N,M,2))
#Calculate the Nth minimum or maximum values.
#When n<=10, running time is linearly increased with N.
#When n>10, function will automatically decide whether using sort.
#See also:
#Xiong Jieyi, March 7, 2015>Jun 6, 2016 > 5 Jan 2024

export nthmin, nthmax, nthminh, nthmaxh
function nthmin(n::Integer,X::AbstractVector{T}) where T<:Real
    if n<=1
        minimum(X)
    elseif n<=10 && n<length(X)*log(length(X))
        Y=collect(X)
        for _=2:n
            deleteat!(Y, argmin(Y))
        end
        minimum(Y)
    else
        sort(X)[n]
    end::T
end
#<v0.6# function nthmin{T}(n::Integer,X::AbstractMatrix{T},dim::Int)
function nthmin(n::Integer,X::AbstractMatrix{T},dim::Int) where {T}
    if dim==1
        T[nthmin(n,X[:,i]) for i=1:size(X,2)]'
    else
        @assert(dim==2,"dim can only be 1 or 2.")
        hcat(T[nthmin(n,vec(X[i,:])) for i=1:size(X,1)])
    end
end
function nthmax(n::Integer,X::AbstractVector{T}) where T<:Real
    if n<=1
        maximum(X)
    elseif n<=10 && n<length(X)*log(length(X))
        Y=collect(X)
        for _=2:n
            deleteat!(Y,argmax(Y))
        end
        maximum(Y)
    else
        sort(X, rev=true)[n]
    end::T
end
#<v0.6# function nthmax{T}(n::Integer,X::AbstractMatrix{T},dim::Int)
function nthmax(n::Integer,X::AbstractMatrix{T},dim::Int) where {T}
    if dim==1
        T[nthmax(n,X[:,i]) for i=1:size(X,2)]'
    else
        @assert(dim==2,"dim can only be 1 or 2.")
        hcat(T[nthmax(n,vec(X[i,:])) for i=1:size(X,1)])
    end
end
nthminh(n::Integer,X::AbstractMatrix)=vec(nthmin(n,X,2))
nthmaxh(n::Integer,X::AbstractMatrix)=vec(nthmax(n,X,2))
#}}

#{{ slidwin
#[[ slidwin ]]
# (fun_out, out_x) = slidwin(fun, Y1, Y2, ...; x::Group=1:N, winlen= ~N/10(oddnum), x_sorted=false)
#Slid-window calculation and output fun(Ys_in_window...) and the mid-x of window. For input Ys and x in row number N, the output will have N-winlen+1 rows. Function output will be pile up by mapr.
#If x is sorted, set x_sorted=true to save time. Note that output will be wrong without any warning if x is not in sort.
#See also: mapr
#Xiong Jieyi, Feb 29, 2016 > 26 Oct 2021

export slidwin
function slidwin(fun, Y::Group...; x::Group=1:rownum(Y), winlen::Integer=floor(Int, length(x)/10)|>x->mod(x,2)==1 ? x : x+1, x_sorted::Bool=x==1:length(x))
    @assert(mod(winlen,2)==1 && winlen>0, "winlen must be a positive odd integer.")
    @assert(iseqrownum(Y..., x), "The 2-last inputs (and x if given) should have the identical row number.")
    if !x_sorted
        x,t=sortr(x)
        Y=getrow(Y,t)
    end
    wing=int((winlen-1)/2)
    mapr(wing+1:length(x)-wing) do i
        (fun(getrow(Y, i-wing:i+wing)...), x[i])
    end
end
#}}

#{{ sumh,medianh,meanh,varh,stdh,maximumh,minimumh,anyh,allh,counth
#[[ sumh medianh meanh varh stdh maximumh maxh minimumh minh anyh allh counth ]]
# sumh(X::AbstractMatrix)=vec(sum(X,2)) #The same for other functions.
# maxh and minh is equal to maximumh and minimumh.
#Just a syntax suger.
#see also: truesbut, falsesbut
#Xiong Jieyi, Jun 24, 2015 >Jan 19, 2016

export sumh,medianh,meanh,varh,stdh,maximumh,maxh,minimumh,minh,anyh,allh,counth
sumh(X::AbstractMatrix)=vec(sum(X, dims=2))
medianh(X::AbstractMatrix)=vec(median(X, dims=2))
meanh(X::AbstractMatrix)=vec(mean(X, dims=2))
varh(X::AbstractMatrix)=vec(var(X, dims=2))
stdh(X::AbstractMatrix)=vec(std(X, dims=2))
maximumh(X::AbstractMatrix)=vec(maximum(X, dims=2))
maxh(X::AbstractMatrix)=vec(maximum(X, dims=2))
minimumh(X::AbstractMatrix)=vec(minimum(X, dims=2))
minh(X::AbstractMatrix)=vec(minimum(X, dims=2))
anyh(X::AbstractMatrix)=vec(any(X, dims=2))
allh(X::AbstractMatrix)=vec(all(X, dims=2))
counth(X::AbstractMatrix)=vec(count(X, dims=2))
#}}

#{{ mx2DataFrame

#[[ mx2DataFrame ]]
# DF =  mx2DataFrame(Matrix; col=vector, row=vector)
#Create a DataFrame from matrix. This funcion is useful to output data tidily on screen.
#See also: pandasDF, tb2DataFrame
#Xiong Jieyi, Feb 13, 2016

export mx2DataFrame
function mx2DataFrame(M::Matrix; col::Vector=map(f"C$1",1:size(M,2)), row::Vector=[])
    df=importpkg(:DataFrames)
    C=splitdimv(M,1)
    if !isempty(row)
        C=Vector[row, C...]
        col=[""; col]
    end
    df.DataFrame(;map((x,y)->Pair(symbol(x),y),col,C)...)
end

#}}

#{{ quantilenorm
#[[ quantilenorm ]]
# M=quantilenorm(Matrix)
# V1,V2,...=quantilenorm(V1,V2,...)
#Quantile normalization for each column of matrix or each inputed vector.
#See also:
#Xiong Jieyi, Aug 25, 2016

export quantilenorm
#<v0.6# function quantilenorm{T<:Real}(M::AbstractMatrix{T})
function quantilenorm(M::AbstractMatrix{T}) where {T<:Real}
    if size(M,2)==1
        return M
    end
    ln=size(M,1)
    IM=hcat(map(sortperm,eachcol(M))...)
    SM=hcat(map((x,l)->x[l],eachcol(M),eachcol(IM))...)
    Sme=rowfun(mean,SM,as_vec=true)
    reshape(Sme[vec(IM)],size(IM))
end
quantilenorm(Vs::AbstractVector...)=splitdim(quantilenorm(hcat(Vs...)),1)
#}}

#{{ princomp
#[[ princomp ]]
#(score_matrix, contribution_precentiages) = princomp(M; method=:pca|:pca_R|:tSNE_py|:mds, tsne_init_pc=0 or 50, mds_n_comp=2, tsne_kw=...)
#Calcuate PCA score matrix and contribution precentages using MultivariateStats package (method=:pca) or R"prcomp" (method=:pca_R).
#This function also support tSNE (method=:tSNE_py) and MDS (method=:mds). In these modes, the 2nd output is always nothing.
#For tSNE,
#    tsne_init_pc assign the initial PC number. By default, tsne_init_pc=50 if size(M, 2)>100, or no pre-PCA otherwise.
#    tsne_kw assigns the arguments of sklearn.manifold.TSNE. See https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
#For MDS, M is a distance matrix (use distvec2mx() to transfer distance vector first). mds_n_comp assign the output column number (default: 2).
#See also: pcaplot, distvec2mx
#Xiong Jieyi, Nov 21, 2016 > 10 Oct 2019

export princomp
function princomp(M::AbstractMatrix{TT}; method::Symbol=:pca, tsne_init_pc::Int=size(M, 2)>100 ? 50 : 0, mds_n_comp::Int=2, tsne_kw=NamedTuple()) where {TT<:Real}
    if method==:pca
        @pkgfun(MultivariateStats, fit=>__mvs_fit, transform=>__mvs_transform, principalvars=>__mvs_principalvars, principalratio=>__mvs_principalratio)
        pca_obj=__mvs_fit(importpkg(:MultivariateStats).PCA, M')
        pcaM=trsp(__mvs_transform(pca_obj,M'))
        pc_contri=__mvs_principalvars(pca_obj)|>x->__mvs_principalratio(pca_obj)*x./sum(x)
    elseif method==:pca_R
        @pkgfun(myR,callr,relemj,relem,r2j)
        pca_Robj=callr(:prcomp,M,center=true,scale=false)
        pcaM=relemj(pca_Robj,"x")
        t=relem(callr(:summary,pca_Robj),"importance")
        rn=r2j(callr(:row!names,t), keepvec=true)
        pc_contri=vec(r2j(t, keepvec=true)[findonly(rn.=="Proportion of Variance"),:])
    elseif method==:tSNE_py
        #According to http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#examples-using-sklearn-manifold-tsne
        if tsne_init_pc>0
            M=princomp(M, method=:pca)[1]
            if size(M, 2)>tsne_init_pc
                M=M[:, 1:tsne_init_pc]
            end
        end
        pcaM = @il importpy("sklearn.manifold").TSNE(; tsne_kw...).fit_transform(M)
        pc_contri=nothing
    elseif method==:mds
        @pkgfun(MultivariateStats, classical_mds)
        pcaM=trsp(classical_mds(M, mds_n_comp))
        pc_contri=nothing
    else
        error("Unsupported method $method.")
    end
    (pcaM, pc_contri)
end
#}}

#{{ fdrcutoff
#[[ fdrcutoff ]]
# pval_cutoff = fdrcutoff(S_test, S_H0; fdr=0.05, max_fdr=min(3*fdr, 1), nowarn=false)
# Calculate the cutoff of S, based on test sets and H0-control sets. count( S_H0 .< cutoff )/count( S_test .< cutoff ) <= given_FDR.
# Function will search from the smallest S until the current FDR reached the max_fdr. If chosen a cutoff is unfeasible, function return zero or the maximum inputed value instead, along with a warning. If all inputed S are below the FDR, function will output the maximum number with a warning. The warning can be repressed by set nowarn=true.
#See also: glm12test
#Xiong Jieyi, Jun 27, 2017 > 21 Oct 2020

export fdrcutoff
function fdrcutoff(P1::Vector{T1}, P0::Vector{T2}; fdr::AbstractFloat=0.05, max_fdr::AbstractFloat=min(3*fdr,1.0), nowarn::Bool=false) where {T1<:Real, T2<:Real}
    (Ps, si)=sortr([P1;P0])
    Ls=[trues(length(P1));falses(length(P0))][si]
    N1=0
    N0=0
    Pcutoff=zero(Ps[1])
    flag=true
    for i=1:length(Ps)
        if Ls[i]
            N1+=1
        else
            N0+=1
        end
        if i<length(Ps) && Ps[i+1]==Ps[i]
            continue
        end
        cP=N0/N1
        if cP<=fdr
            flag=false
            Pcutoff=Ps[i]
        elseif N1>=100 && cP>max_fdr
            break
        end
    end
    if !nowarn
        if flag
            @warn("No any cutoff could satisfy the FDR. fdrcutoff output zero instead.")
        elseif Pcutoff==Ps[end]
            @warn("Any cutoff could satisfy the FDR. fdrcutoff output the maximum inputed value instead.")
        end
    end
    Pcutoff
end
#}}

#{{ zscore
#[[ zscore ]]
# score_mx = zscore( Matrix; dims=0, baseon=X )
# Normalize data to mean=0, std=1. If baseon is assigned, the mean and std will calculated on the assigned data.
# See also:
#Xiong Jieyi, 27 Oct 2017 > 25 Oct 2018

export zscore
function zscore(X::Union{AbstractArray, Base.SkipMissing}, D::Integer=0; baseon::Union{AbstractArray, Base.SkipMissing}=X, dims::Integer=D)
    if dims==0
        (X.-mean(baseon))./std(baseon)
    else
        (X.-mean(baseon, dims=dims))./std(baseon, dims=dims)
    end
end
#}}

#{{ modebyr
#[[ modebyr ]]
# mode_element = modebyr(Group; big=false, sorted=false, sortcheck=true)
# Find the mode, i.e. the most occurred element (by row) in a group-data. If the input is sorted, set sorted=true to speed up.
# If more than one category have the same maximum counts, the one with the smallest group ID (big=false), or the one with the biggest group ID (big=true), will be returned.
# See also: freq, idennum, uninum
# Xiong Jieyi, 5 Nov 2019

export modebyr
function modebyr(X::Group; sorted::Bool=false, sortcheck::Bool=true, big::Bool=false)
    idx=if sorted
        if sortcheck && !issortedr(X)
            error("Input is not sorted.")
        end
        1:rownum(X)
    else
        sortri(X)
    end
    P=1
    maxct=0
    maxI=1
    for i=2:length(idx)
        if !isequalr2(X, idx[i], X, idx[i-1])
            cct=i-P
            P=i
            if big ? (maxct<=cct) : (maxct<cct)
                maxct=cct
                maxI=i-1
            end
        end
    end
    cct=length(idx)-P+1
    if big ? (maxct<=cct) : (maxct<cct)
        maxI=length(idx)
    end
    getrow(X, idx[maxI])
end
#}}

#{{ iselequal
#[[ iselequal ]]
# T|F = iselequal([function, ]Iternation)
#Calculate if all the elements, or all the outputs of fun(element), are equal (isequal() return true). Empty or one-element iternation is always return true.
#See also:
#Xiong Jieyi, 9 Jul 2018

export iselequal
if VERSION<v"0.7.0"
    function iselequal(X)
        st=start(X)
        done(X, st) && return true
        x1, st=next(X, st)
        while !done(X, st)
            cx, st=next(X, st)
            isequal(cx, x1) || return false
        end
        return true
    end
else
    function iselequal(X)
        out=iterate(X)
        if out==nothing
            return true
        else
            x1, st=out
        end
        while true
            out=iterate(X, st)
            out==nothing && return true
            cx, st=out
            isequal(cx, x1) || return false
        end
    end
end
iselequal(fun::Function, X)=iselequal(fun(cx) for cx in X)
#}}

#{{ valid
#[[ valid ]]
# Y1, Y2, ... = valid(X1, X2, ..., Xn; delif=Number:nan|inf or OtherType:missing|nothing|empty, follower|onlyon=index(1-n)_vec|index_set|TF_vec(length=n) )
# Remove the invalid values (missing, NaN and Inf by default, or any others where delif return true) in float vectors. Inputs should be vectors in the same lengths. The element "invalid" in any of input vectors will be removed. The outputs are always in the same lengths.
# The #i input where i in `follower` will not be checked by `delif`, but only be filtered at last.
# onlyon=... is the reverse of follower=.... The two keyword arguments cannot be assigned simultaneously.
# It is useful such as cor(valid(X, Y)...). Avoid to use it in the preformance-critical codes.
#See also:
#Xiong Jieyi, 17 Aug 2018 > 30 Jan 2019>21 Feb 2019>30 Jan 2023

export valid
function valid(X1::AbstractVector, Xs::AbstractVector...; delif::Function=x->isa(x, Number) ? (isnan(x) || isinf(x)) : (ismissing(x) || isnothing(x) || isempty(x)), follower=Set(), onlyon=nothing)
    if !isnothing(onlyon)
        @assert(isempty(follower), "follower=... and onlyon=... cannot be assigned simultaneously.")
        follower=if isa(onlyon, AbstractVector{Bool})
            Set(findall(.!onlyon))
        else
            setdiff!(Set(1:(length(Xs)+1)), onlyon)
        end
    elseif isa(follower, AbstractVector{Bool})
        follower=Set(findall(follower))
    elseif !isa(follower, AbstractSet)
        follower=Set(follower)
    end
    l=in(1, follower) ? falses(length(X1)) : delif.(X1)
    if isempty(Xs)
        return X1[.!l]
    else
        for i in 1:length(Xs)
            if !in(i+1, follower)
                l=l .| delif.(Xs[i])
            end
        end
        return X1[.!l], map(x->x[.!l], Xs)...
    end
end
#}}

#{{ grp2mx
#[[ grp2mx ]]
# * grp2mx() will be disused in the future. Try coo2mx() instead. 8 Oct 2020 *
# RowId, ColId, Mx1, Mx2, ... = grp2mx(RowGrp, ColGrp, V1[=>default], V2[=>default], ...;
#                               (option I)    stable|row_stable|col_stable=false,
#                               (option II)   row_order=..., col_order=...) #Should be row-unique.
# Convert a table like vector to a matrix.
# If the '=>default' is not given, the default value will be the cooresponding emptyval. Note that input as fastgrp object is not supported. default value should be in the same type as V, or be missing.
# Note that this function will not check 1) whether (RowGrp, ColGrp) are row-unique; 2) whether row_order/col_order are row-unique. For 2, the value will NOT be duplicated but just assigned with default values in all-but-one row/columns. Utilize this feature is not recommended. These two problems have been fixed in coo2mx().
# See also: grpvec, fullfactors, coo2mx
# Xiong Jieyi, 3 Jun 2019

export grp2mx
function grp2mx(R::Group, C::Group, Vs...;
                stable::Bool=false, row_stable::Bool=stable, col_stable::Bool=stable,
                row_order::Union{Group, Nothing}=nothing,
                col_order::Union{Group, Nothing}=nothing)
    if row_order==nothing
        ri, rl=uniqr(R; stable=row_stable)
        uR=getrow(R, ri)
        rlen=length(ri)
    else
        uR=row_order
        rl=memberr(R, uR)
        rlen=rownum(uR)
    end
    if col_order==nothing
        ci, cl=uniqr(C; stable=col_stable)
        uC=getrow(C, ci)
        clen=length(ci)
    else
        uC=col_order
        cl=memberr(C, uC)
        clen=rownum(uC)
    end
    Ms=map(Vs) do V
        if isa(V, Pair)
            df=V[2]
            V=V[1]
        else
            df=emptyval(eltype(V))
        end
        M=fill(df, rlen, clen)
        if ismissing(df)
            M=convert(Matrix{Union{eltype(V), Missing}}, M)
        end
        foreach(rl, cl, V) do r, c, v
            if r>0 && c>0
                M[r, c]=v
            end
        end
        M
    end
    (uR, uC, Ms...)
end
#}}

#{{ fullfactors
#[[ fullfactors ]]
# all_factor_comb = fullfactors((factor1, ...)[ =>(full_factor1|:, ...) ])
# all_factor_comb, filled_values = fullfactors(factors[=>full_factors], values, default|nothing; safe=true)
# Table = fullfactors(Tabel, c"fields"[=>all_factors], default; safe=true)
# Full all the combinations of given factors, fill other data-in-row using default value (or empty value when the 3nd input is 'nothing').
# Inputed factors and full_factors MUST be Tuples, even only one factor (in this case, using dtshift() is more recommended). The 'values', 'default' and 'safe=...' have the same meaning as dtshift() #3 and #4 inputs. default could be 'nothing' for emptyval. Note that safe=true is the default, which is different from dtshift(), and the 3rd input is always required unless no any values need to be filled.
# all_factors can have colon (:) element, meaning the unival() of the cooresponding factor.
# For the output, the loop order is the same as eachcombr().
# e.g. (V1, V2) = fullfactors((iV1, iV2))
#      (V1, V2), (V3, V4, V5) = fullfactors((iV1, iV2), (iV3, iV4, iV5), (0, "n.a.", 0.0))
#      Table = fullfactors(T, c"key1, key2"=>(1:4, :), tb(key3=0, key4="NA"))
#
# In below examples, key1 and key2 are considered as a join-group, and its pairs will not be recombind.
#      Table = fullfactors(T, (("key1", "key2"), "key3"), tb(key4=0))
#      Table = fullfactors(T, (c"key1, key2", "key3"), tb(key4=0))
# Also the same as:
#     ((V1, V2), V3), ... = fullfactors(((V1, V2), V3)=>..., ...)
#See also: dtshift, coo2mx, eachcombr, mxexpand
#Xiong Jieyi, 27 Jan 2020

export fullfactors
function fullfactors(iFs::Pair{T1, T2}) where {T1<:Tuple, T2<:Tuple}
    Fs=iFs.first
    Ds=map(iFs.first, iFs.second) do cF, cD
        if cD==(:)
            unival(cF)
        else
            cD
        end
    end
    vcatr(eachcombr(Ds...)...)
end
fullfactors(iFs::Tuple, arg...; kw...)=fullfactors(iFs=>tuple(fill(:, length(iFs))...), arg...; kw...)
function fullfactors(iFs::Pair, V::Group, default; safe::Bool=true)
    function __flatten_2layers(T1, T2)
        Ys=(Any[], Any[])
        for X in zip(T1, T2)
            if isa(X[1], Tuple)
                for xx in zip(X...)
                    for (y, x) in zip(Ys, xx)
                        push!(y, x)
                    end
                end
            else
                for (y, x) in zip(Ys, X)
                    push!(y, x)
                end
            end
        end
        Tuple(Ys[1]), Tuple(Ys[2])
    end
    
    vF=fullfactors(iFs)
    fF, fI=__flatten_2layers(vF, iFs.first)
    (vF, dtshift(fF, fI, V, default; safe=safe))
end
function fullfactors(T::Dict{Tk,}, ky::Pair, arg...; kw...) where {Tk<:AbstractString}
    function __flatten_2layers(T1, T2)
        Ys=(AbstractString[], Any[])
        for X in zip(T1, T2)
            if isa(X[1], Tuple) || isa(X[1], AbstractVector)
                for xx in zip(X...)
                    for (y, x) in zip(Ys, xx)
                        push!(y, x)
                    end
                end
            else
                for (y, x) in zip(Ys, X)
                    push!(y, x)
                end
            end
        end
        Ys[1], Tuple(Ys[2])
    end
    
    vF=fullfactors(dict2tuple(T, ky.first)=>ky.second) #vF could be nested tuple.
    (fK, fF)=__flatten_2layers(ky.first, vF) #fF and fK is not nested.
    dT=fd_i(T, fK)
    oT=isempty(dT) ? tb() : dtshift(fF, dict2tuple(T, fK), dT, arg...; kw...)
    tb(oT, fK=>fF)
end
fullfactors(T::Dict{Tk,}, ky::Union{AbstractVector{Tv}, Tuple}, arg...; kw...) where {Tk<:AbstractString, Tv<:AbstractString} = fullfactors(T, ky=>tuple(fill(:, length(ky))...), arg...; kw...)

#}}

#{{ pysparse2coo pysparse2jl
#[[ pysparse2coo pysparse2jl ]]
# rowIdx, colIdx, value = pysparse2coo( python_sparseMatrix_obj ) # indice are 1-based.
# julia_sparse_mx = pysparse2jl( python_sparseMatrix_obj ) # `using SparseArrays` required.
# it_self = pysparse2jl( julia_mx ) # Just for convenience.
# Convert python sparse matrix to three vectors. The output indice are 1 based. PyCall should be preloaded.
# See also: coo2mx
# Jieyi Xiong 17 Jun 2019

export pysparse2coo, pysparse2jl
function pysparse2coo(X)
    # ri, ci, v=importpy("scipy.sparse").find(X::Main.PyCall.PyObject)
    # (ri.+1, ci.+1, v)
    X=(X::Main.PyCall.PyObject).tocoo()
    (X.row.+1, X.col.+1, X.data)
end
pysparse2jl(X)=invokelatest(importpkg(:SparseArrays).sparse, pysparse2coo(X)..., X.shape...)
pysparse2jl(X::AbstractArray)=X
#}}

#{{ linkgrp
#[[ linkgrp ]]
# grp_no = linkgrp( int_matrix )
# Link pairs to groups. The input is a Int64 matrix, each row is a paired link (it actually supported more than two columns). The function output a Int64 vector with the same length as the row number of inputted matrix. Two rows will has the same output numbers if, but not only if, any of their elements are identical. For two different output groups, all their elements in the input matrix are different.
# The zeors in inputted matrix will be ignored and not be regarded as a pair information.
# See also: rankgrp
# Xiong Jieyi, 27 Mar 2020

export linkgrp
function linkgrp(X::AbstractMatrix{Int})
    L=X.>0
    M=zeros(Int, size(X))
    M[L].=uniqr(X[L])[2]
    maxV=maximum(M[L])
    V=collect(1:maxV)
    for m in Base.eachrow(M)
        m=m[m.>0]
        mm=minimum(m)
        xv=copy(m)
        for x in m
            while V[x]<x
                x=V[x]
                push!(xv, x)
            end
            mm=min(mm, x)
        end
        V[xv].=mm
    end
    B=falses(length(V))
    for x=length(V):-1:1
        B[x] && continue
        xv=[x]
        while V[x]<x
            push!(xv, x)
            x=V[x]
        end
        V[xv].=x
        B[xv].=true
    end
    t=rowfun(M, as_vec=true) do x
        V[x[findfirst(x->x>0, x)]]
    end
    uniqr(t)[2]
end

#Below function should have exact the same beheave as linkgrp(). I kept it here for validation. 27 Mar 2020
function linkgrp_bypkg(X::AbstractMatrix{Int})
    @pkgfun(LightGraphs, SimpleGraph, add_edge!, connected_components)
    L=X.>0
    M=zeros(Int, size(X))
    M[L].=uniqr(X[L])[2]
    maxV=maximum(M[L])
    G=SimpleGraph(maxV)
    for m in Base.eachrow(M)
        m=m[m.>0]
        for mm in m[2:end]
            add_edge!(G, m[1], mm)
        end
    end
    C=connected_components(G)
    Gid, elem=vcatr_with(1:length(C), C...)
    O=zeros(Int, size(M))
    O[L]=dtshift(M[L], elem, Gid)
    rowfun(O, as_vec=true) do x
        for cx in x
            if cx>0
                return cx
            end
        end
        return 0
    end
end

#}}

#{{ setcmpr
#[[ setcmpr ]]
# setcmpr(A::Group, B::Group[, C::Group]; venn=false, label=c"A, B[, C]")
#Print the shared relationship of Row-sets A and B (or among there groups). It has no output value.
#If venn is true, function will also draw a venn figure using matplotlib_venn python package.
#See also: ls, quantilels, eachr, intersectri, ismbr, venn
#Xiong Jieyi, May 31, 2015 > Jan 14, 2016 > 3 Sep 2020 > 7 Oct 2020

export setcmpr
function setcmpr(A::Group, B::Group; venn::Bool=false, label=c"A, B")

    #Added in 22 Sep 2018, only for user error prevention.
    if xor(typeof(A)==Any, typeof(B)==Any)
        @warn("Only one of the inputs is Any[].")
    end

    A=Set(isa(A, AbstractVector) ? A : eachr(A))
    Al=length(A)
    @printf("A uni-num: %d\n", Al)
    B=Set(isa(B, AbstractVector) ? B : eachr(B))
    Ao=length(setdiff(A,B))
    @printf("A only: %d (%f%%)\n", Ao, 100*Ao/Al)
    AB=length(intersect(A,B))
    Bl=length(B)
    @printf("AB shared: %d (%f%% in A, %f%% in B)\n", AB, 100*AB/Al, 100*AB/Bl)
    Bo=length(setdiff(B,A))
    @printf("B only: %d (%f%%)\n", Bo, 100*Bo/Bl)
    @printf("B uni-num: %d\n", length(B))

    if venn
        ob=importpy(:matplotlib_venn)
        ob.venn2(subsets=ds("10"=>Ao, "01"=>Bo, "11"=>AB), set_labels=label)
    else
        nothing
    end
end
function setcmpr(A::Group, B::Group, C::Group; venn::Bool=false, label=c"A, B, C")

    if 1 .<= count([typeof(A)==Any, typeof(B)==Any, typeof(C)==Any]) .<=3
        @warn("Only one or two inputs are Any[].")
    end

    A=Set(isa(A, AbstractVector) ? A : eachr(A))
    Al=length(A)
    @printf(" A uni-num: %d\n", Al)

    B=Set(isa(B, AbstractVector) ? B : eachr(B))
    Bl=length(B)
    @printf(" B uni-num: %d\n", Bl)
    
    C=Set(isa(C, AbstractVector) ? C : eachr(C))
    Cl=length(C)
    @printf(" C uni-num: %d\n", Cl)

    AB=intersect(A, B)
    ABl=length(AB)
    @printf("AB  shared: %d (%f%% in A, %f%% in B)\n", ABl, 100*ABl/Al, 100*ABl/Bl)

    AC=intersect(A, C)
    ACl=length(AC)
    @printf("AC  shared: %d (%f%% in A, %f%% in C)\n", ACl, 100*ACl/Al, 100*ACl/Cl)

    BC=intersect(B, C)
    BCl=length(BC)
    @printf("BC  shared: %d (%f%% in B, %f%% in C)\n", BCl, 100*BCl/Bl, 100*BCl/Cl)

    ABC=intersect(AB, BC)
    ABCl=length(ABC)
    @printf("ABC shared: %d (%f%% in A, %f%% in B, %f%% in C)\n", ABCl, 100*ABCl/Al, 100*ABCl/Bl, 100*ABCl/Cl)

    uA=setdiff(A, union(B, C))
    uAl=length(uA)
    @printf(" A  only  : %d (%f%%)\n", uAl, 100*uAl/Al)

    uB=setdiff(B, union(A, C))
    uBl=length(uB)
    @printf(" B  only  : %d (%f%%)\n", uBl, 100*uBl/Bl)

    uC=setdiff(C, union(A, B))
    uCl=length(uC)
    @printf(" C  only  : %d (%f%%)\n", uCl, 100*uCl/Cl)
    
    if venn
        importpkg(:PyPlot, preloaded=true)
        ob=importpy(:matplotlib_venn)
        ob.venn3(subsets=ds("100"=>uAl, "010"=>uBl, "110"=>ABl-ABCl,
                            "001"=>uCl, "101"=>ACl-ABCl, "011"=>BCl-ABCl,
                            "111"=>ABCl), set_labels=label)
    else
        nothing
    end
end

#}}

#{{ mxasvecdo
#[[ mxasvecdo ]]
# Matrix = mxasvecdo(fun, Matrix)
# Matrix = mxasvecdo(fun, Matrix1, Matrix2, ...) #fun has multiple vector inputs.
# Matrix1, Matrix2, ... = asvecdo(...) #fun has multiple vector outputs.
# The inputed matrix will be convert as vector then as arguments of function. The result will convert back to Matrix. If the function output is a Tuple, every vector in the tuple will be converted.
# Note: in multiple inputs mode, the shape of output matrix(es) is only based on the first inputed matrix. The shapes of other inputs will NOT be checked.
#See also: vcatdo, mxexpand, coo2mx
#Xiong Jieyi, Sep 18, 2015 > 7 Oct 2020

export mxasvecdo
function mxasvecdo(fun, X1::AbstractMatrix, Xs::AbstractMatrix...)
    O=if isempty(Xs)
        fun(vec(X1))
    else
        Xv=map(vec, tuple(X1, Xs...))
        fun(Xv...)
    end
    sz=size(X1)
    if isa(O, Tuple)
        map(O) do x
            @assert(isa(x, AbstractVector), "Function output is not a vector.")
            reshape(x, sz)
        end
    else
        @assert(isa(O, AbstractVector), "Function output is not a vector.")
        reshape(O, sz)
    end
end
#}}

#{{ vcatdo
#[[ vcatdo ]]
# (O1, O2, ...) = vcatdo(function, I1, I2, ...)
# Do as vcatr(O1,O2,...)=function( vcatr(I1,I2,...)), i.e., merge all in-s by row, apply the function, then split the output by row.
#See also: mxasvecdo
#Xiong Jieyi, Sep 20, 2015 > 7 Oct 2020

export vcatdo
function vcatdo(fun, X1::Group, Xs::Group...)
    L=1:length(Xs)+1
    idx, V=vcatr_with(L, X1, Xs...)
    W=fun(V)
    Tuple(grpvec(fastgrp(idx, L), W))
end
#}}

#{{ TaroneZtest
#[[ TaroneZtest ]]
# Z, Pvalue, ratio = TaroneZtest(A, B)
# Perform Tarone's Z test, code are adopted from https://stats.stackexchange.com/a/410376/6378
# Input A and B are vectors with the same lengths. Each of elements are the numbers of success and failed trials in total of A.+B trials. The outputted ratio is sum(A)/sum(A.+B).
# A large Z or a small Pvalue means beta-binomial model is more fitted than binomial model for this data.
# Xiong Jieyi, 10 Jan 2023

export TaroneZtest
function TaroneZtest(M::AbstractVector{<:Integer}, P::AbstractVector{<:Integer})
    (any(P.<0) || any(M.<0)) && error("Inputs should be positive.")
    length(P)==length(M) || error("Inputs have different lengths.")
    N=P.+M
    estimate = sum(M)/sum(N)
    S = (estimate == 1) ? sum(N) : sum((M - N.*estimate).^2 ./ (estimate*(1 - estimate)))
    statistic = (S .- sum(N))./sqrt(2*sum(N.*(N.-1)))
    @pkgfun(myR, callr1)
    pvalue = 2*callr1("pnorm", -abs(statistic), 0, 1)
    (statistic, pvalue, estimate)
end
#}}

#{{ Some common R functions
export wilcoxtest, studenttest, padjust, cortest, fishertest, fishertest_l, binomtest, binomtest_l
#[[ wilcoxtest ]]
# p.value = wilcoxtest(x, y = NULL,
#     alternative = c("two.sided", "less", "greater"),
#     mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
#     conf.int = FALSE, conf.level = 0.95, ...;
#     vb=0, keywordcheck=true)
# Performs one- and two-sample Wilcoxon tests on vectors of data; the latter is also known as 'Mann-Whitney' test.
# x: numeric vector of data values.  Non-finite (e.g., infinite or missing) values will be omitted.
# y: an optional numeric vector of data values: as with 'x' non-finite values will be omitted.
# alternative: a character string specifying the alternative hypothesis, must be one of '"two.sided"' (default), '"greater"' or '"less"'.  You can specify just the initial letter.
# paired: a logical indicating whether you want a paired test.
#For detail, see rhelp("wilcox.test")
#See also: studenttest, padjust, cortest, fishertest, binomtest
#Xiong Jieyi, 27 Oct 2017

function _checkkeyword(arg, kwset)
    for (k, _) in arg
        if !in(k, kwset)
            error("Invalid keyword argument: $k")
        end
    end
end

function wilcoxtest(arg...; vb::Int=0, keywordcheck::Bool=true, warg...)
    keywordcheck && _checkkeyword(warg, Set([:x, :y, :mu, :paired, :exact, :correct, :conf!int, :conf!level, :formula, :data, :subset, :na!action]))
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("wilcox.test", arg...; warg...)
    rshow(o; vb=vb)
    r2j(relem(o, "p.value"), NA=NaN, keepvec=false)
end

#[[ studenttest ]]
# p.value, meanX - meanY, [CI_low, CI_high] = studenttest(x, y = NULL,
#        alternative = c("two.sided", "less", "greater"),
#        mu = 0, paired = FALSE, var.equal = FALSE,
#        conf.level = 0.95, ...; vb=0, keywordcheck=true)
# Performs one and two sample t-tests on vectors of data.
# x: a (non-empty) numeric vector of data values.
# y: an optional (non-empty) numeric vector of data values.
# paired: a logical indicating whether you want a paired t-test.
# For detail, see rhelp("t.test")
# See also: studenttest, padjust, cortest, fishertest, binomtest
# Xiong Jieyi, 27 Oct 2017

function studenttest(arg...; vb::Int=0, keywordcheck::Bool=true, warg...)
    keywordcheck && _checkkeyword(warg, Set([:x, :y, :alternative, :mu, :paired, :var!equal, :conf!level, :formula, :data, :subset, :na!action]))
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("t.test", arg...; warg...)
    rshow(o; vb=vb)
    pval=r2j(relem(o, "p.value"), NA=NaN, keepvec=false)
    mn=r2j(relem(o, "estimate"), NA=NaN, keepvec=true)
    dmn=if length(mn)==1 # if paired=true
        mn[1]
    else
        @assert length(mn)==2
        mn[1]-mn[2]
    end
    cfi=r2j(relem(o, "conf.int"), NA=NaN, keepvec=true)
    (pval, dmn, cfi)
end
#[[ padjust ]]
#  padj = padjust(p, method = p.adjust.methods)
#  Given a set of p-values, returns p-values adjusted using one of several methods.
#  p: numeric vector of p-values (possibly with 'NA's).  Any other R is coerced by 'as.numeric'.
#  p.adjust.methods: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
# For detail, see rhelp("p.adjust")
# See also: wilcoxtest, studenttest, cortest, fishertest, binomtest
# Xiong Jieyi, 27 Oct 2017

function padjust(arg...; warg...)
    #p.adjust does not need to check keyword.
    @pkgfun(myR, callr, r2j)
    r2j(callr("p.adjust", arg...; warg...), NA=NaN, keepvec=true)
end

#[[ cortest ]]
# r, p.value = cortest(x, y,
#          alternative = c("two.sided", "less", "greater"),
#          method = c("pearson", "kendall", "spearman"),
#          exact = NULL, conf.level = 0.95, continuity = FALSE, ...;
#          vb=0, keywordcheck=true)
# Test for association between paired samples, using one of Pearson's product moment correlation coefficient, Kendall's tau or Spearman's rho.
# x, y: numeric vectors of data values.  'x' and 'y' must have the same length.
# method: a character string indicating which correlation coefficient is to be used for the test.  One of '"pearson"', '"kendall"', or '"spearman"', can be abbreviated.
# For detail, see rhelp("cor.test")
# See also: wilcoxtest, studenttest, padjust, fishertest, binomtest
# Xiong Jieyi, 27 Oct 2017

function cortest(arg...; vb::Int=0, keywordcheck::Bool=true, warg...)
    keywordcheck && _checkkeyword(warg, Set([:x, :y, :alternative, :method, :exact, :conf!level, :continuity, :formula, :data, :subset, :na!action]))
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("cor.test", arg...; warg...)
    rshow(o; vb=vb)
    r2j.(relem(o, "estimate", "p.value"), NA=NaN, keepvec=false)
end

#[[ fishertest fishertest_l ]]
# p.value = fishertest(x, y = NULL, workspace = 200000, hybrid = FALSE,
#   control = list(), or = 1, alternative = "two.sided(default)|greater|less"",
#   conf.int = TRUE, conf.level = 0.95,
#   simulate.p.value = FALSE, B = 2000; vb=0)
# ... = fishertest_l(A::bool_vec, B::bool_vec; ...) #Note the different function name.
# Performs Fisher's exact test for testing the null of independence of rows and columns in a contingency table with fixed marginals.
# x: either a two-dimensional contingency table in matrix form, or a factor object.
# y: a factor object; ignored if 'x' is a matrix.
# alternative: indicates the alternative hypothesis and must be one of '"two.sided"', '"greater"' or '"less"'.  You can specify just the initial letter.  Only used in the 2 by 2 case.
# For detail, see rhelp("fisher.test")
# See also: wilcoxtest, studenttest, padjust, cortest, fisherexacttest, binomtest
# Xiong Jieyi, 27 Oct 2017 > 12 Apr 2022

function fishertest(arg...; vb::Int=0, warg...)
    #fisher.test does not need to check keyword.
    @pkgfun(myR, callr, relem, r2j, rshow)
    o=callr("fisher.test", arg...; warg...)
    rshow(o; vb=vb)
    r2j(relem(o, "p.value"), NA=NaN, keepvec=false)
end
fishertest_l(A::AbstractVector{Bool}, B::AbstractVector{Bool}; kw...)=fishertest([count(A) count(B); count(.!A) count(.!B)]; kw...)
fishertest(::AbstractVector{Bool}, arg...; kw...)=error("Invalid inputs. Do you mean fishertest_l(bool_vec, bool_vec; ...)?")


#[[ binomtest binomtest_l ]]
# p.value, conf.int = binomtest(x, n, p = 0.5,
#            alternative = c("two.sided", "less", "greater"),
#            conf.level = 0.95; vb=0)
# ... = binomtest_l(X::bool_vec; ...) #equal to: binomtest(count(X), length(X); ...). Note the different function name.
# binomtest(682, 682 + 243, p = 3/4) is equal to binomtest([682, 243], p = 3/4).
# Performs an exact test of a simple null hypothesis about the probability of success in a Bernoulli experiment.
# x: number of successes, or a vector of length 2 giving the numbers of successes and failures, respectively.
# n: number of trials; ignored if x has length 2.
# p: hypothesized probability of success.
# alternative: indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter.
# conf.level: confidence level for the returned confidence interval.
# For detail, see rhelp("binom.test")
# See also: wilcoxtest, studenttest, padjust, cortest, fishertest
# Xiong Jieyi, 27 Oct 2017 > 12 Apr 2022

function binomtest(arg...; vb::Int=0, warg...)
    #binom.test does not need to check keyword.
    @pkgfun(myR, callrw, relem, r2j, rshow)
    ro=callrw("binom.test", arg...; warg...)
    rshow(ro; vb=vb)
    (r2j(relem(ro, "p.value"), NA=NaN, keepvec=false), r2j(relem(ro, "conf.int"), NA=NaN, keepvec=true))
end
binomtest_l(X::AbstractVector{Bool}; kw...)=binomtest(count(X), length(X); kw...)
binomtest(X::AbstractVector{Bool}, arg...; kw...)=error("Invalid inputs. Do you mean binomtest_l(bool_vec; ...)?")
#}}
