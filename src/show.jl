
#----- TimeExpandedTerm
function Base.show(io::IO, p::TimeExpandedTerm)
    t = times(p.basisfunction)
    tPlot = show_shorten_vector(t)

    tprint(
        io,
        Term.highlight("$(p.basisfunction.name)", :operator) *
        Term.highlight(": timeexpand($(p.term)) for times $(tPlot)"),
    )

end
function Base.show(io::IO, f::FormulaTerm{<:Union{<:InterceptTerm,TimeExpandedTerm}})
    tprint(io, f.rhs)

end

function Base.show(io::IO, f::Vector{<:FormulaTerm}; eventnames = repeat([""], length(f)))
    for k = 1:length(f)
        Term.tprintln(
            io,
            "{bold blue_light} $(eventnames[k]){/bold blue_light}",
            "=>",
            f[k],
        )
    end
end

function show_shorten_vector(t)
    if length(t) > 3
        t = (round.(t[[1, 2, end]]; sigdigits = 4))
        tPlot = "[" * string(t[1]) * ", " * string(t[2]) * " ... " * string(t[3]) * "]"
    else
        tPlot = string(t)
    end
    return tPlot
end

function Base.show(io::IO, obj::BasisFunction)
    print(io, renderable(obj))

end

function renderable(obj::BasisFunction; title = "::BasisFunction")
    d = OrderedDict(
        :name => name(obj),
        :kerneltype => typeof(obj),
        :width => width(obj),
        :height => height(obj),
        :colnames => colnames(obj) |> show_shorten_vector,
        :times => times(obj) |> show_shorten_vector,
        :collabel => collabel(obj),
        :shiftOnset => shiftOnset(obj),
    )# |> x -> Tree(x; title = title)

    str = "{bold blue}::BasisFunction{/bold blue}\n"
    for (key, val) in d
        str =
            str *
            "{bold blue_light}$key: {/bold blue_light}" *
            Term.highlight("$val", :number) *
            "\n"
    end

    io = IOBuffer()
    show(
        IOContext(io, :limit => true, :displaysize => (40, 40)),
        "text/plain",
        kernel(obj, 0.5),
    )

    s = String(take!(io))
    return Panel(str) * s
end


#---- UnfoldModel
function Base.show(io::IO, ::MIME"text/plain", obj::T) where {T<:UnfoldModel}
    Term.tprintln(io, "Unfold-Type: ::$(typeof(obj)){{$(typeof(obj).parameters[1])}}")

    keys = first.(design(obj))
    forms = formulas(obj)


    empty_uf = T()
    is_not_fit = modelfit(obj) == modelfit(empty_uf)



    #@debug typeof(formulas(obj))
    Base.show(io, formulas(obj); eventnames = keys)
    println(io)
    Term.tprintln(
        io,
        is_not_fit ? "{red}❌{/red]} model not fit" :
        "{bold green}✔{/bold green} model is fit. {blue_light} Parameters: $(size(coef(obj))){/blue_light}",
    )
    #println(io, "Design: $(design(obj))")
    println(io)

    tprint(
        "{gray}Useful functions:{/gray} `design(uf)`, `designmatrix(uf)`, `coef(uf)`, `coeftable(uf)`",
    )

end

```
print design in a beautiful way
```
function Base.show(io::IO, design::Vector{<:Pair{<:Any,<:Tuple}})
    basisList = []
    for (first, second) in design
        push!(basisList, Panel(first * "=>" * "($(second[1]),$(second[2]))", fit = false))
    end
    print(io, Term.grid(basisList))
end

# function Base.show(io::IO, obj::UnfoldModel)
#     println(io, "Unfold-Type: $(typeof(obj)) \n")
#     println(io, "formula: $(obj.design)")
#     println(
#         io,
#         "Useful functions:\n 
#     design(uf) \t\t(returns Dict of event => (formula,times/basis))  \n
#     designmatrix(uf) \t(returns DesignMatrix with events) \n
#     modelfit(uf) \t\t(returns modelfit object) \n
#     coeftable(uf) \t\t(returns tidy result dataframe) \n",
#     )
# end
