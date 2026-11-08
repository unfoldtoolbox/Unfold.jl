
# Term.jl always emits ANSI color codes, even into IOs that don't support color
# (repr, Pluto.jl), where they turn into raw escape garbage. All Term-based show
# methods below are gated on `get(io, :color, false)` via `term_show`. Documenter.jl
# does support colors (it converts ANSI to HTML), the REPL and `:color=>true` IOs
# get the pretty version, everything else clean plain text.
function term_show(f, io::IO)
    if get(io, :color, false)
        f(io)
    else
        old = Term.NOCOLOR[]
        Term.NOCOLOR[] = true
        try
            buf = IOBuffer()
            f(buf)
            print(io, Term.cleantext(String(take!(buf))))
        finally
            Term.NOCOLOR[] = old
        end
    end
end


#----- TimeExpandedTerm
function Base.show(io::IO, p::TimeExpandedTerm)
    term_show(io) do io
        t = times(p.basisfunction)
        tPlot = show_shorten_vector(t)

        tprint(
            io,
            Term.highlight("$(p.basisfunction.name)", :operator) *
            Term.highlight(": timeexpand($(p.term)) for times $(tPlot)"),
        )
    end
end
function Base.show(io::IO, f::FormulaTerm{<:Union{<:InterceptTerm,TimeExpandedTerm}})
    term_show(io) do io
        tprint(io, f.rhs)
    end
end

function Base.show(io::IO, f::Vector{<:FormulaTerm}; eventnames = repeat([""], length(f)))
    term_show(io) do io
        for k = 1:length(f)
            Term.tprintln(
                io,
                "{bold blue_light} $(eventnames[k]){/bold blue_light}",
                "=>",
                f[k],
            )
        end
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
    term_show(io) do io
        print(io, renderable(obj))
    end
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
        :shift_onset => shift_onset(obj),
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
        obj.scale_duration == false ? kernel(obj, 0) : kernel(obj, [0, width(obj) ÷ 2]),
    )

    s = String(take!(io))
    return Panel(str) * s
end


#---- UnfoldModel
function Base.show(io::IO, ::MIME"text/plain", obj::T) where {T<:UnfoldModel}
    term_show(io) do io
        Term.tprintln(io, "Unfold-Type: ::$(typeof(obj)){{$(typeof(obj).parameters[1])}}")

        keys = first.(design(obj))
        forms = formulas(obj)

        empty_uf = T()
        is_not_fit = length(coef(obj)) == 0#modelfit(obj) == modelfit(empty_uf)

        #@debug typeof(formulas(obj))
        Base.show(io, formulas(obj); eventnames = keys)
        println(io)
        Term.tprintln(
            io,
            is_not_fit ? "{red}❌{/red} model not fit. Use fit!(uf,data) to fit it." :
            "{bold green}✔{/bold green} model is fit. {blue_light} size(coefs) $(size(coef(obj))){/blue_light}",
        )
        #println(io, "Design: $(design(obj))")
        println(io)

        tprint(
            io,
            "{gray}Useful functions:{/gray} `design(uf)`, `designmatrix(uf)`, `coef(uf)`, `coeftable(uf)`",
        )
    end
end

```
print design in a beautiful way
```
function Base.show(io::IO, design::Vector{<:Pair{<:Any,<:Tuple}})
    term_show(io) do io
        basisList = []
        for (first, second) in design
            push!(
                basisList,
                Panel(string(first) * "=>" * "($(second[1]),$(second[2]))", fit = false),
            )
        end
        print(io, Term.grid(basisList))
    end
end


function Base.show(io::IO, ::MIME"text/plain", d::Vector{<:AbstractDesignMatrix})
    term_show(io) do io
        Term.tprintln(io, "Vector of length $(length(d)), $(unique(typeof.(d)))")
        Term.tprintln(io, "$(formulas(d))")
        Term.tprintln(
            io,
            "\nuseful functions: `formulas(d)`, `modelmatrix(d)/modelmatrices(d)`, `events(d)`",
        )
    end
end
function Base.show(io::IO, d::AbstractDesignMatrix)
    term_show(io) do io
        Term.tprintln(io, "::$(typeof(d))")
        Term.tprintln(io, "@formula: $(formulas(d))")

        sz_evts = isa(d.events, Vector) ? size.(d.events) : size(d.events)
        sz_modelmatrix =
            (isa(d.modelmatrix, Vector) | isa(d.modelmatrix, Tuple)) ? size.(d.modelmatrix) :
            size(d.modelmatrix)

        Term.tprintln(io, "")
        Term.tprintln(io, "- modelmatrix (time/trials x predictors): $sz_modelmatrix")
        Term.tprintln(io, "- $(sz_evts[1]) events with $(sz_evts[2]) columns")
        Term.tprintln(
            io,
            "\nuseful functions: `formulas(d)`, `modelmatrix(d)/modelmatrices(d)`, `events(d)`",
        )
        Term.tprintln(io, "Fields: `.formula`, `.modelmatrix`, `.events`")
    end
end
