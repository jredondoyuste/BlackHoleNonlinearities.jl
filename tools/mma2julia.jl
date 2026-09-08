const SECTORS = ["EEE", "EEO", "EOE", "EOO", "OOE", "OOO"]
const MMA_DIR = joinpath(@__DIR__, "..", "mma")
const OUT_DIR = joinpath(@__DIR__, "..", "jl", "generated")

const FREE_LIST = [(1, -2), (1, -1), (2, -2), (2, -1), (3, -2), (4, -2)]
const FREE_KWARG = Dict((1, -2) => "c1m2", (1, -1) => "c1m1", (2, -2) => "c2m2",
                         (2, -1) => "c2m1", (3, -2) => "c3m2", (4, -2) => "c4m2")

signed_str(n::Integer) = n >= 0 ? string(n) : "m$(-n)"

function split_statements(text::String)
    lines = split(text, '\n')
    stmts = String[]
    cur = String[]
    for line in lines
        if isempty(strip(line))
            continue
        end
        if !startswith(line, " ") && !startswith(line, "\t")
            if !isempty(cur)
                push!(stmts, replace(join(cur, ""), r"\s+" => ""))
            end
            cur = String[line]
        else
            push!(cur, line)
        end
    end
    if !isempty(cur)
        push!(stmts, replace(join(cur, ""), r"\s+" => ""))
    end
    return stmts
end

function apply_renames(s::String)
    for x in ("e", "o"), (n, tgt) in (("1", "ψ1"), ("2", "ψ2"))
        s = replace(s, "Derivative[1][\\[Psi]1l$(n)$(x)][r]" => tgt)
    end
    for x in ("e", "o"), (n, tgt) in (("1", "ϕ1"), ("2", "ϕ2"))
        s = replace(s, "\\[Psi]1l$(n)$(x)[r]" => tgt)
    end
    s = replace(s, "\\[Omega]1" => "ω1")
    s = replace(s, "\\[Omega]2" => "ω2")
    s = replace(s, "\\[Lambda]2l" => "λ2l")
    s = replace(s, "\\[Chi]1" => "χ1")
    s = replace(s, "\\[Chi]2" => "χ2")
    s = replace(s, "\\[Theta]" => "θ")
    return s
end

const RE_REGSOURCE  = r"^RegSource(EEE|EEO|EOE|EOO|OOE|OOO)=(.*)$"
const RE_M          = r"^M=(.*)$"
const RE_CLOG       = r"^c\[(\d+)\]\[log\]=(.*)$"
const RE_C          = r"^c\[(\d+),(-?\d+)\]=(.*)$"
const RE_JSOVERRIDE = r"^JS\[\{l1,(-?\d+)\},\{l2,(-?\d+)\},\{l,(-?\d+)\}\]=(.*)$"
const RE_ZQ         = r"^z([pm])Q\[(\d+),(\d+),(-?\d+)\]=(.*)$"
const RE_LAMBDA     = r"^λ2l=(.*)$"
const RE_PARITY     = r"^parity=(.*)$"

function classify(s::String, sector::String)
    if occursin(":=", s)
        return (:skip, (), "")
    end
    m = match(RE_REGSOURCE, s)
    if m !== nothing
        m.captures[1] == sector || error("RegSource sector mismatch: $(m.captures[1]) vs $sector")
        return (:regsource, (), String(m.captures[2]))
    end
    m = match(RE_M, s)
    if m !== nothing
        return (:M, (), String(m.captures[1]))
    end
    m = match(RE_CLOG, s)
    if m !== nothing
        return (:clog, (parse(Int, m.captures[1]),), String(m.captures[2]))
    end
    m = match(RE_C, s)
    if m !== nothing
        return (:c, (parse(Int, m.captures[1]), parse(Int, m.captures[2])), String(m.captures[3]))
    end
    m = match(RE_JSOVERRIDE, s)
    if m !== nothing
        return (:jsoverride, (parse(Int, m.captures[1]), parse(Int, m.captures[2]), parse(Int, m.captures[3])),
                String(m.captures[4]))
    end
    m = match(RE_ZQ, s)
    if m !== nothing
        return (:zq, (parse(Int, m.captures[2]), parse(Int, m.captures[3]), parse(Int, m.captures[4])),
                String(m.captures[5]))
    end
    m = match(RE_LAMBDA, s)
    if m !== nothing
        return (:lambda2l, (), String(m.captures[1]))
    end
    m = match(RE_PARITY, s)
    if m !== nothing
        return (:parity, (), String(m.captures[1]))
    end
    error("unrecognized statement form: $(s[1:min(end, 200)])")
end

function parse_rhs_fallback(rhs::String)
    chars = collect(rhs)
    n = length(chars)
    terms = Tuple{Char,String}[]
    depth = 0
    start = 1
    sign = '+'
    for i in 1:n
        c = chars[i]
        if c in ('(', '[', '{')
            depth += 1
        elseif c in (')', ']', '}')
            depth -= 1
        elseif depth == 0 && (c == '+' || c == '-') && i > 1
            prevc = chars[i-1]
            if !(prevc in ('*', '/', '^', '(', '[', '{', ',', '+', '-'))
                push!(terms, (sign, String(chars[start:i-1])))
                sign = c
                start = i + 1
            end
        end
    end
    push!(terms, (sign, String(chars[start:end])))
    e = nothing
    for (sgn, t) in terms
        te = Meta.parse(t)
        if e === nothing
            e = sgn == '+' ? te : Expr(:call, :-, te)
        else
            e = Expr(:call, sgn == '+' ? :(+) : :(-), e, te)
        end
    end
    return e
end

function parse_rhs(rhs::String)
    try
        return Meta.parse(rhs)
    catch err
        @warn "direct Meta.parse failed, falling back to top-level +/- split" exception = err
        return parse_rhs_fallback(rhs)
    end
end

# dependency / usage scanning on the raw (pre-translation) parsed Expr
function scan!(e, cdeps::Set{Any}, jsused::Set{NTuple{3,Int}})
    e isa Expr || return
    if e.head == :ref
        f = e.args[1]
        idx = e.args[2:end]
        if f isa Expr && f.head == :ref && length(f.args) == 2 && f.args[1] === :c && idx == Any[:log]
            push!(cdeps, (:clog, f.args[2]))
            return
        elseif f === :c
            i, j = idx
            push!(cdeps, (:c, i, j))
            return
        elseif f === :JS
            a = idx[1].args[2]; b = idx[2].args[2]; cc = idx[3].args[2]
            push!(jsused, (a, b, cc))
            return
        elseif f === :zpQ || f === :zmQ
            for a in idx
                scan!(a, cdeps, jsused)
            end
            return
        else
            for a in idx
                scan!(a, cdeps, jsused)
            end
            return
        end
    elseif e.head == :braces
        return
    else
        for a in e.args
            scan!(a, cdeps, jsused)
        end
    end
end

function brace_second(x)
    x isa Expr && x.head == :braces && length(x.args) == 2 || error("expected {sym, int}, got $x")
    return x.args[2]
end

function translate_expr(e)
    if e isa Expr
        if e.head == :ref
            return translate_ref(e)
        elseif e.head == :braces
            error("unexpected bare {...} node in translated expression: $e")
        else
            return Expr(e.head, map(translate_expr, e.args)...)
        end
    elseif e === :I
        return :im
    else
        return e
    end
end

function translate_ref(e::Expr)
    f = e.args[1]
    idx = e.args[2:end]
    if f isa Expr && f.head == :ref && length(f.args) == 2 && f.args[1] === :c && idx == Any[:log]
        i = f.args[2]
        return Symbol("clog_$(i)")
    elseif f === :Sqrt
        length(idx) == 1 || error("Sqrt with != 1 arg: $e")
        return Expr(:call, :sqrt, translate_expr(idx[1]))
    elseif f === :Cot
        return Expr(:call, :cot, translate_expr(idx[1]))
    elseif f === :Csc
        return Expr(:call, :csc, translate_expr(idx[1]))
    elseif f === :Cos
        return Expr(:call, :cos, translate_expr(idx[1]))
    elseif f === :Log
        return Expr(:call, :log, translate_expr(idx[1]))
    elseif f === :c
        i, j = idx
        return Symbol("c_$(i)_$(signed_str(j))")
    elseif f === :JS
        a = brace_second(idx[1]); b = brace_second(idx[2]); cc = brace_second(idx[3])
        return Symbol("JS_$(signed_str(a))_$(signed_str(b))_$(signed_str(cc))")
    elseif f === :zpQ || f === :zmQ
        i, j, k = idx
        return Symbol("zq_$(i)_$(j)_$(signed_str(k))")
    else
        error("unrecognized ref head $(f) in $e")
    end
end

function toposort(nodes::Vector, edges::Dict)
    order = Any[]
    state = Dict{Any,Symbol}(n => :white for n in nodes)
    function visit(n)
        state[n] == :black && return
        state[n] == :gray && error("dependency cycle detected at $n")
        state[n] = :gray
        for d in get(edges, n, ())
            if haskey(state, d)
                visit(d)
            end
        end
        state[n] = :black
        push!(order, n)
    end
    for n in nodes
        visit(n)
    end
    return order
end

# CSE (common subexpression elimination) on a translated Julia Expr
function count_size(e)
    e isa Expr ? 1 + sum(count_size, e.args; init=0) : 1
end

function cse_hashes!(e, memo::IdDict{Any,UInt64})
    haskey(memo, e) && return memo[e]
    h = if e isa Expr
        hh = hash(e.head)
        for a in e.args
            hh = hash(cse_hashes!(a, memo), hh)
        end
        hh
    else
        hash(e)
    end
    memo[e] = h
    return h
end

function cse(e; prefix="_t", counter::Ref{Int}=Ref(0))
    memo = IdDict{Any,UInt64}()
    counts = Dict{UInt64,Int}()
    function count_calls!(x)
        if x isa Expr
            for a in x.args
                count_calls!(a)
            end
            if x.head == :call
                h = cse_hashes!(x, memo)
                counts[h] = get(counts, h, 0) + 1
            end
        end
    end
    count_calls!(e)

    bindings = Pair{Symbol,Any}[]
    group_var = Dict{UInt64,Symbol}()
    group_repr = Dict{UInt64,Any}()

    function rebuild(x)
        if x isa Expr
            newargs = [rebuild(a) for a in x.args]
            newe = Expr(x.head, newargs...)
            if x.head == :call
                h = cse_hashes!(x, memo)
                if get(counts, h, 0) >= 2
                    if haskey(group_var, h) && group_repr[h] == newe
                        return group_var[h]
                    elseif !haskey(group_var, h)
                        counter[] += 1
                        name = Symbol("$(prefix)$(counter[])")
                        group_var[h] = name
                        group_repr[h] = newe
                        push!(bindings, name => newe)
                        return name
                    end
                    # hash collision with structurally different expr: don't hoist
                end
            end
            return newe
        else
            return x
        end
    end
    final = rebuild(e)
    return bindings, final
end

function process_file(sector::String)
    path = joinpath(MMA_DIR, "Coefficients$(sector).wl")
    text = read(path, String)
    stmts = split_statements(text)

    regsource_raw = nothing
    defs = Dict{Tuple,Any}()
    js_override = Set{NTuple{3,Int}}()

    for raw in stmts
        s = apply_renames(raw)
        kind, args, rhs_str = classify(s, sector)
        kind == :skip && continue
        rhs_expr = parse_rhs(rhs_str)
        if kind == :regsource
            regsource_raw = rhs_expr
        elseif kind == :M
            rhs_expr == 1 || error("expected M=1, got $rhs_expr")
        elseif kind == :c
            defs[(:c, args...)] = rhs_expr
        elseif kind == :clog
            defs[(:clog, args...)] = rhs_expr
        elseif kind == :jsoverride
            push!(js_override, args)
        elseif kind == :zq
            defs[(:zq, args...)] = rhs_expr
        elseif kind == :lambda2l
            defs[(:lambda2l,)] = rhs_expr
        elseif kind == :parity
            defs[(:parity,)] = rhs_expr
        end
    end
    regsource_raw === nothing && error("no RegSource$(sector) statement found")
    return regsource_raw, defs, js_override
end

struct GenResult
    code::String
    stats::Dict{String,Any}
end

function generate_sector(sector::String; cse_enabled::Bool=true)
    regsource_raw, defs, js_override = process_file(sector)

    cdeps = Set{Any}()
    jsused = Set{NTuple{3,Int}}()
    scan!(regsource_raw, cdeps, jsused)
    for (k, v) in defs
        scan!(v, cdeps, jsused)
    end
    for t in js_override
        push!(jsused, t)
    end

    used_c = Set{Tuple{Int,Int}}((k[2], k[3]) for k in cdeps if k[1] == :c)
    defined_c = Set{Tuple{Int,Int}}((k[2], k[3]) for k in keys(defs) if k[1] == :c)
    free_c = setdiff(used_c, defined_c)
    expected_free = Set(FREE_LIST)
    if free_c != expected_free
        error("sector $sector: unexpected free coefficient set $(sort(collect(free_c))); " *
              "expected $(sort(collect(expected_free)))")
    end

    used_clog = Set{Int}(k[2] for k in cdeps if k[1] == :clog)
    defined_clog = Set{Int}(k[2] for k in keys(defs) if k[1] == :clog)

    all_c_nodes = union(defined_c, free_c)
    all_clog_nodes = union(defined_clog, used_clog)

    nodes = Any[]
    for (i, j) in all_c_nodes
        push!(nodes, (:c, i, j))
    end
    for i in all_clog_nodes
        push!(nodes, (:clog, i))
    end

    edges = Dict{Any,Set{Any}}()
    for (i, j) in defined_c
        d = Set{Any}(); js = Set{NTuple{3,Int}}()
        scan!(defs[(:c, i, j)], d, js)
        edges[(:c, i, j)] = Set{Any}((k[1] == :c ? (:c, k[2], k[3]) : (:clog, k[2])) for k in d)
    end
    for i in defined_clog
        d = Set{Any}(); js = Set{NTuple{3,Int}}()
        scan!(defs[(:clog, i)], d, js)
        edges[(:clog, i)] = Set{Any}((k[1] == :c ? (:c, k[2], k[3]) : (:clog, k[2])) for k in d)
    end
    for (i, j) in free_c
        edges[(:c, i, j)] = Set{Any}()
    end
    for i in setdiff(used_clog, defined_clog)
        edges[(:clog, i)] = Set{Any}()
    end

    order = toposort(nodes, edges)

    sector_lc = lowercase(sector)
    lines_code = String[]
    push!(lines_code, "KK = Float64(wigner3j(l1, l2, l, m1, m2, -(m1+m2))) * sqrt((2l+1)*(2l1+1)*(2l2+1)/(4*pi))")

    sorted_js = sort(collect(jsused))
    for (a, b, c) in sorted_js
        name = "JS_$(signed_str(a))_$(signed_str(b))_$(signed_str(c))"
        if (a, b, c) in js_override
            push!(lines_code, "$name = 0.0")
        else
            push!(lines_code, "$name = Float64(wigner3j(l1, l2, l, $a, $b, $c))")
        end
    end

    push!(lines_code, "l1 = Float64(l1)")
    push!(lines_code, "l2 = Float64(l2)")
    push!(lines_code, "l = Float64(l)")
    push!(lines_code, "M = 1.0")

    for (i, j) in FREE_LIST
        push!(lines_code, "c_$(i)_$(signed_str(j)) = $(FREE_KWARG[(i, j)])")
    end

    # c[i,j]/clog[i] RHS's can have hundreds of thousands of nodes with heavy internal
    # repetition (brutal for Julia's compiler) — CSE each one, sharing one temp-var counter
    # across the function body so names don't collide.
    coeff_cse_counter = Ref(0)
    coeff_size_before = 0
    coeff_size_after = 0
    for node in order
        if node[1] == :c
            i, j = node[2], node[3]
            (i, j) in free_c && continue
            rhs_expr = translate_expr(defs[(:c, i, j)])
            coeff_size_before += count_size(rhs_expr)
            if cse_enabled
                cbindings, final_rhs = cse(rhs_expr; prefix="_ct", counter=coeff_cse_counter)
            else
                cbindings, final_rhs = Pair{Symbol,Any}[], rhs_expr
            end
            coeff_size_after += count_size(final_rhs) + sum(count_size(b.second) for b in cbindings; init=0)
            for (nm, rhs) in cbindings
                push!(lines_code, "$(nm) = $(string(rhs))")
            end
            push!(lines_code, "c_$(i)_$(signed_str(j)) = $(string(final_rhs))")
        else
            i = node[2]
            if haskey(defs, (:clog, i))
                rhs_expr = translate_expr(defs[(:clog, i)])
                coeff_size_before += count_size(rhs_expr)
                if cse_enabled
                    cbindings, final_rhs = cse(rhs_expr; prefix="_ct", counter=coeff_cse_counter)
                else
                    cbindings, final_rhs = Pair{Symbol,Any}[], rhs_expr
                end
                coeff_size_after += count_size(final_rhs) + sum(count_size(b.second) for b in cbindings; init=0)
                for (nm, rhs) in cbindings
                    push!(lines_code, "$(nm) = $(string(rhs))")
                end
                push!(lines_code, "clog_$(i) = $(string(final_rhs))")
            else
                push!(lines_code, "clog_$(i) = 0.0")
            end
        end
    end

    if haskey(defs, (:lambda2l,))
        push!(lines_code, "λ2l = $(string(translate_expr(defs[(:lambda2l,)])))")
    end
    if haskey(defs, (:parity,))
        push!(lines_code, "parity = $(string(translate_expr(defs[(:parity,)])))")
    end

    zq_entries = sort([(k[2], k[3], k[4]) for k in keys(defs) if k[1] == :zq])
    zq_lines = String[]
    for (i, j, k) in zq_entries
        rhs_expr = translate_expr(defs[(:zq, i, j, k)])
        push!(zq_lines, "($i, $j, $k) => $(string(rhs_expr))")
    end
    zq_closure = "(χ1, χ2) -> Dict(" * join(zq_lines, ", ") * ")"

    field_list = String[]
    push!(field_list, "l1=l1, l2=l2, l=l, ω1=ω1, ω2=ω2, KK=KK")
    for (a, b, c) in sorted_js
        nm = "JS_$(signed_str(a))_$(signed_str(b))_$(signed_str(c))"
        push!(field_list, "$nm=$nm")
    end
    for (i, j) in sort(collect(all_c_nodes))
        nm = "c_$(i)_$(signed_str(j))"
        push!(field_list, "$nm=$nm")
    end
    for i in sort(collect(all_clog_nodes))
        nm = "clog_$(i)"
        push!(field_list, "$nm=$nm")
    end
    push!(field_list, "zq=$(zq_closure)")

    coeff_fn = """
function coefficients_$(sector_lc)(l1::Int, l2::Int, l::Int, m1::Int, m2::Int, ω1, ω2;
                          c1m2=0, c1m1=0, c2m2=0, c2m1=0, c3m2=0, c4m2=0)
    $(join(lines_code, "\n    "))
    return ($(join(field_list, ",\n            ")))
end
"""

    regsource_translated = translate_expr(regsource_raw)
    size_before = count_size(regsource_translated)
    if cse_enabled
        bindings, final_expr = cse(regsource_translated)
    else
        bindings, final_expr = Pair{Symbol,Any}[], regsource_translated
    end
    size_after = count_size(final_expr) + sum(count_size(b.second) for b in bindings; init=0)

    unpack_lines = String[]
    push!(unpack_lines, "l1 = cf.l1"); push!(unpack_lines, "l2 = cf.l2"); push!(unpack_lines, "l = cf.l")
    push!(unpack_lines, "ω1 = cf.ω1"); push!(unpack_lines, "ω2 = cf.ω2"); push!(unpack_lines, "KK = cf.KK")
    for (a, b, c) in sorted_js
        nm = "JS_$(signed_str(a))_$(signed_str(b))_$(signed_str(c))"
        push!(unpack_lines, "$nm = cf.$nm")
    end
    for (i, j) in sort(collect(all_c_nodes))
        nm = "c_$(i)_$(signed_str(j))"
        push!(unpack_lines, "$nm = cf.$nm")
    end
    for i in sort(collect(all_clog_nodes))
        nm = "clog_$(i)"
        push!(unpack_lines, "$nm = cf.$nm")
    end

    cse_lines = [ "$(nm) = $(string(rhs))" for (nm, rhs) in bindings ]

    source_fn = """
function source_$(sector_lc)(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    $(join(unpack_lines, "\n    "))
    M = 1.0
    θ = pi/4
    $(join(cse_lines, "\n    "))
    return $(string(final_expr))
end
"""

    header = "using WignerSymbols\n\n"
    code = header * coeff_fn * "\n" * source_fn

    stats = Dict{String,Any}(
        "n_statements" => length(sorted_js) + length(all_c_nodes) + length(all_clog_nodes),
        "n_js" => length(sorted_js),
        "n_c" => length(all_c_nodes),
        "n_free_c" => length(free_c),
        "n_clog" => length(all_clog_nodes),
        "n_zq" => length(zq_entries),
        "cse_bindings" => length(bindings),
        "regsource_size_before_cse" => size_before,
        "regsource_size_after_cse" => size_after,
        "coeff_cse_bindings" => coeff_cse_counter[],
        "coeff_size_before_cse" => coeff_size_before,
        "coeff_size_after_cse" => coeff_size_after,
    )
    return GenResult(code, stats)
end

function main(; cse_enabled::Bool=true, sectors=SECTORS)
    mkpath(OUT_DIR)
    for sector in sectors
        t0 = time()
        result = generate_sector(sector; cse_enabled=cse_enabled)
        dt = time() - t0
        outpath = joinpath(OUT_DIR, "Source$(sector).jl")
        write(outpath, result.code)
        st = result.stats
        shrink = st["regsource_size_before_cse"] == 0 ? 0.0 :
                 100 * (1 - st["regsource_size_after_cse"] / st["regsource_size_before_cse"])
        coeff_shrink = st["coeff_size_before_cse"] == 0 ? 0.0 :
                 100 * (1 - st["coeff_size_after_cse"] / st["coeff_size_before_cse"])
        println("[$sector] wrote $outpath in $(round(dt, digits=2))s -- ",
                "JS=$(st["n_js"]) c=$(st["n_c"]) (free=$(st["n_free_c"])) clog=$(st["n_clog"]) zq=$(st["n_zq"]) ",
                "-- source CSE bindings=$(st["cse_bindings"]), tree nodes $(st["regsource_size_before_cse"]) -> ",
                "$(st["regsource_size_after_cse"]) ($(round(shrink, digits=1))% smaller)",
                " -- coeff CSE bindings=$(st["coeff_cse_bindings"]), tree nodes $(st["coeff_size_before_cse"]) -> ",
                "$(st["coeff_size_after_cse"]) ($(round(coeff_shrink, digits=1))% smaller)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    cse_flag = !("--no-cse" in ARGS)
    main(; cse_enabled=cse_flag)
end
