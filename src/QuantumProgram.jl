using MacroTools: postwalk, prewalk, rmlines, @capture

const QProg = Expr

QEmpty = rmlines(quote end)

M_counts = [0]

all_symbols(x::Symbol) = [x]
all_symbols(e::Expr) = e.head == :tuple ? length(e.args) == 0 ? Symbol[] : vcat(map(all_symbols, e.args)...) : nothing

macro qprog(name, args, e::Expr)
    args = all_symbols(args)
    e = postwalk(
        x -> @capture(x, M(b_)) ? :(M($(b), Expr(:string, $("M_$(name)_"), $(b), $("_$(M_counts[1]+=1)")))) : x,
        e
    )
    body = Expr(:quote, postwalk(x -> x in args ? Expr(:$, x) : x, e))
    fn = quote
        function $(esc(name))($(args...))
            $body
        end
    end
    fn = prewalk(rmlines, fn)
    return fn
end