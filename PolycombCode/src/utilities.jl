flatten(a::Array{T,1}) where T = any(map(x->isa(x,Array),a)) ? flatten(vcat(map(flatten,a)...)) : a
flatten(a::Array{T}) where T = reshape(a,prod(size(a)))
flatten(a::Number)=a

function gentimestamp()
    return Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
end

sigmoid(x, a) = 1/(1+exp(-a*x))
sigmoid(x) = 1/(1+exp(-5*x))
