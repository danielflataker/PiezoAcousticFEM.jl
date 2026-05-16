"""
    DirichletDofs(indices, values)

Essential boundary conditions on global DOFs. `indices` are 1-based indices in
the unreduced system, and `values` are the prescribed DOF values.
"""
struct DirichletDofs{T}
    indices::Vector{Int}
    values::Vector{T}

    function DirichletDofs(indices, values)
        i = collect(Int, indices)
        g = collect(values)
        length(i) == length(g) ||
            throw(DimensionMismatch("Dirichlet indices and values must have the same length"))
        length(unique(i)) == length(i) ||
            throw(ArgumentError("Dirichlet indices contain duplicates"))

        return new{eltype(g)}(i, g)
    end
end


"""
    DirichletReduction

Reduced linear system produced by eliminating constrained DOFs from `A*x=b`.
`A` and `b` contain the free-free block and corrected RHS. Use
[`reconstruct_solution`](@ref) to expand a free solution back to the full
unreduced DOF vector.
"""
struct DirichletReduction{T,A,B}
    A::A
    b::B
    free::Vector{Int}
    constrained::Vector{Int}
    values::Vector{T}
    full_length::Int
end


"""
    apply_dirichlet(A, b, constraints)

Eliminate Dirichlet DOFs from `A*x=b`. If `x_c=g`, the returned reduced system
solves

```math
A_{ff} x_f = b_f - A_{fc} g.
```
"""
function apply_dirichlet(A::AbstractMatrix, b::AbstractVector, constraints::DirichletDofs)
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("A must be square"))
    length(b) == size(A, 1) || throw(DimensionMismatch("b length must match A"))

    n = length(b)
    c = constraints.indices
    all(1 .<= c .<= n) || throw(ArgumentError("Dirichlet index outside 1:$n"))

    free = setdiff(collect(1:n), c)
    A_ff = A[free, free]
    A_fc = A[free, c]
    b_f = b[free] - A_fc * constraints.values

    return DirichletReduction(A_ff, b_f, free, c, constraints.values, n)
end


"""
    reconstruct_solution(reduction, x_free)

Expand a reduced free-DOF solution back to the full DOF vector, inserting the
prescribed Dirichlet values at constrained indices.
"""
function reconstruct_solution(reduction::DirichletReduction, x_free::AbstractVector)
    length(x_free) == length(reduction.free) ||
        throw(DimensionMismatch("x_free length must match reduced free DOFs"))

    T = promote_type(eltype(x_free), eltype(reduction.values))
    x = zeros(T, reduction.full_length)
    x[reduction.free] = x_free
    x[reduction.constrained] = reduction.values

    return x
end
