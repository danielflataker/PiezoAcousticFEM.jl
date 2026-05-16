"""
    KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)

Assembled K-form system for a piezoelectric structure in vacuum,

```text
-ω² [Muu  0] [u] + [Kuu  Kuϕ] [u] = [ F]
    [ 0   0] [ϕ]   [Kϕu  Kϕϕ] [ϕ]   [-Q]
```

The field names use a simple variant of Kocbach's notation, while avoiding
Unicode subscripts in the public API.
"""
struct KFormSystem{
    KUU<:AbstractMatrix,
    KUP<:AbstractMatrix,
    KPU<:AbstractMatrix,
    KPP<:AbstractMatrix,
    MUU<:AbstractMatrix,
}
    Kuu::KUU
    Kuϕ::KUP
    Kϕu::KPU
    Kϕϕ::KPP
    Muu::MUU
end
