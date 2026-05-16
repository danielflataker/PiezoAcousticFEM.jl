"""
    KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)

Assembled K-form-system for en piezoelektrisk struktur i vakuum,

```text
-ω² [Muu  0] [u] + [Kuu  Kuϕ] [u] = [ F]
    [ 0   0] [ϕ]   [Kϕu  Kϕϕ] [ϕ]   [-Q]
```

Feltnavnene bruker en enkel variant av Kocbachs notasjon, men unngår
Unicode-subskript i public API.
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
