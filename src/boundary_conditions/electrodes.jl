"""
    FacetElectrode(facetset)

Eksplisitt elektrodeobjekt som peker på et Ferrite facetset eller et
facetset-navn. Strengnavn holdes dermed i adapterlaget, ikke som kjerne-API.
"""
abstract type AbstractElectrode end

struct FacetElectrode{F} <: AbstractElectrode
    facetset::F
end


"""
    TwoTerminalElectrodes(signal, reference)

To-terminal elektrodeoppsett. Analysen bestemmer om `signal` er drevet,
kortsluttet, flytende osv.; `reference` er referanseelektroden.
"""
struct TwoTerminalElectrodes{D<:AbstractElectrode,G<:AbstractElectrode}
    signal::D
    reference::G
end
