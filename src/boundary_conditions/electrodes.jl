"""
    FacetElectrode(facetset)

Explicit electrode object pointing to a Ferrite facetset or facetset name.
String names therefore stay in the adapter layer, not as the core API.
"""
abstract type AbstractElectrode end

struct FacetElectrode{F} <: AbstractElectrode
    facetset::F
end


"""
    TwoTerminalElectrodes(signal, reference)

Two-terminal electrode setup. The analysis decides whether `signal` is driven,
short-circuited, floating, etc.; `reference` is the reference electrode.
"""
struct TwoTerminalElectrodes{D<:AbstractElectrode,G<:AbstractElectrode}
    signal::D
    reference::G
end
