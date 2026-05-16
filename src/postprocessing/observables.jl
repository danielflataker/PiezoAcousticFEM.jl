"""
    AbstractObservable

Typed request for a scalar or vector quantity derived from an analysis result.
"""
abstract type AbstractObservable end

"""
    AdmittanceObservable()

Extract the driven-electrode admittance from a harmonic voltage result.
"""
struct AdmittanceObservable <: AbstractObservable end

"""
    ChargeObservable()

Extract the driven-electrode charge from a harmonic voltage result.
"""
struct ChargeObservable <: AbstractObservable end

"""
    CurrentObservable()

Extract the driven-electrode current from a harmonic voltage result.
"""
struct CurrentObservable <: AbstractObservable end

"""
    ModalFrequenciesObservable()

Extract modal frequencies in Hz from a short-circuit modal result.
"""
struct ModalFrequenciesObservable <: AbstractObservable end

"""
    evaluate(observable, result)

Evaluate an observable against an analysis result.
"""
function evaluate end

evaluate(::AdmittanceObservable, result::HarmonicVoltageResult) = result.solution.admittance
evaluate(::ChargeObservable, result::HarmonicVoltageResult) = result.solution.charge
evaluate(::CurrentObservable, result::HarmonicVoltageResult) = result.solution.current
evaluate(::ModalFrequenciesObservable, result::ShortCircuitModalResult) = result.frequencies
