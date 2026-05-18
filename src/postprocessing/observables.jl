"""
    AbstractObservable

Typed request for a scalar or vector quantity derived from an analysis result.
"""
abstract type AbstractObservable end

"""
    AdmittanceObservable()

Extract the driven-terminal admittance from a harmonic voltage result, in
siemens. The value follows the result's harmonic time convention.
"""
struct AdmittanceObservable <: AbstractObservable end

"""
    ChargeObservable()

Extract the driven-terminal total charge from a harmonic voltage result, in
coulombs. The axisymmetric integration factor is already included.
"""
struct ChargeObservable <: AbstractObservable end

"""
    CurrentObservable()

Extract the driven-terminal current from a harmonic voltage result, in amperes.
The value follows the result's harmonic time convention.
"""
struct CurrentObservable <: AbstractObservable end

"""
    ModalFrequenciesObservable()

Extract modal frequencies from a short-circuit modal result, in hertz.
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
