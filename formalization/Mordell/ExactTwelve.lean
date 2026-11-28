import Mordell.TwelveSecant
import Mordell.TwelvePoint

/-! The archived twelve-point root supplies an attained lower bound. -/

noncomputable section

namespace Mordell

def exactTwelveConfiguration : Configuration 12 := gaugeConfiguration TwelveRootBox.isolatedCoordinates

def exactTwelveValue : ℝ := distanceProduct exactTwelveConfiguration

theorem exactTwelveConfiguration_mass : mass exactTwelveConfiguration = 12 := TwelveRootBox.isolated_mass

theorem exactTwelveConfiguration_injective : Function.Injective exactTwelveConfiguration :=
  TwelveRootBox.injective_on_box TwelveRootBox.isolatedCoordinates_spec.1

theorem exactTwelveValue_le_of_maximizer {z : Configuration 12} (hz : IsMaximizer 12 z) :
    exactTwelveValue ≤ distanceProduct z :=
  hz.2 exactTwelveConfiguration exactTwelveConfiguration_mass

def downScaleTwelve : ℂ := (Real.sqrt (11 / 12 : ℝ) : ℝ)

def massElevenTwelve : Configuration 12 := scale downScaleTwelve exactTwelveConfiguration

theorem normSq_downScaleTwelve : Complex.normSq downScaleTwelve = (11 / 12 : ℝ) := by
  rw [downScaleTwelve, Complex.normSq_ofReal, ← sq, Real.sq_sqrt (by norm_num)]

theorem massElevenTwelve_mass : mass massElevenTwelve = 11 := by
  rw [massElevenTwelve, mass_scale, exactTwelveConfiguration_mass, normSq_downScaleTwelve]
  norm_num

theorem massElevenTwelve_stationary : Stationary 6 massElevenTwelve := by
  have h := TwelveRootBox.isolatedCoordinates_spec.2.scale downScaleTwelve
  rw [normSq_downScaleTwelve] at h
  norm_num at h
  exact h

end Mordell
