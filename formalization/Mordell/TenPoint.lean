import Mordell.RootBox
import Mordell.Kernel

/-! The ten-point sharp bound reduced to its mass-nine logarithmic certificate. -/

noncomputable section

namespace Mordell

abbrev upScaleTen : ℂ := (Real.sqrt (10 / 9 : ℝ) : ℝ)

theorem normSq_upScaleTen : Complex.normSq upScaleTen = (10 / 9 : ℝ) := by
  rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt (by norm_num)]

theorem pairCount_ten : pairCount 10 = 45 := by decide

theorem ten_value_scaling (p : Parameters) :
    distanceProduct (scale upScaleTen (candidate p)) =
      (10 / 9 : ℝ) ^ 45 * distanceProduct (candidate p) := by
  rw [distanceProduct_scale, pairCount_ten, normSq_upScaleTen]

theorem sharp_ten_of_log_certificate (p : Parameters) (hp : p.IsRoot)
    (hcert : LogCertificate 9 (candidate p)) :
    IsMaximizer 10 (scale upScaleTen (candidate p)) ∧
      (∀ z : Configuration 10, IsMaximizer 10 z → SameOrbit z (scale upScaleTen (candidate p))) := by
  have hs : upScaleTen ≠ 0 := Complex.normSq_pos.mp (by rw [normSq_upScaleTen]; norm_num)
  have h := sharp_of_log_certificate (mass_candidate_of_root hp)
    (candidate_injective hp.1.admissible) hs hcert
  rw [normSq_upScaleTen, show (10 / 9 : ℝ) * 9 = 10 by norm_num] at h
  exact h

def exactTenConfiguration : Configuration 10 :=
  scale upScaleTen (candidate isolatedParameters)

def exactTenValue : ℝ := distanceProduct exactTenConfiguration

theorem exactTenConfiguration_mass : mass exactTenConfiguration = 10 := by
  rw [exactTenConfiguration, mass_scale, normSq_upScaleTen,
    mass_candidate_of_root isolatedParameters_isRoot]
  norm_num

/-- The isolated algebraic configuration gives an attained lower bound
    on the ten-point maximum, independently of the pair cover. -/
theorem exactTenValue_le_of_maximizer {z : Configuration 10} (hz : IsMaximizer 10 z) :
    exactTenValue ≤ distanceProduct z :=
  hz.2 exactTenConfiguration exactTenConfiguration_mass

theorem sharp_ten_of_pair_majorants
    (h : ∀ z : Configuration 10, IsMaximizer 9 z → Function.Injective z →
      Nonempty (PairMajorant z (candidate isolatedParameters))) :
    IsMaximizer 10 exactTenConfiguration ∧
      (∀ z : Configuration 10, IsMaximizer 10 z → SameOrbit z exactTenConfiguration) := by
  have hc := log_certificate_of_pair_majorants (by decide : 1 < 10)
    (candidate isolatedParameters) (by simpa only [Nat.cast_ofNat, show (10 : ℝ) - 1 = 9 by norm_num] using h)
  exact sharp_ten_of_log_certificate isolatedParameters isolatedParameters_isRoot
    (by simpa only [Nat.cast_ofNat, show (10 : ℝ) - 1 = 9 by norm_num] using hc)

end Mordell
