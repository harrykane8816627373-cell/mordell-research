import Mordell.Kernel
import Mathlib.Logic.Equiv.Fin.Basic

/-! The mass-eleven four-triangle model and the twelve-point variational reduction. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

def omega : ℂ := ⟨-1 / 2, Real.sqrt 3 / 2⟩

theorem omega_sq : omega ^ 2 = conj omega := by
  have hs := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 3)
  apply Complex.ext <;>
    simp only [pow_two, Complex.mul_re, Complex.mul_im, Complex.conj_re, Complex.conj_im, omega]
  · nlinarith
  · ring

theorem normSq_omega : Complex.normSq omega = 1 := by
  simp only [Complex.normSq_apply, omega]
  nlinarith [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 3)]

theorem omega_cube : omega ^ 3 = 1 := by
  rw [show (3 : ℕ) = 2 + 1 by decide, pow_add, pow_one, omega_sq,
    mul_comm, Complex.mul_conj, normSq_omega, Complex.ofReal_one]

theorem omega_triangle : 1 + omega + omega ^ 2 = 0 := by
  rw [omega_sq]
  apply Complex.ext <;>
    simp only [omega, Complex.add_re, Complex.add_im,
      Complex.one_re, Complex.one_im, Complex.conj_re, Complex.conj_im,
      Complex.zero_re, Complex.zero_im] <;> ring

private theorem omega_powers_injective : Function.Injective (fun j : Fin 3 => omega ^ j.val) := by
  have hi : 0 < Real.sqrt 3 / 2 := by positivity
  intro i j h
  have he := congrArg Complex.im h
  fin_cases i <;> fin_cases j
  all_goals first | rfl | skip
  all_goals simp only [pow_zero, pow_one, omega_sq, Complex.conj_im, Complex.one_im] at he
  all_goals { simp only [omega] at he; exfalso; linarith }

private theorem omega_power_cube (j : ℕ) : (omega ^ j) ^ 3 = 1 := by
  rw [← pow_mul, mul_comm j, pow_mul, omega_cube, one_pow]

structure TriangleParameters where
  inner : ℝ
  middle : ℝ
  outerRe : ℝ
  outerIm : ℝ

namespace TriangleParameters

def outer (p : TriangleParameters) : ℂ := ⟨p.outerRe, p.outerIm⟩

def seeds (p : TriangleParameters) : Fin 4 → ℂ :=
  ![(p.inner : ℂ), -(p.middle : ℂ), p.outer, conj p.outer]

def RadiusMass (p : TriangleParameters) : Prop :=
  3 * p.inner ^ 2 + 3 * p.middle ^ 2 + 6 * (p.outerRe ^ 2 + p.outerIm ^ 2) = 11

end TriangleParameters

def triangleOrbits (a : Fin 4 → ℂ) : Configuration 12 :=
  fun k =>
    let ij : Fin 4 × Fin 3 := finProdFinEquiv.symm k
    a ij.1 * omega ^ ij.2.val

def candidateTwelve (p : TriangleParameters) : Configuration 12 :=
  triangleOrbits p.seeds

theorem mass_triangleOrbits (a : Fin 4 → ℂ) :
    mass (triangleOrbits a) = 3 * ∑ j, Complex.normSq (a j) := by
  classical
  have hpow (j : Fin 3) : Complex.normSq (omega ^ j.val) = 1 := by
    rw [map_pow, normSq_omega, one_pow]
  calc
    mass (triangleOrbits a) =
        ∑ ij : Fin 4 × Fin 3, Complex.normSq (triangleOrbits a (finProdFinEquiv ij)) :=
      (Equiv.sum_comp (finProdFinEquiv (m := 4) (n := 3))
        (fun k : Fin 12 => Complex.normSq (triangleOrbits a k))).symm
    _ = _ := by
      simp only [triangleOrbits, Equiv.symm_apply_apply, Fintype.sum_prod_type, Complex.normSq_mul,
        hpow, mul_one, Finset.sum_const,
        Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, ← Finset.mul_sum]
      norm_num

theorem mass_candidateTwelve (p : TriangleParameters) :
    mass (candidateTwelve p) =
      3 * p.inner ^ 2 + 3 * p.middle ^ 2 + 6 * (p.outerRe ^ 2 + p.outerIm ^ 2) := by
  rw [candidateTwelve, mass_triangleOrbits]
  simp only [TriangleParameters.seeds, Fin.sum_univ_succ, Matrix.cons_val_zero,
    Matrix.cons_val_succ, Finset.univ_eq_empty, Finset.sum_empty, Complex.normSq_neg,
    Complex.normSq_conj, Complex.normSq_ofReal, TriangleParameters.outer, Complex.normSq_mk]
  ring

theorem triangleOrbits_injective {a : Fin 4 → ℂ}
    (hnz : ∀ j, a j ≠ 0) (hc : Function.Injective (fun j => a j ^ 3)) :
    Function.Injective (triangleOrbits a) := by
  intro k l h
  let i : Fin 4 × Fin 3 := finProdFinEquiv.symm k
  let j : Fin 4 × Fin 3 := finProdFinEquiv.symm l
  change a i.1 * omega ^ i.2.val = a j.1 * omega ^ j.2.val at h
  have hp : (a i.1 * omega ^ i.2.val) ^ 3 = (a j.1 * omega ^ j.2.val) ^ 3 :=
    congrArg (fun z : ℂ => z ^ 3) h
  simp only [mul_pow, omega_power_cube, mul_one] at hp
  have hfirst : i.1 = j.1 := hc hp
  rw [hfirst] at h
  have hsecond := omega_powers_injective (mul_left_cancel₀ (hnz _) h)
  exact finProdFinEquiv.symm.injective (Prod.ext hfirst hsecond)

theorem triangle_first_moment (w : ℂ) :
    (∑ j : Fin 3, w * omega ^ j.val) = 0 := by
  rw [← Finset.mul_sum]
  suffices (∑ j : Fin 3, omega ^ j.val) = 0 by rw [this, mul_zero]
  simpa [Fin.sum_univ_succ, add_assoc] using omega_triangle

theorem triangle_third_moment (w : ℂ) :
    (∑ j : Fin 3, (w * omega ^ j.val) ^ 3) = 3 * w ^ 3 := by
  simp only [mul_pow, omega_power_cube, mul_one,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  norm_num

abbrev upScaleTwelve : ℂ := (Real.sqrt (12 / 11 : ℝ) : ℝ)

theorem normSq_upScaleTwelve : Complex.normSq upScaleTwelve = (12 / 11 : ℝ) := by
  rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt (by norm_num)]

theorem pairCount_twelve : pairCount 12 = 66 := by decide

theorem twelve_value_scaling (p : TriangleParameters) :
    distanceProduct (scale upScaleTwelve (candidateTwelve p)) =
      (12 / 11 : ℝ) ^ 66 * distanceProduct (candidateTwelve p) := by
  rw [distanceProduct_scale, pairCount_twelve, normSq_upScaleTwelve]

theorem sharp_twelve_of_log_certificate (p : TriangleParameters)
    (hm : p.RadiusMass) (hnz : ∀ j, p.seeds j ≠ 0)
    (hc : Function.Injective (fun j => p.seeds j ^ 3))
    (hcert : LogCertificate 11 (candidateTwelve p)) :
    IsMaximizer 12 (scale upScaleTwelve (candidateTwelve p)) ∧
      (∀ z : Configuration 12, IsMaximizer 12 z →
        SameOrbit z (scale upScaleTwelve (candidateTwelve p))) := by
  have hm' : mass (candidateTwelve p) = 11 := (mass_candidateTwelve p).trans hm
  have hs : upScaleTwelve ≠ 0 := Complex.normSq_pos.mp (by rw [normSq_upScaleTwelve]; norm_num)
  have h := sharp_of_log_certificate hm'
    (triangleOrbits_injective hnz hc) hs hcert
  rw [normSq_upScaleTwelve, show (12 / 11 : ℝ) * 11 = 12 by norm_num] at h
  exact h

theorem sharp_twelve_of_pair_majorants (p : TriangleParameters)
    (hm : p.RadiusMass) (hnz : ∀ j, p.seeds j ≠ 0)
    (hc : Function.Injective (fun j => p.seeds j ^ 3))
    (h : ∀ z : Configuration 12, IsMaximizer 11 z → Function.Injective z →
      Nonempty (PairMajorant z (candidateTwelve p))) :
    IsMaximizer 12 (scale upScaleTwelve (candidateTwelve p)) ∧
      (∀ z : Configuration 12, IsMaximizer 12 z →
        SameOrbit z (scale upScaleTwelve (candidateTwelve p))) := by
  have hlog := log_certificate_of_pair_majorants (by decide : 1 < 12)
    (candidateTwelve p) (by simpa only [Nat.cast_ofNat, show (12 : ℝ) - 1 = 11 by norm_num] using h)
  exact sharp_twelve_of_log_certificate p hm hnz hc
    (by simpa only [Nat.cast_ofNat, show (12 : ℝ) - 1 = 11 by norm_num] using hlog)

end Mordell
