import Mordell.Stationarity

/-! Pair invariants and the real radial force row. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

def realDot (z w : ℂ) : ℝ := (z * conj w).re

theorem distanceSq_invariants {n : ℕ} (z : Configuration n) (i j : Fin n) :
    distanceSq z i j =
      Complex.normSq (z i) + Complex.normSq (z j) - 2 * realDot (z i) (z j) :=
  Complex.normSq_sub _ _

theorem gram_inequality (z w : ℂ) :
    realDot z w ^ 2 ≤ Complex.normSq z * Complex.normSq w := by
  have hs := Complex.normSq_mul z (conj w)
  rw [Complex.normSq_conj] at hs
  rw [← hs]
  unfold realDot
  rw [Complex.normSq_apply]
  nlinarith [sq_nonneg (z * conj w).im]

private theorem real_force_term (z w : ℂ) :
    (z * (z - w)⁻¹).re = (Complex.normSq z - realDot z w) / Complex.normSq (z - w) := by
  rw [← div_eq_mul_inv, Complex.div_re, ← add_div]
  congr 1
  simp only [Complex.normSq_apply, realDot, Complex.mul_re, Complex.conj_re,
    Complex.conj_im, Complex.sub_re, Complex.sub_im]
  ring

theorem radial_force_sum {n : ℕ} {z : Configuration n} {s : ℝ}
    (hz : Stationary (s : ℂ) z) (i : Fin n) :
    (∑ j, if i = j then 0 else
      (Complex.normSq (z i) - realDot (z i) (z j)) / distanceSq z i j) =
      s * Complex.normSq (z i) := by
  have hf := hz i
  have he : z i * force z i = (s : ℂ) * Complex.normSq (z i) := by
    rw [hf, ← Complex.mul_conj]
    ring
  have hreal : (z i * force z i).re = s * Complex.normSq (z i) := by
    simpa only [← Complex.ofReal_mul, Complex.ofReal_re] using congrArg Complex.re he
  rw [force, Finset.mul_sum] at hreal
  have hsum (f : Fin n → ℂ) : (∑ j, f j).re = ∑ j, (f j).re :=
    map_sum Complex.reAddGroupHom f Finset.univ
  rw [hsum] at hreal
  rw [← hreal]
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases h : i = j
  · simp [h]
  · simp only [if_neg h, distanceSq]
    exact (real_force_term (z i) (z j)).symm

end Mordell
