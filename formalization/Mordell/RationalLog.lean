import Mordell.LogBounds

/-! Rational bounds for logarithms from a finite odd-power series. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def logSeriesRat (terms : ℕ) (y : ℚ) : ℚ :=
  let t := (y - 1) / (y + 1)
  2 * (List.ofFn fun k : Fin terms => t ^ (2 * k.val + 1) / (2 * (k.val : ℚ) + 1)).sum

def logSeriesErrorRat (terms : ℕ) (y : ℚ) : ℚ :=
  let t := (y - 1) / (y + 1)
  2 * t ^ (2 * terms + 1) / (1 - t ^ 2)

private theorem logSeriesRat_cast (terms : ℕ) (y : ℚ) :
    (logSeriesRat terms y : ℝ) =
      2 * ∑ k ∈ Finset.range terms, (((y : ℝ) - 1) / ((y : ℝ) + 1)) ^ (2 * k + 1) / (2 * (k : ℝ) + 1) := by
  simp only [logSeriesRat, List.sum_ofFn, Rat.cast_mul, Rat.cast_ofNat, Rat.cast_sum, Rat.cast_div,
    Rat.cast_sub, Rat.cast_add, Rat.cast_one, Rat.cast_pow, Rat.cast_natCast]
  congr 1
  exact Fin.sum_univ_eq_sum_range (fun k : ℕ =>
    (((y : ℝ) - 1) / ((y : ℝ) + 1)) ^ (2 * k + 1) / (2 * (k : ℝ) + 1)) terms

theorem logSeriesRat_bounds {y : ℚ} (hy : 1 ≤ y) (terms : ℕ) :
    (logSeriesRat terms y : ℝ) ≤ Real.log (y : ℝ) ∧
      Real.log (y : ℝ) ≤ ((logSeriesRat terms y + logSeriesErrorRat terms y : ℚ) : ℝ) := by
  have hy' : (1 : ℝ) ≤ (y : ℝ) := by exact_mod_cast hy
  let t : ℝ := ((y : ℝ) - 1) / ((y : ℝ) + 1)
  have hp : 0 < (y : ℝ) + 1 := by linarith
  have ht : 0 ≤ t := div_nonneg (sub_nonneg.mpr hy') hp.le
  have ht' : t < 1 := (div_lt_one hp).mpr (by linarith)
  have he : (1 + t) / (1 - t) = (y : ℝ) := by
    have ht'' : 1 - t ≠ 0 := (sub_pos.mpr ht').ne'
    apply (div_eq_iff ht'').mpr
    dsimp [t]
    field_simp
    ring
  have hl := Real.sum_range_le_log_div ht ht' terms
  have hu := Real.log_div_le_sum_range_add ht ht' terms
  rw [he] at hl hu
  rw [Rat.cast_add, logSeriesRat_cast]
  simp only [logSeriesErrorRat, Rat.cast_div, Rat.cast_sub, Rat.cast_add, Rat.cast_ofNat, Rat.cast_one,
    Rat.cast_mul, Rat.cast_pow]
  change 2 * (∑ k ∈ Finset.range terms, t ^ (2 * k + 1) / (2 * (k : ℝ) + 1)) ≤ Real.log (y : ℝ) ∧
    Real.log (y : ℝ) ≤
      2 * (∑ k ∈ Finset.range terms, t ^ (2 * k + 1) / (2 * (k : ℝ) + 1)) +
        2 * t ^ (2 * terms + 1) / (1 - t ^ 2)
  rw [mul_div_assoc]
  constructor <;> linarith only [hl, hu]

def logUpperRat (terms : ℕ) (y : ℚ) : ℚ :=
  if y ≤ 1 then -logSeriesRat terms y⁻¹ else logSeriesRat terms y + logSeriesErrorRat terms y

theorem logUpperRat_sound {y : ℚ} (hy : 0 < y) (terms : ℕ) :
    Real.log (y : ℝ) ≤ (logUpperRat terms y : ℝ) := by
  by_cases h : y ≤ 1
  · have h' : 1 ≤ y⁻¹ := (one_le_inv₀ hy).2 h
    have hl := (logSeriesRat_bounds h' terms).1
    simp only [Rat.cast_inv, Real.log_inv] at hl
    simp only [logUpperRat, if_pos h, Rat.cast_neg]
    linarith
  · exact (logSeriesRat_bounds (le_of_lt (lt_of_not_ge h)) terms).2.trans_eq
      (by rw [logUpperRat, if_neg h])

end Mordell
