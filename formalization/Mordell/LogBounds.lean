import Mordell.Basic
import Mathlib.Analysis.Convex.Deriv
import Mathlib.Analysis.SpecialFunctions.Log.Deriv

/-! The cubic upper bound for the logarithm after clearing the distance denominator. -/

noncomputable section

namespace Mordell

private def cubicGap (t : ℝ) : ℝ :=
  (3 / 2) * t ^ 2 - 2 * t + 1 / 2 + (t - 1) ^ 3 / 3 - t ^ 2 * Real.log t

private def cubicGapDeriv (t : ℝ) : ℝ := t ^ 2 - 1 - 2 * t * Real.log t

private theorem cubicGap_hasDerivAt {t : ℝ} (ht : 0 < t) :
    HasDerivAt cubicGap (cubicGapDeriv t) t := by
  have h := ((((((hasDerivAt_id t).pow 2).const_mul (3 / 2 : ℝ)).sub
    ((hasDerivAt_id t).const_mul (2 : ℝ))).add_const (1 / 2 : ℝ)).add
    ((((hasDerivAt_id t).sub_const (1 : ℝ)).pow 3).div_const (3 : ℝ))).sub
    (((hasDerivAt_id t).pow 2).mul (Real.hasDerivAt_log ht.ne'))
  convert h using 1
  simp only [cubicGapDeriv, Pi.pow_apply, id_eq]
  field_simp [ht.ne']
  ring

private theorem cubicGapDeriv_hasDerivAt {t : ℝ} (ht : 0 < t) :
    HasDerivAt cubicGapDeriv (2 * (t - 1 - Real.log t)) t := by
  have h := (((hasDerivAt_id t).pow 2).sub_const (1 : ℝ)).sub
    (((hasDerivAt_id t).const_mul (2 : ℝ)).mul (Real.hasDerivAt_log ht.ne'))
  convert h using 1
  simp only [id_eq]
  field_simp [ht.ne']
  ring

private theorem cubicGap_nonneg {t : ℝ} (ht : 0 < t) : 0 ≤ cubicGap t := by
  have hconv : ConvexOn ℝ (Set.Ioi 0) cubicGap := by
    apply convexOn_of_hasDerivWithinAt2_nonneg (convex_Ioi 0)
      (f' := cubicGapDeriv) (f'' := fun x => 2 * (x - 1 - Real.log x))
    · exact fun _ hx => (cubicGap_hasDerivAt hx).continuousAt.continuousWithinAt
    · exact fun _ hx => (cubicGap_hasDerivAt (interior_subset hx)).hasDerivWithinAt
    · exact fun _ hx => (cubicGapDeriv_hasDerivAt (interior_subset hx)).hasDerivWithinAt
    · intro x hx
      have h := Real.log_le_sub_one_of_pos (show 0 < x from interior_subset hx)
      linarith
  have hzero : cubicGap 1 = 0 := by norm_num [cubicGap]
  have hderiv : HasDerivAt cubicGap 0 1 := by
    simpa [cubicGapDeriv] using cubicGap_hasDerivAt (by norm_num : (0 : ℝ) < 1)
  have h1 : (1 : ℝ) ∈ Set.Ioi 0 := by norm_num
  have ht' : t ∈ Set.Ioi (0 : ℝ) := ht
  rcases lt_trichotomy t 1 with h | rfl | h
  · have hs := hconv.slope_le_of_hasDerivAt ht' h1 h hderiv
    rw [slope_def_field, hzero] at hs
    rcases div_nonpos_iff.mp hs with hwrong | hright
    · linarith [hwrong.2]
    · linarith [hright.1]
  · rw [hzero]
  · have hs := hconv.deriv_le_slope h1 ht' h hderiv.differentiableAt
    rw [hderiv.deriv, slope_def_field, hzero] at hs
    rcases div_nonneg_iff.mp hs with hright | hwrong
    · linarith [hright.1]
    · linarith [hwrong.2]

theorem cubic_log_majorant {d A : ℝ} (hd : 0 < d) (hA : 0 < A) :
    d ^ 2 * Real.log d ≤ (Real.log A + 3 / 2) * d ^ 2 - 2 * A * d + A ^ 2 / 2 +
      (d - A) ^ 3 / (3 * A) := by
  let t := d / A
  have ht : 0 < t := div_pos hd hA
  have hdt : d = A * t := by dsimp [t]; field_simp
  rw [hdt, Real.log_mul hA.ne' ht.ne']
  have he :
      (Real.log A + 3 / 2) * (A * t) ^ 2 - 2 * A * (A * t) + A ^ 2 / 2 +
        (A * t - A) ^ 3 / (3 * A) -
        (A * t) ^ 2 * (Real.log A + Real.log t) = A ^ 2 * cubicGap t := by
    unfold cubicGap
    field_simp [hA.ne']
    ring
  exact sub_nonneg.mp (he.symm ▸ mul_nonneg (sq_nonneg A) (cubicGap_nonneg ht))

theorem log_tangent_majorant {d A : ℝ} (hd : 0 < d) (hA : 0 < A) :
    Real.log d ≤ Real.log A + (d - A) / A := by
  have h := Real.log_le_sub_one_of_pos (div_pos hd hA)
  rw [Real.log_div hd.ne' hA.ne'] at h
  rw [sub_div, div_self hA.ne']
  linarith

private theorem quadratic_gap_hasDerivAt {t : ℝ} (ht : t < 1) :
    HasDerivAt (fun x : ℝ => -x - x ^ 2 / 2 - Real.log (1 - x)) (t ^ 2 / (1 - t)) t := by
  have hdiff := (hasDerivAt_const t (1 : ℝ)).sub (hasDerivAt_id t)
  have hnonzero : 1 - t ≠ 0 := (sub_pos.mpr ht).ne'
  have h := ((hasDerivAt_id t).neg.sub (((hasDerivAt_id t).pow 2).div_const (2 : ℝ))).sub
    ((Real.hasDerivAt_log hnonzero).comp t hdiff)
  convert h using 1
  simp only [id_eq]
  field_simp [hnonzero]
  ring

theorem log_quadratic_majorant {x : ℝ} (hx : 0 ≤ x) (h1 : x < 1) :
    Real.log (1 - x) ≤ -x - x ^ 2 / 2 := by
  let g := fun x : ℝ => -x - x ^ 2 / 2 - Real.log (1 - x)
  have hm : MonotoneOn g (Set.Iio 1) := by
    apply monotoneOn_of_deriv_nonneg (convex_Iio 1)
    · exact fun _ hy => (quadratic_gap_hasDerivAt hy).continuousAt.continuousWithinAt
    · intro y hy
      have hy' : y ∈ Set.Iio (1 : ℝ) := interior_subset hy
      exact (quadratic_gap_hasDerivAt (show y < 1 from hy')).differentiableAt.differentiableWithinAt
    · intro y hy
      have hy' : y ∈ Set.Iio (1 : ℝ) := interior_subset hy
      have hlt : y < 1 := hy'
      rw [(quadratic_gap_hasDerivAt hlt).deriv]
      exact div_nonneg (sq_nonneg _) (sub_nonneg.mpr hlt.le)
  have he := hm (show (0 : ℝ) ∈ Set.Iio 1 by norm_num) (show x ∈ Set.Iio 1 from h1) hx
  dsimp [g] at he
  norm_num at he
  linarith

end Mordell
