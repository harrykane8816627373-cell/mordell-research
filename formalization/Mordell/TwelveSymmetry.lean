import Mordell.ExactTwelve

/-! The archived rotation and reflection maps preserve the isolated twelve-point root. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

private theorem stationary_relabel {n : ℕ} {z : Configuration n} {s : ℂ}
    (hs : Stationary s z) (e : Fin n ≃ Fin n) : Stationary s (z ∘ e) := by
  intro i
  have h : force (z ∘ e) i = force z (e i) := by
    unfold force
    convert Equiv.sum_comp e (fun j => if e i = j then (0 : ℂ) else (z (e i) - z j)⁻¹) using 1
    simp [e.injective.eq_iff]
  exact h.trans (hs (e i))

private theorem stationary_conj {n : ℕ} {z : Configuration n} {s : ℝ}
    (hs : Stationary (s : ℂ) z) : Stationary (s : ℂ) (fun i => conj (z i)) := by
  intro i
  have hf : force (fun j => conj (z j)) i = conj (force z i) := by
    simp only [force, map_sum]
    exact Finset.sum_congr rfl fun j _ => by split_ifs <;> simp
  rw [hf, hs i, map_mul, Complex.conj_ofReal]

private theorem configuration_projection {n : ℕ} {z : Configuration (n + 1)} (h : (z 0).im = 0) :
    gaugeConfiguration (gaugeProjection z) = z := by
  funext i
  refine Fin.cases ?_ (fun j => ?_) i
  · apply Complex.ext
    · rfl
    · exact h.symm
  · apply Complex.ext <;> simp only [gaugeConfiguration, gaugeRead, Fin.cases_succ,
      gaugeProjection_apply_pair, complexComponent, if_true, Fin.one_eq_zero_iff,
      Nat.reduceEqDiff, if_false]

namespace TwelveRootBox

def rotationPerm : Fin 12 ≃ Fin 12 :=
  ⟨![2, 8, 4, 10, 0, 3, 9, 6, 11, 7, 5, 1], ![4, 11, 0, 5, 2, 10, 7, 9, 1, 6, 3, 8],
    by intro i; fin_cases i <;> decide, by intro i; fin_cases i <;> decide⟩

def reflectionPerm : Fin 12 ≃ Fin 12 :=
  ⟨![11, 4, 8, 3, 1, 10, 9, 7, 2, 6, 5, 0], ![11, 4, 8, 3, 1, 10, 9, 7, 2, 6, 5, 0],
    by intro i; fin_cases i <;> decide, by intro i; fin_cases i <;> decide⟩

private def actionPerm (reflection : Bool) : Fin 12 ≃ Fin 12 :=
  if reflection then reflectionPerm else rotationPerm

def symmetryImage (reflection : Bool) : Configuration 12 :=
  let w := exactTwelveConfiguration ∘ actionPerm reflection
  if reflection then (fun i => conj (w i)) else w

def regaugedSymmetry (reflection : Bool) : Configuration 12 :=
  let w := symmetryImage reflection
  scale (conj (w 0) / (‖w 0‖ : ℝ)) w

def rootComponentBall (i : Fin 12) (a : Fin 2) : RationalBall :=
  ⟨gaugeRead centerRat i a, 1 / 2 ^ (205 : ℕ)⟩

theorem rootComponentBall_contains (i : Fin 12) (a : Fin 2) :
    (rootComponentBall i a).Contains (complexComponent (exactTwelveConfiguration i) a) := by
  have hd := stationary_close_to_center isolatedCoordinates_spec.1 isolatedCoordinates_spec.2
  have hcomp (p : Fin 23) : |isolatedCoordinates p - center p| ≤ (1 : ℝ) / 2 ^ (205 : ℕ) := by
    have h := (dist_pi_le_iff (by positivity : (0 : ℝ) ≤ 1 / 2 ^ (205 : ℕ))).1 hd p
    simpa only [Real.dist_eq] using h
  have hc : gaugeRead center i a = ((gaugeRead (S := ℚ) centerRat i a : ℚ) : ℝ) := gaugeRead_ratCast centerRat i a
  simp only [RationalBall.Contains, rootComponentBall, exactTwelveConfiguration, component_gaugeConfiguration,
    ← hc, Rat.cast_div, Rat.cast_pow, Rat.cast_ofNat, Rat.cast_one]
  refine Fin.cases ?_ (fun j => ?_) i
  · by_cases ha : a = 0
    · simpa only [gaugeRead, Fin.cases_zero, ha, if_true] using hcomp 0
    · simp [gaugeRead, ha]
  · exact hcomp _

private def sourceBall (f : Bool) (i : Fin 12) (a : Fin 2) : RationalBall :=
  if f && decide (a = 1) then (rootComponentBall (actionPerm f i) a).neg else rootComponentBall (actionPerm f i) a

private theorem sourceBall_contains (f : Bool) (i : Fin 12) (a : Fin 2) :
    (sourceBall f i a).Contains (complexComponent (symmetryImage f i) a) := by
  cases f <;> fin_cases a
  · exact rootComponentBall_contains _ _
  · exact rootComponentBall_contains _ _
  · simpa [sourceBall, symmetryImage, Function.comp_def, complexComponent] using
      rootComponentBall_contains (actionPerm true i) 0
  · simpa [sourceBall, symmetryImage, Function.comp_def, complexComponent] using
      RationalBall.neg_sound (rootComponentBall_contains (actionPerm true i) 1)

private def squareRadius (f : Bool) : RationalBall :=
  ((sourceBall f 0 0).mul (sourceBall f 0 0)).add ((sourceBall f 0 1).mul (sourceBall f 0 1))

private def normBall : RationalBall := ⟨centerRat 0, 1 / 2 ^ (140 : ℕ)⟩

private theorem normBall_check (f : Bool) :
    |(squareRadius f).mid - (normBall.mid)^2| + (squareRadius f).rad ≤ normBall.rad := by
  cases f <;> decide +kernel

private theorem normBall_contains (f : Bool) : normBall.Contains ‖symmetryImage f 0‖ := by
  let r := ‖symmetryImage f 0‖
  let m : ℝ := normBall.mid
  let e : ℝ := normBall.rad
  have h₁ : 1 ≤ m := by
    simpa only [m, Rat.cast_one] using (Rat.cast_le (K := ℝ)).2 (by decide +kernel : 1 ≤ normBall.mid)
  have hr : 0 ≤ r := norm_nonneg _
  have hs : (squareRadius f).Contains (r ^ 2) := by
    have h0 := sourceBall_contains f 0 0
    have h1 := sourceBall_contains f 0 1
    change (squareRadius f).Contains (‖symmetryImage f 0‖ ^ 2)
    rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply]
    simpa only [squareRadius, complexComponent, if_true, Fin.one_eq_zero_iff, Nat.reduceEqDiff, if_false] using
      RationalBall.add_sound (RationalBall.mul_sound h0 h0) (RationalBall.mul_sound h1 h1)
  have hc : |((squareRadius f).mid : ℝ) - m ^ 2| + (squareRadius f).rad ≤ e := by
    have h := (Rat.cast_le (K := ℝ)).2 (normBall_check f)
    simpa only [m, e, Rat.cast_abs, Rat.cast_sub, Rat.cast_add, Rat.cast_pow] using h
  have hsq : |r ^ 2 - m ^ 2| ≤ e := by
    have ht : |r ^ 2 - ((squareRadius f).mid : ℝ)| ≤ (squareRadius f).rad := hs
    rw [show r ^ 2 - m ^ 2 = (r ^ 2 - ((squareRadius f).mid : ℝ)) + ((squareRadius f).mid - m ^ 2) by ring]
    exact (abs_add_le _ _).trans (by linarith)
  have hm : 1 ≤ r + m := by linarith
  rw [show r ^ 2 - m ^ 2 = (r - m) * (r + m) by ring, abs_mul,
    abs_of_nonneg (show 0 ≤ r + m by linarith)] at hsq
  change |r - m| ≤ e
  simpa only [mul_one] using (mul_le_mul_of_nonneg_left hm (abs_nonneg (r - m))).trans hsq

private theorem symmetry_pivot_nonzero (f : Bool) : symmetryImage f 0 ≠ 0 := by
  have h : |‖symmetryImage f 0‖ - (normBall.mid : ℝ)| ≤ (normBall.rad : ℝ) := normBall_contains f
  intro hz
  rw [hz, norm_zero, zero_sub, abs_neg] at h
  have hm : (normBall.rad : ℝ) < |(normBall.mid : ℝ)| := by
    simpa only [Rat.cast_abs] using (Rat.cast_lt (K := ℝ)).2
      (by decide +kernel : normBall.rad < |normBall.mid|)
  exact (not_lt_of_ge h) hm

private def imageBall (f : Bool) (i : Fin 12) (a : Fin 2) : RationalBall :=
  let x := sourceBall f 0 0
  let y := sourceBall f 0 1
  let u := sourceBall f i 0
  let v := sourceBall f i 1
  (if a = 0 then (x.mul u).add (y.mul v) else (x.mul v).sub (y.mul u)).mul normBall.inv

private theorem imageBall_contains (f : Bool) (i : Fin 12) (a : Fin 2) :
    (imageBall f i a).Contains (complexComponent (regaugedSymmetry f i) a) := by
  have hn := RationalBall.inv_sound (normBall_contains f) (by decide +kernel : normBall.rad < |normBall.mid|)
  have x := sourceBall_contains f 0 0
  have y := sourceBall_contains f 0 1
  have u := sourceBall_contains f i 0
  have v := sourceBall_contains f i 1
  have he : regaugedSymmetry f i = (conj (symmetryImage f 0) * symmetryImage f i) / (‖symmetryImage f 0‖ : ℝ) := by
    simp only [regaugedSymmetry, scale, div_eq_mul_inv]
    ring
  rw [he]
  fin_cases a
  · have h := RationalBall.mul_sound (RationalBall.add_sound (RationalBall.mul_sound x u)
      (RationalBall.mul_sound y v)) hn
    rw [← div_eq_mul_inv] at h
    simpa [imageBall, complexComponent, Complex.mul_re, sub_eq_add_neg] using h
  · have h := RationalBall.mul_sound (RationalBall.sub_sound (RationalBall.mul_sound x v)
      (RationalBall.mul_sound y u)) hn
    rw [← div_eq_mul_inv] at h
    simpa [imageBall, complexComponent, Complex.mul_im, sub_eq_add_neg] using h

private theorem image_check (f : Bool) (i : Fin 12) : ∀ a : Fin 2,
    |(imageBall f i a).mid - gaugeRead centerRat i a| + (imageBall f i a).rad ≤ radiusRat := by
  cases f <;> fin_cases i <;> decide +kernel

private theorem image_close (f : Bool) (i : Fin 12) (a : Fin 2) :
    |complexComponent (regaugedSymmetry f i) a - gaugeRead center i a| ≤ radius := by
  let b := imageBall f i a
  have h : |complexComponent (regaugedSymmetry f i) a - (b.mid : ℝ)| ≤ (b.rad : ℝ) := imageBall_contains f i a
  have hb : |(b.mid : ℝ) - gaugeRead center i a| + (b.rad : ℝ) ≤ radius := by
    have hc : gaugeRead center i a = ((gaugeRead (S := ℚ) centerRat i a : ℚ) : ℝ) := gaugeRead_ratCast centerRat i a
    rw [hc]
    simpa only [b, Rat.cast_abs, Rat.cast_add, Rat.cast_sub, radius] using
      (Rat.cast_le (K := ℝ)).2 (image_check f i a)
  rw [show complexComponent (regaugedSymmetry f i) a - gaugeRead center i a =
    (complexComponent (regaugedSymmetry f i) a - (b.mid : ℝ)) + ((b.mid : ℝ) - gaugeRead center i a) by ring]
  exact (abs_add_le _ _).trans (by linarith)

/-- In the original coordinate order, exact regauging after either archived symmetry fixes the root. -/
theorem regauged_symmetry_fixed (f : Bool) : regaugedSymmetry f = exactTwelveConfiguration := by
  let w := symmetryImage f
  let s : ℂ := conj (w 0) / (‖w 0‖ : ℝ)
  have hnorm : 0 < ‖w 0‖ := norm_pos_iff.mpr (symmetry_pivot_nonzero f)
  have hunit : Complex.normSq s = 1 := by
    rw [show s = conj (w 0) / (‖w 0‖ : ℝ) from rfl, Complex.normSq_div, Complex.normSq_conj,
      Complex.normSq_ofReal, Complex.normSq_eq_norm_sq, ← sq, div_self (pow_ne_zero 2 hnorm.ne')]
  have hw : Stationary (11 / 2 : ℂ) w := by
    have h := stationary_relabel isolatedCoordinates_spec.2 (actionPerm f)
    cases f
    · exact h
    · have hh : Stationary ((11 / 2 : ℝ) : ℂ) (exactTwelveConfiguration ∘ actionPerm true) := by
        simpa using h
      simpa [w, symmetryImage] using stationary_conj hh
  have hs : Stationary (11 / 2 : ℂ) (regaugedSymmetry f) := by
    simpa only [regaugedSymmetry, hunit, Complex.ofReal_one, div_one] using hw.scale s
  have h0 : (regaugedSymmetry f 0).im = 0 := by
    change (((conj (w 0) / (‖w 0‖ : ℝ)) * w 0)).im = 0
    rw [div_mul_eq_mul_div, mul_comm (conj (w 0)), Complex.mul_conj, ← Complex.ofReal_div, Complex.ofReal_im]
  let x := gaugeProjection (regaugedSymmetry f)
  have hx : x ∈ box := by
    apply Metric.mem_closedBall.mpr
    apply (dist_pi_le_iff radius_pos.le).2
    intro p
    rw [Real.dist_eq]
    refine Fin.cases ?_ (fun q => ?_) p
    · simpa only [x, gaugeProjection, Fin.cases_zero, complexComponent, if_true, gaugeRead] using image_close f 0 0
    · obtain ⟨⟨j, a⟩, rfl⟩ := finProdFinEquiv.surjective q
      simpa only [x, gaugeProjection_apply_pair, gaugeRead, Fin.cases_succ] using image_close f j.succ a
  have hconfig : gaugeConfiguration x = regaugedSymmetry f := configuration_projection h0
  have hstat : Stationary (11 / 2 : ℂ) (gaugeConfiguration x) := hconfig.symm ▸ hs
  obtain ⟨m, _, hm⟩ := existsUnique_stationary
  have he : x = isolatedCoordinates := (hm x ⟨hx, hstat⟩).trans (hm _ isolatedCoordinates_spec).symm
  exact hconfig.symm.trans (congrArg gaugeConfiguration he)

end TwelveRootBox
end Mordell
