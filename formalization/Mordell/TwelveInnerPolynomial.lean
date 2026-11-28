import Mordell.InnerInitChecks
import Mordell.TwelveInnerCover

/-! The archived rational power lower envelope implies the signed inner pair comparison. -/

open scoped BigOperators

namespace Mordell.InnerCoverData

noncomputable section

private def axisControl (j : Fin 3) (i k : Fin (degree j + 1)) : ℝ :=
  axisPowerControl ℝ (degree j) (root.lower j) (root.upper j) i.val k.val

private theorem transform_sound (axis : Fin 3) (before after scale : ℕ)
    (m : Fin (degree axis + 1) → Fin (degree axis + 1) → ℤ) (src dst : GridIndex degree → ℤ)
    (hden : after = before * scale) (hnum : ∀ i, dst i = tensorAxisMap axis m src i)
    (hm : ∀ i j, (m i j : ℚ) / (scale : ℚ) =
      axisPowerControl ℚ (degree axis) (root.lower axis) (root.upper axis) i.val j.val) :
    realTensor after dst = tensorAxisMap axis (axisControl axis) (realTensor before src) := by
  have hm' : (fun i j => (m i j : ℝ) / (scale : ℝ)) = axisControl axis := by
    funext i j
    have h := congrArg (fun q : ℚ => (q : ℝ)) (hm i j)
    simpa only [Rat.cast_div, Rat.cast_intCast, Rat.cast_natCast, axisPowerControl_ratCast, axisControl] using h
  funext i
  calc
    _ = realTensor (before * scale) (tensorAxisMap axis m src) i := by
      rw [hden]
      simp only [realTensor, hnum]
    _ = _ := by rw [realTensor_axisMap, hm']

private theorem final_stage_sound : realTensor stageDen3 stage3 = rootPowerControls degree root (realTensor stageDen0 stage0) := by
  have h1 := transform_sound 0 stageDen0 stageDen1 axisDen0 axisMatrix0 stage0 stage1
    (by decide +kernel) stage1_checked matrix0_checked
  have h2 := transform_sound 1 stageDen1 stageDen2 axisDen1 axisMatrix1 stage1 stage2
    (by decide +kernel) stage2_checked matrix1_checked
  have h3 := transform_sound 2 stageDen2 stageDen3 axisDen2 axisMatrix2 stage2 stage3
    (by decide +kernel) stage3_checked matrix2_checked
  rw [h3, h2, h1]
  funext i
  exact tensorAxisMap_three axisControl (realTensor stageDen0 stage0) i

theorem archived_root_envelope : GridEnvelope degree root grid0 (densePowerEval (realTensor stageDen0 stage0)) := by
  apply gridEnvelope_of_densePower
  intro i
  rw [← final_stage_sound]
  simp only [gridControls, realTensor, Rat.cast_div, Rat.cast_intCast, Rat.cast_pow, Rat.cast_ofNat]
  apply (div_le_div_iff₀ (by positivity) (by norm_num [stageDen3])).2
  exact_mod_cast rounding_checked i

private theorem zeroIndex_iff (i : GridIndex degree) : i 0 = 0 ∧ i 1 = 0 ∧ i 2 = 0 ↔ i = zeroIndex := by
  constructor
  · rintro ⟨h0, h1, h2⟩
    funext j
    fin_cases j
    · exact h0
    · exact h1
    · exact h2
  · rintro rfl
    exact ⟨rfl, rfl, rfl⟩

private theorem trimmed_coefficient (i : GridIndex degree) :
    realTensor stageDen0 stage0 i = realTensor sourceDen sourceInside i +
      if i = zeroIndex then (discardedLower sourceCap sourceOutside : ℝ) else 0 := by
  classical
  have he : (stageDen0 : ℝ) = (inputMultiplier : ℝ) * (sourceDen : ℝ) := by
    exact_mod_cast (by decide +kernel : (stageDen0 : ℤ) = inputMultiplier * (sourceDen : ℤ))
  have hm : (inputMultiplier : ℝ) ≠ 0 := by norm_num [inputMultiplier]
  have hd : (sourceDen : ℝ) ≠ 0 := by norm_num [sourceDen]
  have hdiv : ((inputMultiplier : ℝ) * (sourceInside i : ℝ)) / (stageDen0 : ℝ) = realTensor sourceDen sourceInside i := by
    rw [he]
    unfold realTensor
    field_simp [hm, hd]
  have ht : ((tailNumerator : ℝ) / (stageDen0 : ℝ)) = (discardedLower sourceCap sourceOutside : ℝ) := by
    have h := congrArg (fun q : ℚ => (q : ℝ)) tail_checked
    simpa only [Rat.cast_div, Rat.cast_intCast, Rat.cast_natCast] using h
  unfold realTensor
  simp only [input_checked, zeroIndex_iff, Int.cast_add, Int.cast_mul, add_div]
  change ((inputMultiplier : ℝ) * (sourceInside i : ℝ)) / (stageDen0 : ℝ) + _ = _
  rw [hdiv]
  split_ifs <;> simp only [realTensor, Int.cast_zero, zero_div, ht]

private theorem trimmed_power_eval (x : Fin 3 → ℝ) :
    densePowerEval (realTensor stageDen0 stage0) x =
      densePowerEval (realTensor sourceDen sourceInside) x + (discardedLower sourceCap sourceOutside : ℝ) := by
  classical
  simp only [densePowerEval, trimmed_coefficient, add_mul, Finset.sum_add_distrib]
  congr 1
  rw [Finset.sum_eq_single zeroIndex]
  · simp [zeroIndex]
  · intro i _ hi
    simp [hi]
  · simp

private theorem in_source_cap {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    ∀ j, 0 ≤ x j ∧ x j ≤ (sourceCap j : ℝ) := by
  intro j
  fin_cases j
  · exact ⟨hd.1, (hx 0).2⟩
  · exact ⟨hd.2.1, (hx 1).2⟩
  · change 0 ≤ x 2 ∧ x 2 ≤ ((3 / 40 : ℚ) : ℝ)
    simpa using signed_inner_angular_cap hx hd

/-- The exported power polynomial, before discarding any monomials. -/
def archivedPowerPolynomial (x : Fin 3 → ℝ) : ℝ :=
  densePowerEval (realTensor sourceDen sourceInside) x + powerEval (castPowerTerms sourceOutside) x

/-- The signed inner pair comparison reaches the full archived rational lower power polynomial. -/
theorem signed_inner_archive_power_bound {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    Real.log (gapDistance x) ≤ archivedPowerPolynomial x := by
  have h := signed_inner_log_bound archived_root_envelope hx hd
  rw [trimmed_power_eval] at h
  exact h.trans (add_le_add le_rfl (discardedLower_sound sourceCap (in_source_cap hx hd) sourceOutside))

/-- It suffices to enclose the real kernel's coefficients in the archived basis from below. -/
theorem signed_inner_coefficients_bound (c : GridIndex degree → ℝ) (a : (Fin 3 → ℕ) × ℚ → ℝ)
    (hc : ∀ i, realTensor sourceDen sourceInside i ≤ c i) (ha : ∀ t ∈ sourceOutside, (t.2 : ℝ) ≤ a t)
    {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    Real.log (gapDistance x) ≤ densePowerEval c x + powerEval (sourceOutside.map fun t => (t.1, a t)) x := by
  have hn (j : Fin 3) : 0 ≤ x j := (in_source_cap hx hd j).1
  have h : densePowerEval (realTensor sourceDen sourceInside) x ≤ densePowerEval c x := by
    unfold densePowerEval
    exact Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_right (hc i)
      (Finset.prod_nonneg fun j _ => pow_nonneg (hn j) _)
  exact (signed_inner_archive_power_bound hx hd).trans
    (add_le_add h (powerEval_mono_same_exponents sourceOutside a hn ha))

end
end Mordell.InnerCoverData
