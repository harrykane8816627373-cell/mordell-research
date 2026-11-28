import Mordell.BernsteinTransform
import Mordell.BoxReplay
import Mordell.RationalLog
import Mathlib.Data.Fin.VecNotation

/-! Coefficientwise certification of the rational logarithmic upper bounds at box leaves. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def castPowerTerms {n : ℕ} (p : PowerTerms n ℚ) : PowerTerms n ℝ := p.map fun (e, a) => (e, (a : ℝ))

theorem axisPowerControl_ratCast (d : ℕ) (lo hi : ℚ) (i k : ℕ) :
    (axisPowerControl ℚ d lo hi i k : ℝ) = axisPowerControl ℝ d lo hi i k := by
  simp only [axisPowerControl, Rat.cast_sum, Rat.cast_mul, Rat.cast_pow, Rat.cast_natCast, Rat.cast_div, Rat.cast_sub]

theorem powerControls_ratCast {n : ℕ} (d : Fin n → ℕ) (b : RationalBox n) (p : PowerTerms n ℚ)
    (i : ∀ j, Fin (d j + 1)) :
    (powerControls (R := ℚ) d b.lower b.upper p i : ℝ) =
      powerControls d b.toBox.lower b.toBox.upper (castPowerTerms p) i := by
  induction p with
  | nil => simp [powerControls, castPowerTerms]
  | cons t p ih =>
    have h : ((∏ j, axisPowerControl ℚ (d j) (b.lower j) (b.upper j) (i j).val (t.1 j) : ℚ) : ℝ) =
        ∏ j, axisPowerControl ℝ (d j) (b.toBox.lower j) (b.toBox.upper j) (i j).val (t.1 j) := by
      simp only [Rat.cast_prod, axisPowerControl_ratCast, RationalBox.toBox]
    simpa only [powerControls, castPowerTerms, List.map_cons, List.sum_cons, Rat.cast_add,
      Rat.cast_mul, h] using congrArg ((t.2 : ℝ) * (∏ j, axisPowerControl ℝ (d j)
        (b.toBox.lower j) (b.toBox.upper j) (i j).val (t.1 j)) + ·) ih

theorem power_le_of_lower_controls {n : ℕ} (d : Fin n → ℕ) (b : RationalBox n)
    (p : PowerTerms n ℚ) (c : (∀ j, Fin (d j + 1)) → ℚ) (H : (Fin n → ℝ) → ℝ)
    (hp : ∀ t ∈ p, ∀ j, t.1 j ≤ d j)
    (hc : ∀ i, powerControls d b.lower b.upper p i ≤ c i)
    (hH : ∀ u, (∀ j, 0 ≤ u j ∧ u j ≤ 1) →
      tensorBernstein d (fun i => (c i : ℝ)) u ≤
        H (fun j => (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u j)) :
    ∀ x, b.Contains x → powerEval (castPowerTerms p) x ≤ H x := by
  intro x hx
  obtain ⟨u, hu, he⟩ := RationalBox.affineUnit_coordinates hx
  have hx' : x = fun j => (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u j := funext he
  rw [hx']
  have hpd : ∀ t ∈ castPowerTerms p, ∀ j, t.1 j ≤ d j := by
    intro t ht j
    obtain ⟨q, hq, rfl⟩ := List.mem_map.mp ht
    exact hp q hq j
  change powerEval (castPowerTerms p) (fun j => b.toBox.lower j + (b.toBox.upper j - b.toBox.lower j) * u j) ≤ _
  rw [← powerControls_sound d b.toBox.lower b.toBox.upper u _ hpd]
  refine (tensorBernstein_mono hu ?_).trans (hH u hu)
  intro i
  rw [← powerControls_ratCast]
  exact (Rat.cast_le (K := ℝ)).2 (hc i)

def affineLogTerms (a upper : ℚ) : PowerTerms 3 ℚ :=
  [(![0, 0, 0], upper - 1), (![1, 0, 0], 1 / a), (![0, 1, 0], 1 / a), (![0, 0, 1], -2 / a)]

abbrev tangentTerms (terms : ℕ) (a : ℚ) : PowerTerms 3 ℚ := affineLogTerms a (logUpperRat terms a)

def affineKnot (d : Fin 3 → ℕ) (b : RationalBox 3) (i : ∀ j, Fin (d j + 1)) (j : Fin 3) : ℚ :=
  b.lower j + (b.upper j - b.lower j) * ((i j).val : ℚ) / (d j : ℚ)

private theorem axis_zero (d i : ℕ) (lo hi : ℚ) : axisPowerControl ℚ d lo hi i 0 = 1 := by
  simp [axisPowerControl]

private theorem axis_one (d i : ℕ) (lo hi : ℚ) :
    axisPowerControl ℚ d lo hi i 1 = lo + (hi - lo) * (i : ℚ) / (d : ℚ) := by
  simp [axisPowerControl, Finset.sum_range_succ]
  ring

theorem affineLogControls_eq (d : Fin 3 → ℕ) (b : RationalBox 3) (a upper : ℚ)
    (i : ∀ j, Fin (d j + 1)) :
    powerControls d b.lower b.upper (affineLogTerms a upper) i =
      upper - 1 + gapDistance (affineKnot d b i) / a := by
  simp [powerControls, affineLogTerms, Fin.prod_univ_three, axis_zero, axis_one, gapDistance, affineKnot]
  ring

theorem tangentTerms_degree (terms : ℕ) (a : ℚ) (d : Fin 3 → ℕ) (hd : ∀ j, 1 ≤ d j) :
    ∀ t ∈ tangentTerms terms a, ∀ j, t.1 j ≤ d j := by
  intro t ht j
  simp only [tangentTerms, affineLogTerms, List.mem_cons, List.not_mem_nil, or_false] at ht
  have hj := hd j
  rcases ht with rfl | rfl | rfl | rfl <;> fin_cases j <;>
    first | exact hj | exact Nat.zero_le _

theorem tangentTerms_eval (terms : ℕ) (a : ℚ) (ha : 0 < a) (x : Fin 3 → ℝ) :
    powerEval (castPowerTerms (tangentTerms terms a)) x =
      (logUpperRat terms a : ℝ) + (gapDistance x - a) / a := by
  have h : (a : ℝ) ≠ 0 := by exact_mod_cast ha.ne'
  simp [powerEval, castPowerTerms, tangentTerms, affineLogTerms, gapDistance, Fin.prod_univ_three]
  field_simp [h]
  ring

theorem log_tangent_controls_sound (terms : ℕ) {a : ℚ} (ha : 0 < a) (d : Fin 3 → ℕ)
    (hd : ∀ j, 1 ≤ d j) (b : RationalBox 3) (c : (∀ j, Fin (d j + 1)) → ℚ) (H : (Fin 3 → ℝ) → ℝ)
    (hc : ∀ i, powerControls d b.lower b.upper (tangentTerms terms a) i ≤ c i)
    (hH : ∀ u, (∀ j, 0 ≤ u j ∧ u j ≤ 1) →
      tensorBernstein d (fun i => (c i : ℝ)) u ≤
        H (fun j => (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u j)) :
    ∀ x, b.Contains x → 0 < gapDistance x → Real.log (gapDistance x) ≤ H x := by
  intro x hx hdx
  have ha' : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have h := power_le_of_lower_controls d b (tangentTerms terms a) c H (tangentTerms_degree terms a d hd)
    hc hH x hx
  rw [tangentTerms_eval terms a ha] at h
  exact (log_tangent_majorant hdx ha').trans
    ((add_le_add (logUpperRat_sound ha terms) le_rfl).trans h)

theorem log_deficit_quadratic_rat {d : ℝ} {a : ℚ} (ha : 0 < a) (hd : 0 < d) (hle : d ≤ (a : ℝ))
    (terms : ℕ) : Real.log d ≤ (logUpperRat terms a : ℝ) - (1 - d / a) - (1 - d / a) ^ 2 / 2 := by
  have ha' : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have h₀ : 0 ≤ 1 - d / (a : ℝ) := sub_nonneg.mpr ((div_le_one ha').mpr hle)
  have h₁ : 1 - d / (a : ℝ) < 1 := by linarith [div_pos hd ha']
  have h := log_quadratic_majorant h₀ h₁
  rw [show 1 - (1 - d / (a : ℝ)) = d / a by ring, Real.log_div hd.ne' ha'.ne'] at h
  linarith [logUpperRat_sound ha terms]

theorem cleared_cubic_log_rat {d : ℝ} {a : ℚ} (ha : 0 < a) (hd : 0 < d) (terms : ℕ) :
    d ^ 2 * Real.log d ≤ ((logUpperRat terms a : ℝ) + 3 / 2) * d ^ 2 - 2 * (a : ℝ) * d +
      (a : ℝ) ^ 2 / 2 + (d - a) ^ 3 / (3 * a) := by
  have ha' : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have h := mul_le_mul_of_nonneg_right (logUpperRat_sound ha terms) (sq_nonneg d)
  linarith [cubic_log_majorant hd ha']

end Mordell
