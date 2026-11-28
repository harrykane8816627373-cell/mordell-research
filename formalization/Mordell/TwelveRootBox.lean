import Mordell.BoxRoot
import Mordell.ForceSecant
import Mordell.RationalBall
import Mordell.TwelveRootData
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Data.Rat.Floor

/-! Kernel replay of the collision and midpoint tests for the archived twelve-point root box. -/

open scoped BigOperators Matrix NNReal

noncomputable section

namespace Mordell
namespace TwelveRootBox

def centerRat (i : Fin 23) : ℚ := (TwelveRootData.centerNumerators i : ℚ) / 2 ^ TwelveRootData.centerBits

def center (i : Fin 23) : ℝ := centerRat i

def radiusRat : ℚ := 1 / 2 ^ TwelveRootData.radiusBits

def radius : ℝ := radiusRat

def preconditionerRat (i j : Fin 23) : ℚ :=
  (TwelveRootData.inverseNumerators i j : ℚ) / 2 ^ TwelveRootData.inverseBits

def preconditioner : Matrix (Fin 23) (Fin 23) ℝ := fun i j => preconditionerRat i j

def box : Set (GaugeCoordinates 11) := Metric.closedBall center radius

private theorem centerComponent (i : Fin 12) (a : Fin 2) :
    gaugeRead center i a = (gaugeRead (S := ℚ) centerRat i a : ℚ) :=
  gaugeRead_ratCast centerRat i a

theorem radius_pos : 0 < radius := by
  norm_num [radius, radiusRat, TwelveRootData.radiusBits]

private theorem midpoint_separation_i (i j : Fin 12) : i = j ∨
    2 * radiusRat < |gaugeRead centerRat i 0 - gaugeRead centerRat j 0| ∨
    2 * radiusRat < |gaugeRead centerRat i 1 - gaugeRead centerRat j 1| := by
  fin_cases i <;> fin_cases j <;> decide +kernel

private theorem coordinate_bounds {x : GaugeCoordinates 11} (hx : x ∈ box) (i : Fin 12) (a : Fin 2) :
    |gaugeRead x i a - gaugeRead center i a| ≤ radius := by
  have h (j : Fin 23) : |x j - center j| ≤ radius := by
    simpa only [Real.dist_eq] using (dist_pi_le_iff radius_pos.le).1 (Metric.mem_closedBall.mp hx) j
  refine Fin.cases ?_ (fun j => ?_) i
  · by_cases ha : a = 0
    · simpa [gaugeRead, ha] using h 0
    · simpa [gaugeRead, ha] using radius_pos.le
  · exact h _

theorem injective_on_box {x : GaugeCoordinates 11} (hx : x ∈ box) :
    Function.Injective (gaugeConfiguration x) := by
  intro i j hij
  have hd (a : Fin 2) : |gaugeRead center i a - gaugeRead center j a| ≤ 2 * radius := by
    have he : gaugeRead x i a = gaugeRead x j a := by
      simpa only [component_gaugeConfiguration] using congrArg (fun z => complexComponent z a) hij
    calc
      _ = |-(gaugeRead x i a - gaugeRead center i a) + (gaugeRead x j a - gaugeRead center j a)| := by
        rw [← he]; ring_nf
      _ ≤ _ := by
        have hh := abs_add_le (-(gaugeRead x i a - gaugeRead center i a))
          (gaugeRead x j a - gaugeRead center j a)
        rw [abs_neg] at hh
        linarith [coordinate_bounds hx i a, coordinate_bounds hx j a]
  have hrat (a : Fin 2) : |gaugeRead centerRat i a - gaugeRead centerRat j a| ≤ 2 * radiusRat := by
    have h := hd a
    rw [centerComponent, centerComponent, radius] at h
    exact_mod_cast h
  rcases midpoint_separation_i i j with he | h0 | h1
  · exact he
  · exact (not_lt_of_ge (hrat 0) h0).elim
  · exact (not_lt_of_ge (hrat 1) h1).elim

theorem pivot_pos_on_box {x : GaugeCoordinates 11} (hx : x ∈ box) : 0 < x 0 := by
  have ha : radiusRat < centerRat 0 := by decide +kernel
  have h : radius < center 0 := (Rat.cast_lt (K := ℝ)).2 ha
  have he := (abs_le.mp (coordinate_bounds hx 0 0)).1
  simp only [gaugeRead, Fin.cases_zero, if_true] at he
  linarith

theorem stationary_iff_on_box {x : GaugeCoordinates 11} (hx : x ∈ box) :
    gaugeForce (11 / 2) x = 0 ↔ Stationary (11 / 2 : ℂ) (gaugeConfiguration x) := by
  simpa using gaugeForce_eq_zero_iff (11 / 2) x (pivot_pos_on_box hx).ne' (injective_on_box hx)

private def qCoord (i : Fin 12) (a : Fin 2) : ℚ := gaugeRead centerRat i a

private def qDistSq (i j : Fin 12) : ℚ :=
  (qCoord i 0 - qCoord j 0) ^ 2 + (qCoord i 1 - qCoord j 1) ^ 2

private def qForceTerm (i j : Fin 12) (a : Fin 2) : ℚ :=
  if i = j then 0 else (qCoord i a - qCoord j a) / qDistSq i j

private def qPenalty (i : Fin 12) (a : Fin 2) : ℚ := (11 / 2) * qCoord i a

private theorem qForceTerm_cast (i j : Fin 12) (a : Fin 2) :
    (qForceTerm i j a : ℝ) =
      if i = j then 0 else
        (complexComponent (gaugeConfiguration center i) a - complexComponent (gaugeConfiguration center j) a) /
        distanceSq (gaugeConfiguration center) i j := by
  unfold qForceTerm
  split_ifs with h
  · norm_num
  · rw [component_gaugeConfiguration, component_gaugeConfiguration]
    simp only [qDistSq, qCoord, distanceSq, Complex.normSq_apply, Complex.sub_re, Complex.sub_im,
      gaugeConfiguration, centerComponent, pow_two, Rat.cast_div, Rat.cast_mul, Rat.cast_add, Rat.cast_sub]

private def roundingScale : ℚ := 2 ^ (224 : ℕ)

private def roundedLower (q : ℚ) : ℤ := ⌊q * roundingScale⌋

private def scaledResidualLower (i : Fin 12) (a : Fin 2) : ℤ :=
  (List.ofFn (fun j => roundedLower (qForceTerm i j a))).sum - roundedLower (qPenalty i a) - 1

private theorem midpoint_integer_test (i : Fin 12) (a : Fin 2) :
    -256 ≤ scaledResidualLower i a ∧ scaledResidualLower i a + 13 ≤ 256 := by
  fin_cases i <;> fin_cases a <;> decide +kernel

theorem midpoint_component_small (i : Fin 12) (a : Fin 2) :
    |complexComponent (cartesianResidual (11 / 2) (gaugeConfiguration center) i) a| ≤
      (256 : ℝ) / (roundingScale : ℝ) := by
  have hs : 0 < roundingScale := by norm_num [roundingScale]
  have hlo := Finset.sum_le_sum fun j (_ : j ∈ Finset.univ) => Int.floor_le (qForceTerm i j a * roundingScale)
  have hhi := Finset.sum_le_sum fun j (_ : j ∈ Finset.univ) =>
    (Int.lt_floor_add_one (qForceTerm i j a * roundingScale)).le
  have hfloor (q : ℚ) : ⌊q * roundingScale⌋ = roundedLower q := rfl
  simp only [← Finset.sum_mul, Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul, mul_one, hfloor, Nat.cast_ofNat] at hlo hhi
  let lo := scaledResidualLower i a
  have hqlo : (lo : ℚ) ≤ ((∑ j, qForceTerm i j a) - qPenalty i a) * roundingScale := by
    simp only [lo, scaledResidualLower, List.sum_ofFn, Int.cast_sub, Int.cast_sum, Int.cast_one, sub_mul]
    have hp : qPenalty i a * roundingScale ≤ (roundedLower (qPenalty i a) : ℚ) + 1 :=
      (Int.lt_floor_add_one _).le
    linarith only [hlo, hp]
  have hqhi : ((∑ j, qForceTerm i j a) - qPenalty i a) * roundingScale ≤ (lo : ℚ) + 13 := by
    simp only [lo, scaledResidualLower, List.sum_ofFn, Int.cast_sub, Int.cast_sum, Int.cast_one, sub_mul]
    have hp : (roundedLower (qPenalty i a) : ℚ) ≤ qPenalty i a * roundingScale := Int.floor_le _
    linarith only [hhi, hp]
  have hnum := midpoint_integer_test i a
  have hlow : -(256 : ℚ) ≤ (lo : ℚ) := by exact_mod_cast hnum.1
  have hhigh : (lo : ℚ) + 13 ≤ (256 : ℚ) := by exact_mod_cast hnum.2
  have ht : |(∑ j, qForceTerm i j a) - qPenalty i a| ≤ (256 : ℚ) / roundingScale := by
    rw [abs_le, ← neg_div]
    exact ⟨(div_le_iff₀ hs).2 (hlow.trans hqlo), (le_div_iff₀ hs).2 (hqhi.trans hhigh)⟩
  have hreal := (Rat.cast_le (K := ℝ)).2 ht
  have hpen : (qPenalty i a : ℝ) = (11 / 2 : ℝ) * complexComponent (gaugeConfiguration center i) a := by
    simp [qPenalty, qCoord, component_gaugeConfiguration, centerComponent]
  simpa only [Rat.cast_abs, Rat.cast_sub, Rat.cast_sum, Rat.cast_div, Rat.cast_ofNat, qForceTerm_cast,
    hpen, cartesianResidual_component] using hreal

private theorem preconditioner_rows_rat (i : Fin 23) :
    (List.ofFn (fun j => |preconditionerRat i j|)).sum ≤ 1024 := by
  fin_cases i <;> decide +kernel

theorem preconditioner_rows (i : Fin 23) : ∑ j, |preconditioner i j| ≤ 1024 := by
  have hr := (Rat.cast_le (K := ℝ)).2 (preconditioner_rows_rat i)
  simpa only [preconditioner, List.sum_ofFn, Rat.cast_sum, Rat.cast_abs, Rat.cast_ofNat] using hr

theorem midpoint_preconditioned_small :
    dist (preconditioner *ᵥ gaugeForce (n := 11) (11 / 2) center) 0 ≤ (1 : ℝ) / 2 ^ (206 : ℕ) := by
  have hn : (0 : ℝ) ≤ (256 : ℝ) / (roundingScale : ℝ) := by norm_num [roundingScale]
  have hm : dist (gaugeForce (n := 11) (11 / 2) center) 0 ≤ (256 : ℝ) / (roundingScale : ℝ) := by
    apply (dist_pi_le_iff hn).2
    intro i
    simp only [gaugeForce, Pi.zero_apply, Real.dist_eq, sub_zero]
    refine Fin.cases ?_ (fun j => ?_) i
    · simpa [gaugeProjection, complexComponent] using midpoint_component_small 0 0
    · obtain ⟨⟨k, a⟩, rfl⟩ := finProdFinEquiv.surjective j
      rw [gaugeProjection_apply_pair]
      exact midpoint_component_small k.succ a
  calc
    _ ≤ 1024 * dist (gaugeForce (n := 11) (11 / 2) center) 0 :=
      dist_mulVec_zero_le (by norm_num) preconditioner_rows _
    _ ≤ 1024 * ((256 : ℝ) / (roundingScale : ℝ)) := by gcongr
    _ = _ := by norm_num [roundingScale]

theorem existsUnique_stationary_of_slope
    (hslope : ∀ x ∈ box, ∀ y ∈ box,
      ∀ i, ∑ j, |(1 - preconditioner * gaugeSecant (11 / 2) x y) i j| ≤ 1 / 4) :
    ∃! x ∈ box, Stationary (11 / 2 : ℂ) (gaugeConfiguration x) := by
  have hres : dist (preconditioner *ᵥ gaugeForce (n := 11) (11 / 2) center) 0 ≤
      (1 - ((1 / 4 : ℝ≥0) : ℝ)) * radius := by
    apply midpoint_preconditioned_small.trans
    norm_num [radius, radiusRat, TwelveRootData.radiusBits]
  have hs : ∀ x ∈ box, ∀ y ∈ box,
      ∃ D : Matrix (Fin 23) (Fin 23) ℝ,
        gaugeForce (11 / 2) x - gaugeForce (11 / 2) y = D *ᵥ (x - y) ∧
        ∀ i, ∑ j, |(1 - preconditioner * D) i j| ≤ 1 / 4 := by
    intro x hx y hy
    exact ⟨gaugeSecant (11 / 2) x y,
      gaugeForce_sub_eq_secant (11 / 2) x y (injective_on_box hx) (injective_on_box hy), hslope x hx y hy⟩
  obtain ⟨x, ⟨hx, hf⟩, hu⟩ := existsUnique_root_of_box_slope (gaugeForce (n := 11) (11 / 2)) preconditioner center
    radius_pos.le (by norm_num : (1 / 4 : ℝ≥0) < 1) hres (by simpa [box] using hs)
  exact ⟨x, ⟨hx, (stationary_iff_on_box hx).1 hf⟩, fun y ⟨hy, h⟩ =>
    hu y ⟨hy, (stationary_iff_on_box hy).2 h⟩⟩

private def differenceBall (i j : Fin 12) (a : Fin 2) : RationalBall :=
  ⟨qCoord i a - qCoord j a, 2 * radiusRat⟩

private theorem differenceBall_contains (i j : Fin 12) (a : Fin 2) {x : GaugeCoordinates 11} (hx : x ∈ box) :
    (differenceBall i j a).Contains (pairDifferences x i j a) := by
  have he : pairDifferences x i j a - ((differenceBall i j a).mid : ℝ) =
      (gaugeRead x i a - gaugeRead center i a) - (gaugeRead x j a - gaugeRead center j a) := by
    simp only [pairDifferences, differenceBall, qCoord, Rat.cast_sub, centerComponent]
    ring
  rw [RationalBall.Contains, he, sub_eq_add_neg]
  have h := (abs_add_le _ _).trans (add_le_add (coordinate_bounds hx i a)
    (show |-(gaugeRead x j a - gaugeRead center j a)| ≤ radius by
      simpa only [abs_neg] using coordinate_bounds hx j a))
  change _ ≤ ((2 * radiusRat : ℚ) : ℝ)
  push_cast
  change _ ≤ 2 * radius
  linarith only [h]

private def distanceBall (i j : Fin 12) : RationalBall :=
  ((differenceBall i j 0).mul (differenceBall i j 0)).add
    ((differenceBall i j 1).mul (differenceBall i j 1))

private theorem distanceBall_contains (i j : Fin 12) {x : GaugeCoordinates 11} (hx : x ∈ box) :
    (distanceBall i j).Contains (distanceSq (gaugeConfiguration x) i j) := by
  have h0 := differenceBall_contains i j 0 hx
  have h1 := differenceBall_contains i j 1 hx
  simpa only [distanceBall, pairDifferences, distanceSq, gaugeConfiguration, Complex.normSq_apply,
    Complex.sub_re, Complex.sub_im] using RationalBall.add_sound (RationalBall.mul_sound h0 h0)
    (RationalBall.mul_sound h1 h1)

private theorem distanceBall_safe (i j : Fin 12) : i = j ∨ (distanceBall i j).rad < |(distanceBall i j).mid| := by
  fin_cases i <;> fin_cases j <;> decide +kernel

def pairSecantBall (i j : Fin 12) (a b : Fin 2) : RationalBall :=
  if i = j then RationalBall.point 0 else
    let inv := (distanceBall i j).inv
    let diag := if a = b then inv else RationalBall.point 0
    diag.sub (((differenceBall i j a).mul ((differenceBall i j b).add (differenceBall i j b))).mul (inv.mul inv))

theorem pairSecantBall_contains (i j : Fin 12) (a b : Fin 2) {x y : GaugeCoordinates 11}
    (hx : x ∈ box) (hy : y ∈ box) : (pairSecantBall i j a b).Contains (pairSecant x y i j a b) := by
  by_cases hij : i = j
  · simpa [pairSecantBall, pairSecant, hij] using RationalBall.point_sound 0
  have hs := (distanceBall_safe i j).resolve_left hij
  have hix := RationalBall.inv_sound (distanceBall_contains i j hx) hs
  have hiy := RationalBall.inv_sound (distanceBall_contains i j hy) hs
  have hdiag : (if a = b then (distanceBall i j).inv else RationalBall.point 0).Contains
      (if a = b then (distanceSq (gaugeConfiguration x) i j)⁻¹ else 0) := by
    split_ifs <;> first | exact hix | simpa using RationalBall.point_sound 0
  have hmul := RationalBall.mul_sound (RationalBall.mul_sound (differenceBall_contains i j a hy)
    (RationalBall.add_sound (differenceBall_contains i j b hx) (differenceBall_contains i j b hy)))
    (RationalBall.mul_sound hix hiy)
  simpa only [pairSecantBall, pairSecant, if_neg hij, div_eq_mul_inv, mul_inv] using
    RationalBall.sub_sound hdiag hmul

end TwelveRootBox
end Mordell
