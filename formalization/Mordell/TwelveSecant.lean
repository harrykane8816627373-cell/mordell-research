import Mordell.TwelveRootBox

/-! The rational, finite row check for the archived twelve-point preconditioner. -/

open scoped BigOperators Matrix

noncomputable section

namespace Mordell.TwelveRootBox

private def refPair (i j : Fin 12) (a b : Fin 2) : ℚ :=
  (TwelveRootData.pairSecantNumerators i j a b : ℚ) / 2 ^ (40 : ℕ)

set_option maxRecDepth 16000 in
set_option maxHeartbeats 0 in
private theorem pair_check (i : Fin 12) : ∀ (j : Fin 12) (a b : Fin 2),
    |(pairSecantBall i j a b).mid - refPair i j a b| + (pairSecantBall i j a b).rad ≤ 1 / 2 ^ (32 : ℕ) := by
  fin_cases i <;> decide +kernel

private theorem pair_approx (i j : Fin 12) (a b : Fin 2) {x y : GaugeCoordinates 11}
    (hx : x ∈ box) (hy : y ∈ box) :
    |pairSecant x y i j a b - (refPair i j a b : ℝ)| ≤ 1 / 2 ^ (32 : ℕ) := by
  let B := pairSecantBall i j a b
  have hs : |pairSecant x y i j a b - (B.mid : ℝ)| ≤ B.rad := pairSecantBall_contains i j a b hx hy
  have hq : |(B.mid : ℝ) - (refPair i j a b : ℝ)| + B.rad ≤ 1 / 2 ^ (32 : ℕ) := by
    have h := (Rat.cast_le (K := ℝ)).2 (pair_check i j a b)
    simpa only [B, Rat.cast_abs, Rat.cast_add, Rat.cast_sub, Rat.cast_div, Rat.cast_pow, Rat.cast_ofNat,
      Rat.cast_one] using h
  rw [show pairSecant x y i j a b - (refPair i j a b : ℝ) =
    (pairSecant x y i j a b - (B.mid : ℝ)) + ((B.mid : ℝ) - (refPair i j a b : ℝ)) by ring]
  exact (abs_add_le _ _).trans (by linarith)

private def sumInt {n : ℕ} (f : Fin n → ℤ) : ℤ := (List.ofFn f).sum

private theorem sumInt_eq {n : ℕ} (f : Fin n → ℤ) : sumInt f = ∑ i, f i := List.sum_ofFn

private theorem sum_div_scalar {n : ℕ} (f : Fin n → ℝ) (d : ℝ) :
    (∑ i, f i) / d = ∑ i, f i / d := by simp only [div_eq_mul_inv, Finset.sum_mul]

private def refSecantComputed (p q : Fin 23) : ℤ :=
  let i := (gaugeAxis (n := 11) p).1
  let a := (gaugeAxis (n := 11) p).2
  (sumInt fun (j : Fin 12) => sumInt fun (b : Fin 2) => TwelveRootData.pairSecantNumerators i j a b *
    (gaugeIncidence i b q - gaugeIncidence j b q)) - 11 * 2 ^ (39 : ℕ) * gaugeIncidence i a q

set_option maxRecDepth 16000 in
set_option maxHeartbeats 0 in
private theorem refSecant_matches (p : Fin 23) : ∀ q, TwelveRootData.secantNumerators p q = refSecantComputed p q := by
  fin_cases p <;> decide +kernel

private def refSecant : Matrix (Fin 23) (Fin 23) ℝ :=
  fun p q => (TwelveRootData.secantNumerators p q : ℝ) / 2 ^ (40 : ℕ)

private theorem refSecant_apply (p q : Fin 23) :
    refSecant p q =
      let i := (gaugeAxis (n := 11) p).1
      let a := (gaugeAxis (n := 11) p).2
      (∑ j, ∑ b, (refPair i j a b : ℝ) * (gaugeIncidence i b q - gaugeIncidence j b q)) -
        (11 / 2 : ℝ) * gaugeIncidence i a q := by
  simp only [refSecant, refSecant_matches, refSecantComputed, sumInt_eq, Int.cast_sub, Int.cast_mul,
    Int.cast_pow, Int.cast_ofNat, Int.cast_sum, gaugeIncidence_cast, sub_div, sum_div_scalar,
    refPair, Rat.cast_div, Rat.cast_intCast, Rat.cast_pow, Rat.cast_ofNat]
  congr 1
  · exact Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => by ring
  · norm_num
    ring

private theorem secant_approx {x y : GaugeCoordinates 11} (hx : x ∈ box) (hy : y ∈ box) (p q : Fin 23) :
    |gaugeSecant (11 / 2) x y p q - refSecant p q| ≤ (24 : ℝ) / 2 ^ (32 : ℕ) := by
  let i := (gaugeAxis (n := 11) p).1
  let a := (gaugeAxis (n := 11) p).2
  let h (j : Fin 12) (b : Fin 2) : ℝ := gaugeIncidence i b q - gaugeIncidence j b q
  have he : gaugeSecant (11 / 2) x y p q - refSecant p q =
      ∑ j, ∑ b, (pairSecant x y i j a b - (refPair i j a b : ℝ)) * h j b := by
    rw [refSecant_apply]
    simp only [gaugeSecant, i, a, h, sub_mul, Finset.sum_sub_distrib]
    ring
  rw [he]
  calc
    _ ≤ ∑ j, ∑ b, |(pairSecant x y i j a b - (refPair i j a b : ℝ)) * h j b| :=
      (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun _ _ => Finset.abs_sum_le_sum_abs _ _))
    _ ≤ ∑ _j : Fin 12, ∑ _b : Fin 2, (1 : ℝ) / 2 ^ (32 : ℕ) := by
      refine Finset.sum_le_sum fun j _ => Finset.sum_le_sum fun b _ => ?_
      rw [abs_mul, ← mul_one ((1 : ℝ) / 2 ^ (32 : ℕ))]
      exact mul_le_mul (pair_approx i j a b hx hy) (abs_gaugeIncidence_sub_le_one i j b q)
        (abs_nonneg _) (by positivity)
    _ = _ := by norm_num

private def refErrorNum (i j : Fin 23) : ℤ :=
  (if i = j then 2 ^ (TwelveRootData.inverseBits + 40) else 0) -
    sumInt fun k => TwelveRootData.inverseNumerators i k * TwelveRootData.secantNumerators k j

set_option maxRecDepth 16000 in
set_option maxHeartbeats 0 in
private theorem ref_rows_check (i : Fin 23) :
    8 * sumInt (fun j => |refErrorNum i j|) ≤ 2 ^ (TwelveRootData.inverseBits + 40) := by
  fin_cases i <;> decide +kernel

private theorem refError_cast (i j : Fin 23) :
    (1 - preconditioner * refSecant) i j =
      (refErrorNum i j : ℝ) / 2 ^ (TwelveRootData.inverseBits + 40) := by
  simp only [refErrorNum, Matrix.sub_apply, Matrix.mul_apply, Matrix.one_apply, preconditioner,
    preconditionerRat, refSecant, Int.cast_sub, sumInt_eq, Int.cast_sum, Int.cast_mul, Rat.cast_div,
    Rat.cast_intCast, Rat.cast_pow, Rat.cast_ofNat, sub_div, sum_div_scalar]
  congr 1
  · split_ifs <;> simp
  · exact Finset.sum_congr rfl fun _ _ => by rw [pow_add]; ring

private theorem ref_rows (i : Fin 23) : ∑ j, |(1 - preconditioner * refSecant) i j| ≤ 1 / 8 := by
  have h : (8 : ℝ) * (∑ j, |(refErrorNum i j : ℝ)|) ≤ 2 ^ (TwelveRootData.inverseBits + 40) := by
    have hh := (Int.cast_le (R := ℝ)).2 (ref_rows_check i)
    simpa only [Int.cast_mul, Int.cast_pow, Int.cast_abs, Int.cast_ofNat, sumInt_eq, Int.cast_sum] using hh
  simp only [refError_cast, abs_div, abs_of_nonneg (by positivity : 0 ≤ (2 : ℝ) ^ (TwelveRootData.inverseBits + 40)),
    ← sum_div_scalar]
  apply (div_le_iff₀ (by positivity)).2
  linarith

theorem uniform_secant_rows {x y : GaugeCoordinates 11} (hx : x ∈ box) (hy : y ∈ box) (i : Fin 23) :
    ∑ j, |(1 - preconditioner * gaugeSecant (11 / 2) x y) i j| ≤ 1 / 4 := by
  let D := gaugeSecant (11 / 2) x y
  have herr (j : Fin 23) : |(preconditioner * refSecant) i j - (preconditioner * D) i j| ≤
      (1024 : ℝ) * (24 / 2 ^ (32 : ℕ)) := by
    simp only [Matrix.mul_apply, ← Finset.sum_sub_distrib, ← mul_sub]
    calc
      _ ≤ ∑ k, |preconditioner i k * (refSecant k j - D k j)| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ k, |preconditioner i k| * (24 / 2 ^ (32 : ℕ)) := by
        refine Finset.sum_le_sum fun k _ => ?_
        rw [abs_mul, abs_sub_comm (refSecant k j)]
        exact mul_le_mul_of_nonneg_left (secant_approx hx hy k j) (abs_nonneg _)
      _ ≤ _ := by
        rw [← Finset.sum_mul]
        exact mul_le_mul_of_nonneg_right (preconditioner_rows i) (by positivity)
  calc
    _ ≤ ∑ j, (|(1 - preconditioner * refSecant) i j| + 1024 * (24 / 2 ^ (32 : ℕ))) := by
      refine Finset.sum_le_sum fun j _ => ?_
      rw [show (1 - preconditioner * D) i j = (1 - preconditioner * refSecant) i j +
        ((preconditioner * refSecant) i j - (preconditioner * D) i j) by simp only [Matrix.sub_apply]; ring]
      exact (abs_add_le _ _).trans (add_le_add le_rfl (herr j))
    _ ≤ _ := by
      rw [Finset.sum_add_distrib]
      have h := ref_rows i
      norm_num [Finset.sum_const, Fintype.card_fin, Matrix.sub_apply] at h ⊢
      linarith

/-- The original 23-coordinate dyadic box contains exactly one stationary solution. -/
theorem existsUnique_stationary : ∃! x ∈ box, Stationary (11 / 2 : ℂ) (gaugeConfiguration x) :=
  existsUnique_stationary_of_slope fun _ hx _ hy => uniform_secant_rows hx hy

def isolatedCoordinates : GaugeCoordinates 11 := Classical.choose existsUnique_stationary.exists

theorem isolatedCoordinates_spec : isolatedCoordinates ∈ box ∧
    Stationary (11 / 2 : ℂ) (gaugeConfiguration isolatedCoordinates) :=
  Classical.choose_spec existsUnique_stationary.exists

theorem stationary_close_to_center {x : GaugeCoordinates 11} (hx : x ∈ box)
    (hs : Stationary (11 / 2 : ℂ) (gaugeConfiguration x)) : dist x center ≤ (1 : ℝ) / 2 ^ (205 : ℕ) := by
  let D := gaugeSecant (11 / 2) x center
  let u := x - center
  let r := preconditioner *ᵥ gaugeForce (n := 11) (11 / 2) center
  have hm : center ∈ box := Metric.mem_closedBall_self radius_pos.le
  have hf := gaugeForce_sub_eq_secant (11 / 2) x center (injective_on_box hx) (injective_on_box hm)
  rw [(stationary_iff_on_box hx).2 hs, zero_sub] at hf
  have he : (1 - preconditioner * D) *ᵥ u = u + r := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, ← Matrix.mulVec_mulVec]
    change u - preconditioner *ᵥ (D *ᵥ u) = u + r
    rw [show D *ᵥ u = -gaugeForce (n := 11) (11 / 2) center from hf.symm, Matrix.mulVec_neg, sub_neg_eq_add]
  have hn := dist_mulVec_zero_le (by norm_num) (uniform_secant_rows hx hm) u
  rw [he] at hn
  have ht : dist u 0 ≤ dist (u + r) 0 + dist r 0 := by
    have h := norm_add_le (u + r) (-r)
    simpa only [add_neg_cancel_right, norm_neg, dist_zero_right] using h
  have hres : dist r 0 ≤ (1 : ℝ) / 2 ^ (206 : ℕ) := midpoint_preconditioned_small
  have hh : dist u 0 ≤ (1 : ℝ) / 2 ^ (205 : ℕ) := by linarith only [ht, hn, hres]
  simpa only [dist_eq_norm, sub_zero, u] using hh

theorem isolated_mass : mass (gaugeConfiguration isolatedCoordinates) = 12 :=
  mass_twelve_of_stationary (injective_on_box isolatedCoordinates_spec.1) isolatedCoordinates_spec.2

end Mordell.TwelveRootBox
