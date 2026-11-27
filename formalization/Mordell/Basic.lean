import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Push
import Mathlib.Tactic.Ring

/-! The squared-distance product and unordered pair sums. -/

open scoped BigOperators

noncomputable section

namespace Mordell

abbrev Configuration (n : ℕ) := Fin n → ℂ

def mass {n : ℕ} (z : Configuration n) : ℝ :=
  ∑ i, Complex.normSq (z i)

def distanceSq {n : ℕ} (z : Configuration n) (i j : Fin n) : ℝ :=
  Complex.normSq (z i - z j)

def distanceProduct {n : ℕ} (z : Configuration n) : ℝ :=
  ∏ i, ∏ j, if i < j then distanceSq z i j else 1

def IsMaximizer {n : ℕ} (q : ℝ) (z : Configuration n) : Prop :=
  mass z = q ∧ ∀ w : Configuration n, mass w = q → distanceProduct w ≤ distanceProduct z

def pairSum {n : ℕ} (f : Fin n → Fin n → ℝ) : ℝ :=
  ∑ i, ∑ j, if i < j then f i j else 0

theorem distanceSq_nonneg {n : ℕ} (z : Configuration n) (i j : Fin n) :
    0 ≤ distanceSq z i j :=
  Complex.normSq_nonneg _

theorem distanceSq_pos {n : ℕ} {z : Configuration n}
    (hz : Function.Injective z) {i j : Fin n} (hij : i ≠ j) :
    0 < distanceSq z i j := by
  exact Complex.normSq_pos.mpr (sub_ne_zero.mpr (hz.ne hij))

theorem distanceProduct_nonneg {n : ℕ} (z : Configuration n) :
    0 ≤ distanceProduct z := by
  classical
  unfold distanceProduct
  refine Finset.prod_nonneg fun i _ => Finset.prod_nonneg fun j _ => ?_
  split_ifs
  · exact distanceSq_nonneg z i j
  · norm_num

theorem distanceProduct_pos {n : ℕ} {z : Configuration n} (hz : Function.Injective z) :
    0 < distanceProduct z := by
  classical
  unfold distanceProduct
  refine Finset.prod_pos fun i _ => Finset.prod_pos fun j _ => ?_
  split_ifs with hij
  · exact distanceSq_pos hz hij.ne
  · norm_num

theorem injective_of_distanceProduct_pos {n : ℕ} {z : Configuration n}
    (hz : 0 < distanceProduct z) : Function.Injective z := by
  intro i j he
  by_contra hne
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · have hzero : distanceSq z i j = 0 := by simp [distanceSq, he]
    have hrow : (∏ k : Fin n, if i < k then distanceSq z i k else 1) = 0 := by
      apply Finset.prod_eq_zero (Finset.mem_univ j)
      simp [hlt, hzero]
    have : distanceProduct z = 0 := Finset.prod_eq_zero (Finset.mem_univ i) hrow
    exact hz.ne' this
  · have hzero : distanceSq z j i = 0 := by simp [distanceSq, he]
    have hrow : (∏ k : Fin n, if j < k then distanceSq z j k else 1) = 0 := by
      apply Finset.prod_eq_zero (Finset.mem_univ i)
      simp [hgt, hzero]
    have : distanceProduct z = 0 := Finset.prod_eq_zero (Finset.mem_univ j) hrow
    exact hz.ne' this

theorem log_distanceProduct {n : ℕ} {z : Configuration n}
    (hz : Function.Injective z) :
    Real.log (distanceProduct z) = pairSum (fun i j => Real.log (distanceSq z i j)) := by
  classical
  have hp (i j : Fin n) : 0 < (if i < j then distanceSq z i j else 1) := by
    split_ifs with h
    · exact distanceSq_pos hz h.ne
    · norm_num
  unfold distanceProduct pairSum
  rw [Real.log_prod (fun i _ => (Finset.prod_pos fun j _ => hp i j).ne')]
  apply Finset.sum_congr rfl
  intro i _
  rw [Real.log_prod (fun j _ => (hp i j).ne')]
  apply Finset.sum_congr rfl
  intro j _
  split_ifs <;> simp [*]

theorem pairSum_mono {n : ℕ} {f g : Fin n → Fin n → ℝ}
    (h : ∀ i j, i < j → f i j ≤ g i j) :
    pairSum f ≤ pairSum g := by
  classical
  unfold pairSum
  refine Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => ?_
  split_ifs with hij
  · exact h i j hij
  · exact le_rfl

theorem pairSum_add {n : ℕ} (f g : Fin n → Fin n → ℝ) :
    pairSum (fun i j => f i j + g i j) = pairSum f + pairSum g := by
  classical
  simp only [pairSum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  split_ifs <;> simp

theorem pairSum_sub {n : ℕ} (f g : Fin n → Fin n → ℝ) :
    pairSum (fun i j => f i j - g i j) = pairSum f - pairSum g := by
  classical
  simp only [pairSum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  split_ifs <;> simp

theorem twice_pairSum {n : ℕ} (f : Fin n → Fin n → ℝ) (hf : ∀ i j, f i j = f j i) :
    2 * pairSum f = (∑ i, ∑ j, f i j) - ∑ i, f i i := by
  classical
  have hswap :
      (∑ i : Fin n, ∑ j : Fin n, if j < i then f i j else 0) = pairSum f := by
    rw [Finset.sum_comm]
    unfold pairSum
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [hf j i]
  have hsplit (i j : Fin n) :
      f i j = (if i < j then f i j else 0) + (if j < i then f i j else 0) +
        (if j = i then f i i else 0) := by
    rcases lt_trichotomy i j with hh | rfl | hh
    · simp [hh, hh.not_gt, hh.ne', add_zero]
    · simp
    · simp [hh, hh.not_gt, hh.ne, add_zero]
  conv_rhs =>
    lhs
    arg 2
    ext i
    arg 2
    ext j
    rw [hsplit i j]
  simp only [Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  rw [hswap]
  unfold pairSum
  ring

end Mordell
