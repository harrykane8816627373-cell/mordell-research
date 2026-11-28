import Mordell.Basic
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Tactic.GCongr

/-! A rational congruence and a strict row bound certify positive definiteness throughout a matrix box. -/

open scoped BigOperators Matrix

namespace Mordell

private theorem posDef_of_row_bounds {n : ℕ} {A : Matrix (Fin n) (Fin n) ℝ}
    (hs : A.IsHermitian) (hd : ∀ i, (∑ j, if i = j then 0 else |A i j|) < A i i) :
    A.PosDef := by
  let b (i j : Fin n) := if i = j then 0 else |A i j|
  have hsym (i j : Fin n) : b i j = b j i := by
    have h : A j i = A i j := hs.apply i j
    simp [b, eq_comm, h]
  apply Matrix.PosDef.of_dotProduct_mulVec_pos hs
  intro x hx
  have hpt (i j : Fin n) : (if i = j then A i i * x i ^ 2 else 0) -
      b i j * (x i ^ 2 + x j ^ 2) / 2 ≤ x i * (A i j * x j) := by
    by_cases hij : i = j
    · subst j
      simp only [b, if_true, zero_mul, zero_div, sub_zero]
      ring_nf
      rfl
    · simp only [b, if_neg hij, zero_sub]
      have hab : |x i * x j| ≤ (x i ^ 2 + x j ^ 2) / 2 := by
        rw [abs_mul]
        nlinarith [sq_nonneg (|x i| - |x j|), sq_abs (x i), sq_abs (x j)]
      calc
        -(|A i j| * (x i ^ 2 + x j ^ 2) / 2) ≤ -(|A i j| * |x i * x j|) := by
          rw [mul_div_assoc]
          exact neg_le_neg (mul_le_mul_of_nonneg_left hab (abs_nonneg _))
        _ = -|x i * (A i j * x j)| := by simp only [abs_mul]; ring
        _ ≤ _ := neg_abs_le _
  have hsum := Finset.sum_le_sum (fun i (_ : i ∈ Finset.univ) =>
    Finset.sum_le_sum (fun j (_ : j ∈ Finset.univ) => hpt i j))
  have hdiag : (∑ i, ∑ j, if i = j then A i i * x i ^ 2 else 0) = ∑ i, A i i * x i ^ 2 := by
    simp
  have hsplit : (∑ i, ∑ j, b i j * (x i ^ 2 + x j ^ 2) / 2) =
      ∑ i, (∑ j, b i j) * x i ^ 2 := by
    have hswap : (∑ i, ∑ j, b i j * x j ^ 2) = ∑ i, ∑ j, b i j * x i ^ 2 := by
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by rw [hsym j i]
    calc
      _ = (∑ i, ∑ j, b i j * x i ^ 2) * 2⁻¹ + (∑ i, ∑ j, b i j * x j ^ 2) * 2⁻¹ := by
        simp only [mul_add, div_eq_mul_inv, add_mul, Finset.sum_add_distrib, Finset.sum_mul]
      _ = _ := by
        rw [hswap]
        simp only [← Finset.sum_mul]
        ring
  simp only [Finset.sum_sub_distrib, hdiag, hsplit] at hsum
  have hq : (∑ i, (A i i - ∑ j, b i j) * x i ^ 2) ≤ star x ⬝ᵥ (A *ᵥ x) := by
    simpa only [sub_mul, Finset.sum_sub_distrib, dotProduct, Matrix.mulVec, Finset.mul_sum] using hsum
  apply lt_of_lt_of_le _ hq
  obtain ⟨i, hi⟩ : ∃ i, x i ≠ 0 := by simpa [funext_iff] using hx
  exact Finset.sum_pos' (fun j _ => mul_nonneg (sub_nonneg.mpr (hd j).le) (sq_nonneg _))
    ⟨i, Finset.mem_univ i, mul_pos (sub_pos.mpr (hd i)) (sq_pos_of_ne_zero hi)⟩

private theorem posDef_of_congruence {n : ℕ} {A P : Matrix (Fin n) (Fin n) ℝ}
    (hs : A.IsHermitian) (hP : P.BlockTriangular id) (hdiag : ∀ i, P i i = 1)
    (hpos : (Pᴴ * A * P).PosDef) : A.PosDef := by
  have hdet : P.det = 1 := by rw [Matrix.det_of_upperTriangular hP]; simp [hdiag]
  have hunit : IsUnit P := (Matrix.isUnit_iff_isUnit_det _).mpr (hdet ▸ isUnit_one)
  apply Matrix.PosDef.of_dotProduct_mulVec_pos hs
  intro x hx
  obtain ⟨y, hy⟩ := Matrix.mulVec_surjective_iff_isUnit.mpr hunit x
  have hy0 : y ≠ 0 := by
    intro he
    apply hx
    simpa [he] using hy.symm
  have hh := hpos.dotProduct_mulVec_pos hy0
  have hform : star y ⬝ᵥ ((Pᴴ * A * P) *ᵥ y) = star (P *ᵥ y) ⬝ᵥ (A *ᵥ (P *ᵥ y)) := by
    simp only [Matrix.star_mulVec, Matrix.dotProduct_mulVec, Matrix.vecMul_vecMul]
  rwa [hform, hy] at hh

structure GramCertificate (n : ℕ) where
  lower : Fin n → Fin n → ℚ
  upper : Fin n → Fin n → ℚ
  radius : ℚ
  upperEntries : Fin n → Fin n → ℚ

def gramFromDyadic {n : ℕ} (endpointBits basisBits : ℕ) (radius : ℤ)
    (lower upper basis : Fin n → Fin n → ℤ) : GramCertificate n where
  lower i j := (lower i j : ℚ) / 2 ^ endpointBits
  upper i j := (upper i j : ℚ) / 2 ^ endpointBits
  radius := (radius : ℚ) / 2 ^ endpointBits
  upperEntries i j := (basis i j : ℚ) / 2 ^ basisBits

namespace GramCertificate

def basis {n : ℕ} (c : GramCertificate n) (i j : Fin n) : ℚ :=
  if i < j then c.upperEntries i j else if i = j then 1 else 0

def center {n : ℕ} (c : GramCertificate n) (i j : Fin n) : ℚ := (c.lower i j + c.upper i j) / 2

private def sumRat {n : ℕ} (f : Fin n → ℚ) : ℚ := (List.ofFn f).sum

private theorem sumRat_eq {n : ℕ} (f : Fin n → ℚ) : sumRat f = ∑ i, f i := by
  simp [sumRat, List.sum_ofFn]

def colSum {n : ℕ} (c : GramCertificate n) (i : Fin n) : ℚ := sumRat fun k => |c.basis k i|

def entry {n : ℕ} (c : GramCertificate n) (i j : Fin n) : ℚ :=
  sumRat fun k => sumRat fun l => c.basis k i * c.center k l * c.basis l j

def widthsOK {n : ℕ} (c : GramCertificate n) : Prop :=
  0 ≤ c.radius ∧ ∀ i j, c.upper i j - c.lower i j ≤ 2 * c.radius

instance {n : ℕ} (c : GramCertificate n) : Decidable c.widthsOK :=
  inferInstanceAs (Decidable (0 ≤ c.radius ∧ ∀ i j, c.upper i j - c.lower i j ≤ 2 * c.radius))

def rowOK {n : ℕ} (c : GramCertificate n) (i : Fin n) : Prop :=
  (sumRat fun j => if i = j then 0 else |c.entry i j|) +
    c.radius * c.colSum i * (sumRat c.colSum) < c.entry i i

instance {n : ℕ} (c : GramCertificate n) (i : Fin n) : Decidable (c.rowOK i) :=
  inferInstanceAs (Decidable ((sumRat fun j => if i = j then 0 else |c.entry i j|) +
    c.radius * c.colSum i * (sumRat c.colSum) < c.entry i i))

/-- All real symmetric matrices enclosed by a checked rational box are positive definite. -/
theorem sound {n : ℕ} (c : GramCertificate n) (hw : c.widthsOK) (hr : ∀ i, c.rowOK i)
    {A : Matrix (Fin n) (Fin n) ℝ} (hs : A.IsHermitian)
    (hbox : ∀ i j, (c.lower i j : ℝ) ≤ A i j ∧ A i j ≤ (c.upper i j : ℝ)) : A.PosDef := by
  let P : Matrix (Fin n) (Fin n) ℝ := fun i j => c.basis i j
  let R : Matrix (Fin n) (Fin n) ℝ := Pᴴ * A * P
  let e : ℝ := c.radius
  let d (i : Fin n) := ∑ k, |P k i|
  have hd' (i : Fin n) : d i = c.colSum i := by
    simp [d, colSum, sumRat_eq, Rat.cast_sum, Rat.cast_abs, P]
  have hcenter (k l : Fin n) : |A k l - (c.center k l : ℝ)| ≤ e := by
    have hb := hbox k l
    have hwidth : (c.upper k l : ℝ) - c.lower k l ≤ 2 * (c.radius : ℝ) := by
      exact_mod_cast hw.2 k l
    rw [abs_le]
    push_cast [center]
    exact ⟨by linarith, by linarith⟩
  have hentry (i j : Fin n) : |R i j - (c.entry i j : ℝ)| ≤ e * d i * d j := by
    have hdiff : R i j - (c.entry i j : ℝ) =
        ∑ k, ∑ l, P k i * (A k l - (c.center k l : ℝ)) * P l j := by
      have hR : R i j = ∑ k, ∑ l, P k i * A k l * P l j := by
        simp only [R, Matrix.mul_apply, Matrix.conjTranspose_apply, star_trivial, Finset.sum_mul]
        exact Finset.sum_comm
      rw [hR]
      simp [entry, sumRat_eq, Rat.cast_sum, ← Finset.sum_sub_distrib, mul_sub, sub_mul, P]
    rw [hdiff]
    calc
      _ ≤ ∑ k, ∑ l, |P k i * (A k l - (c.center k l : ℝ)) * P l j| :=
        (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun _ _ => Finset.abs_sum_le_sum_abs _ _))
      _ ≤ ∑ k, ∑ l, e * |P k i| * |P l j| := by
        refine Finset.sum_le_sum fun k _ => Finset.sum_le_sum fun l _ => ?_
        simp only [abs_mul]
        rw [mul_comm e]
        gcongr
        exact hcenter k l
      _ = _ := by
        simp only [d, Finset.mul_sum, Finset.sum_mul]
        exact Finset.sum_comm
  apply posDef_of_congruence hs (P := P)
  · intro i j h
    change j < i at h
    simp [P, basis, not_lt_of_ge (le_of_lt h), ne_of_gt h]
  · intro i
    simp [P, basis]
  · apply posDef_of_row_bounds (Matrix.isHermitian_conjTranspose_mul_mul P hs)
    intro i
    have hrow (j : Fin n) : (if i = j then 0 else |R i j|) ≤
        (if i = j then 0 else |(c.entry i j : ℝ)|) + (if i = j then 0 else e * d i * d j) := by
      split_ifs with hij
      · simp
      · calc
          _ = |(R i j - (c.entry i j : ℝ)) + (c.entry i j : ℝ)| := by ring_nf
          _ ≤ _ := (abs_add_le _ _).trans (by linarith [hentry i j])
    have hrow' := Finset.sum_le_sum fun j (_ : j ∈ Finset.univ) => hrow j
    rw [Finset.sum_add_distrib] at hrow'
    have htest : (∑ j, if i = j then 0 else |(c.entry i j : ℝ)|) + e * d i * (∑ j, d j) <
        (c.entry i i : ℝ) := by
      have h := (Rat.cast_lt (K := ℝ)).2 (hr i)
      have hite (j : Fin n) : ((if i = j then 0 else |c.entry i j| : ℚ) : ℝ) =
          if i = j then 0 else |(c.entry i j : ℝ)| := by split_ifs <;> simp
      simpa [rowOK, sumRat_eq, Rat.cast_sum, Rat.cast_add, Rat.cast_mul, hite, ← hd', e] using h
    have herr : (∑ j, if i = j then 0 else e * d i * d j) ≤
        e * d i * (∑ j, d j) - e * d i * d i := by
      have h : (∑ j, if i = j then e * d i * d j else 0) = e * d i * d i := by simp
      apply le_of_eq
      calc
        _ = ∑ j, (e * d i * d j - if i = j then e * d i * d j else 0) :=
          Finset.sum_congr rfl fun j _ => by split_ifs <;> simp
        _ = _ := by rw [Finset.sum_sub_distrib, h, ← Finset.mul_sum]
    change (∑ j, if i = j then 0 else |R i j|) < R i i
    linarith [neg_le_of_abs_le (hentry i i), herr, hrow', htest]

end GramCertificate
end Mordell
