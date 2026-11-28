import Mordell.Kernel
import Mordell.PairGeometry
import Mathlib.LinearAlgebra.Matrix.PosDef

/-! Polynomial angular factors and positive definite blocks of the two-point kernel. -/

open scoped BigOperators ComplexConjugate Matrix

noncomputable section

namespace Mordell

def harmonicPolynomial (uv t : ℝ) : ℕ → ℝ
  | 0 => 1
  | 1 => t
  | k + 2 => 2 * t * harmonicPolynomial uv t (k + 1) - uv * harmonicPolynomial uv t k

theorem harmonicPolynomial_eq_re_pow (w : ℂ) (k : ℕ) :
    harmonicPolynomial (Complex.normSq w) w.re k = (w ^ k).re := by
  have hw : w ^ 2 = ((2 * w.re : ℝ) : ℂ) * w - (Complex.normSq w : ℂ) := by
    apply Complex.ext <;>
      simp only [pow_two, Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im,
        Complex.ofReal_re, Complex.ofReal_im, Complex.normSq_apply] <;> ring
  induction k using Nat.twoStepInduction with
  | zero => simp [harmonicPolynomial]
  | one => simp [harmonicPolynomial]
  | more k h₀ h₁ =>
    rw [harmonicPolynomial, h₀, h₁, show w ^ (k + 2) = w ^ k * w ^ 2 by rw [pow_add], hw, mul_sub]
    rw [show w ^ k * ((↑(2 * w.re) : ℂ) * w) = (↑(2 * w.re) : ℂ) * w ^ (k + 1) by
      rw [pow_succ]; ring, mul_comm (w ^ k) (Complex.normSq w : ℂ)]
    simp

theorem harmonicPolynomial_pair (z w : ℂ) (k : ℕ) :
    harmonicPolynomial (Complex.normSq z * Complex.normSq w) (realDot z w) k =
      realDot (z ^ k) (w ^ k) := by
  rw [← Complex.normSq_conj w, ← Complex.normSq_mul, realDot, harmonicPolynomial_eq_re_pow,
    realDot, mul_pow, map_pow]

def matrixForm {κ : Type*} [Fintype κ] (Q : Matrix κ κ ℝ) (x y : κ → ℝ) : ℝ :=
  ∑ a, ∑ b, Q a b * (x a * y b)

private theorem matrixForm_diagonal {κ : Type*} [Fintype κ] (Q : Matrix κ κ ℝ) (x : κ → ℝ) :
    matrixForm Q x x = x ⬝ᵥ (Q *ᵥ x) := by
  simp only [matrixForm, dotProduct, Matrix.mulVec, Finset.mul_sum]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

private theorem sum_matrixForm {κ : Type*} [Fintype κ] {n : ℕ}
    (Q : Matrix κ κ ℝ) (x : Fin n → κ → ℝ) :
    (∑ i, ∑ j, matrixForm Q (x i) (x j)) = matrixForm Q (fun a => ∑ i, x i a) (fun a => ∑ i, x i a) := by
  unfold matrixForm
  calc
    _ = ∑ i, ∑ a, ∑ b, ∑ j, Q a b * (x i a * x j b) := by
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ a, ∑ b, ∑ i, ∑ j, Q a b * (x i a * x j b) := by
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = _ := by
      refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
      simp only [Finset.mul_sum, Finset.sum_mul]
      exact Finset.sum_comm

def matrixGram {κ : Type*} [Fintype κ] {n : ℕ} (Q : Matrix κ κ ℝ)
    (f : Fin n → κ → ℂ) (i j : Fin n) : ℝ :=
  matrixForm Q (fun a => (f i a).re) (fun a => (f j a).re) +
    matrixForm Q (fun a => (f i a).im) (fun a => (f j a).im)

theorem sum_matrixGram {κ : Type*} [Fintype κ] {n : ℕ} (Q : Matrix κ κ ℝ)
    (f : Fin n → κ → ℂ) :
    (∑ i, ∑ j, matrixGram Q f i j) =
      matrixForm Q (fun a => (∑ i, f i a).re) (fun a => (∑ i, f i a).re) +
      matrixForm Q (fun a => (∑ i, f i a).im) (fun a => (∑ i, f i a).im) := by
  simp only [matrixGram, Finset.sum_add_distrib, sum_matrixForm]
  have hRe (a : κ) : (∑ i, f i a).re = ∑ i, (f i a).re :=
    map_sum Complex.reAddGroupHom (fun i => f i a) Finset.univ
  have hIm (a : κ) : (∑ i, f i a).im = ∑ i, (f i a).im :=
    map_sum Complex.imAddGroupHom (fun i => f i a) Finset.univ
  simp only [hRe, hIm]

theorem matrixGram_symm {κ : Type*} [Fintype κ] {n : ℕ} {Q : Matrix κ κ ℝ}
    (hQ : Q.IsHermitian) (f : Fin n → κ → ℂ) (i j : Fin n) :
    matrixGram Q f i j = matrixGram Q f j i := by
  have hs (x y : κ → ℝ) : matrixForm Q x y = matrixForm Q y x := by
    unfold matrixForm
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    have h : Q b a = Q a b := hQ.apply a b
    rw [h, mul_comm (x b)]
  simp [matrixGram, hs]

/-- Positive definiteness of a real block controls both components of each complex moment. -/
theorem matrixGram_nonneg_and_zero {κ : Type*} [Fintype κ] {n : ℕ} {Q : Matrix κ κ ℝ}
    (hQ : Q.PosDef) (f : Fin n → κ → ℂ) :
    0 ≤ (∑ i, ∑ j, matrixGram Q f i j) ∧
      ((∑ i, ∑ j, matrixGram Q f i j) = 0 → ∀ a, ∑ i, f i a = 0) := by
  let re : κ → ℝ := fun a => ∑ i, (f i a).re
  let im : κ → ℝ := fun a => ∑ i, (f i a).im
  have he : (∑ i, ∑ j, matrixGram Q f i j) = re ⬝ᵥ (Q *ᵥ re) + im ⬝ᵥ (Q *ᵥ im) := by
    simp only [matrixGram, Finset.sum_add_distrib, sum_matrixForm, matrixForm_diagonal, re, im]
  have hr : 0 ≤ re ⬝ᵥ (Q *ᵥ re) := hQ.posSemidef.dotProduct_mulVec_nonneg re
  have hi : 0 ≤ im ⬝ᵥ (Q *ᵥ im) := hQ.posSemidef.dotProduct_mulVec_nonneg im
  rw [he]
  refine ⟨add_nonneg hr hi, fun hz a => ?_⟩
  have hr0 : re = 0 := by
    by_contra hn
    have h : 0 < re ⬝ᵥ (Q *ᵥ re) := hQ.dotProduct_mulVec_pos hn
    linarith
  have hi0 : im = 0 := by
    by_contra hn
    have h : 0 < im ⬝ᵥ (Q *ᵥ im) := hQ.dotProduct_mulVec_pos hn
    linarith
  apply Complex.ext
  · exact (map_sum Complex.reAddGroupHom (fun i => f i a) Finset.univ).trans (congrFun hr0 a)
  · exact (map_sum Complex.imAddGroupHom (fun i => f i a) Finset.univ).trans (congrFun hi0 a)

theorem matrixGram_harmonic {κ : Type*} [Fintype κ] {n : ℕ} (Q : Matrix κ κ ℝ)
    (z : Configuration n) (r : Fin n → κ → ℝ) (k : ℕ) (i j : Fin n) :
    matrixGram Q (fun l a => (r l a : ℂ) * z l ^ k) i j =
      harmonicPolynomial (Complex.normSq (z i) * Complex.normSq (z j)) (realDot (z i) (z j)) k *
        matrixForm Q (r i) (r j) := by
  rw [harmonicPolynomial_pair]
  simp only [matrixGram, matrixForm, realDot, Complex.mul_re, Complex.conj_re, Complex.conj_im,
    Complex.ofReal_re, Complex.ofReal_im, zero_mul, sub_zero, Complex.mul_im, add_zero,
    sub_neg_eq_add, mul_neg, Finset.mul_sum, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

def blockGram {ι : Type*} [Fintype ι] {κ : ι → Type*} [∀ k, Fintype (κ k)] {n : ℕ}
    (Q : ∀ k, Matrix (κ k) (κ k) ℝ) (f : ∀ k, Fin n → κ k → ℂ) (i j : Fin n) : ℝ :=
  ∑ k, matrixGram (Q k) (f k) i j

/-- Each positive definite block contributes a nonnegative term. Equality
    of the pair bound forces all complex feature moments to vanish. -/
theorem blockGram_majorant {ι : Type*} [Fintype ι] {κ : ι → Type*} [∀ k, Fintype (κ k)]
    {n : ℕ} {z : Configuration n} (hn : 1 < n) (hz : Function.Injective z)
    (hq : mass z = (n : ℝ) - 1) (B : ℝ) (Q : ∀ k, Matrix (κ k) (κ k) ℝ)
    (hQ : ∀ k, (Q k).PosDef) (f : ∀ k, Fin n → κ k → ℂ) (W V : Fin n → Fin n → ℝ)
    (hW : pairSum W = 0) (hV : 0 ≤ pairSum V)
    (hmajor : ∀ i j, i < j → Real.log (distanceSq z i j) ≤
      baseKernel B (fun k => Complex.normSq (z k)) (blockGram Q f) i j + W i j - V i j) :
    Real.log (distanceProduct z) ≤ B ∧
      (Real.log (distanceProduct z) = B → ∀ k a, ∑ i, f k i a = 0) := by
  have hs (i j : Fin n) : blockGram Q f i j = blockGram Q f j i :=
    Finset.sum_congr rfl (fun k _ => matrixGram_symm (hQ k).isHermitian (f k) i j)
  have he : (∑ i, ∑ j, blockGram Q f i j) = ∑ k, ∑ i, ∑ j, matrixGram (Q k) (f k) i j := by
    unfold blockGram
    calc
      _ = ∑ i, ∑ k, ∑ j, matrixGram (Q k) (f k) i j :=
        Finset.sum_congr rfl (fun _ _ => Finset.sum_comm)
      _ = _ := Finset.sum_comm
  have hpos (k : ι) : 0 ≤ ∑ i, ∑ j, matrixGram (Q k) (f k) i j :=
    (matrixGram_nonneg_and_zero (hQ k) (f k)).1
  have hp : 0 ≤ ∑ i, ∑ j, blockGram Q f i j := by
    rw [he]
    exact Finset.sum_nonneg fun k _ => hpos k
  obtain ⟨hb, hzero⟩ := kernel_majorant hn hz hq B _ W V hs hp hW hV hmajor
  refine ⟨hb, fun h k => (matrixGram_nonneg_and_zero (hQ k) (f k)).2 ?_⟩
  have ht := hzero h
  rw [he] at ht
  exact (Finset.sum_eq_zero_iff_of_nonneg (fun k _ => hpos k)).mp ht k (Finset.mem_univ k)

structure BlockBound {n : ℕ} (z w : Configuration n) where
  blockCount : ℕ
  width : Fin blockCount → ℕ
  weights : ∀ k, Matrix (Fin (width k)) (Fin (width k)) ℝ
  positive : ∀ k, (weights k).PosDef
  feature : ∀ k, Fin n → Fin (width k) → ℂ
  ward : Fin n → Fin n → ℝ
  variation : Fin n → Fin n → ℝ
  ward_zero : pairSum ward = 0
  variation_nonneg : 0 ≤ pairSum variation
  majorizes : ∀ i j, i < j → Real.log (distanceSq z i j) ≤
    baseKernel (Real.log (distanceProduct w)) (fun k => Complex.normSq (z k))
      (blockGram weights feature) i j + ward i j - variation i j

theorem BlockBound.log_le_and_moments_zero {n : ℕ} {z w : Configuration n}
    (cert : BlockBound z w) (hn : 1 < n) (hz : Function.Injective z) (hq : mass z = (n : ℝ) - 1) :
    Real.log (distanceProduct z) ≤ Real.log (distanceProduct w) ∧
      (Real.log (distanceProduct z) = Real.log (distanceProduct w) →
        ∀ k a, ∑ i, cert.feature k i a = 0) :=
  blockGram_majorant hn hz hq _ cert.weights cert.positive cert.feature cert.ward cert.variation
    cert.ward_zero cert.variation_nonneg cert.majorizes

end Mordell
