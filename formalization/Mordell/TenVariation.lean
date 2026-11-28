import Mordell.Variational
import Mordell.HarmonicKernel
import Mathlib.Analysis.Matrix.PosDef

/-! The two radial second-variation blocks of the ten-point certificate. -/

open scoped BigOperators ComplexConjugate Matrix

noncomputable section

namespace Mordell

def tenRadialVariation (g : ℝ → ℝ) (sign u v t : ℝ) : ℝ :=
  (5 / 9) * (u * g u ^ 2 + v * g v ^ 2) +
    sign * (((u - t) * g u + (v - t) * g v) ^ 2 - (u * v - t ^ 2) * (g u - g v) ^ 2) /
      (u + v - 2 * t) ^ 2

private theorem radialQuotient_sq_re (z w : ℂ) (a b : ℝ) :
    ((((a : ℂ) * z - (b : ℂ) * w) / (z - w)) ^ 2).re =
      (((Complex.normSq z - realDot z w) * a + (Complex.normSq w - realDot z w) * b) ^ 2 -
        (Complex.normSq z * Complex.normSq w - realDot z w ^ 2) * (a - b) ^ 2) /
        (Complex.normSq z + Complex.normSq w - 2 * realDot z w) ^ 2 := by
  have hRe : (((a : ℂ) * z - (b : ℂ) * w) / (z - w)).re =
      ((Complex.normSq z - realDot z w) * a + (Complex.normSq w - realDot z w) * b) /
        Complex.normSq (z - w) := by
    rw [Complex.div_re, ← add_div]
    congr 1
    simp only [Complex.normSq_apply, realDot, Complex.mul_re, Complex.mul_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.sub_re, Complex.sub_im, Complex.conj_re, Complex.conj_im,
      zero_mul, sub_zero, add_zero]
    ring
  have hIm : (((a : ℂ) * z - (b : ℂ) * w) / (z - w)).im =
      (b - a) * (z * conj w).im / Complex.normSq (z - w) := by
    rw [Complex.div_im, ← sub_div]
    congr 1
    simp only [Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
      Complex.sub_re, Complex.sub_im, Complex.conj_re, Complex.conj_im, zero_mul, sub_zero, add_zero]
    ring
  have hGram := Complex.normSq_mul z (conj w)
  rw [Complex.normSq_conj, Complex.normSq_apply] at hGram
  have hGram' : (z * conj w).im ^ 2 = Complex.normSq z * Complex.normSq w - realDot z w ^ 2 := by
    unfold realDot
    linear_combination hGram
  have hDist : Complex.normSq (z - w) = Complex.normSq z + Complex.normSq w - 2 * realDot z w :=
    Complex.normSq_sub _ _
  rw [pow_two, Complex.mul_re, ← pow_two, ← pow_two, hRe, hIm, div_pow, div_pow, mul_pow,
    hGram', ← sub_div, hDist]
  ring

private theorem pairSum_tenVariationNorm (u a : Fin 10 → ℝ) :
    pairSum (fun i j => (5 / 9) * (u i * a i ^ 2 + u j * a j ^ 2)) = 5 * ∑ i, u i * a i ^ 2 := by
  let l (i : Fin 10) := (5 / 9) * (u i * a i ^ 2)
  have he : (fun i j => (5 / 9) * (u i * a i ^ 2 + u j * a j ^ 2)) = fun i j => l i + l j := by
    funext i j
    exact mul_add _ _ _
  rw [he]
  have h := twice_pairSum (fun i j => l i + l j) (fun _ _ => add_comm _ _)
  simp only [Finset.sum_add_distrib, Finset.sum_const, ← Finset.mul_sum, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul] at h
  have hl : (∑ i, l i) = (5 / 9) * ∑ i, u i * a i ^ 2 := (Finset.mul_sum _ _ _).symm
  rw [hl] at h
  norm_num at h
  linear_combination h / 2

/-- Each signed scalar radial second-variation term has nonnegative edge sum. -/
theorem tenRadialVariation_pairSum_nonneg {z : Configuration 10} (hz : IsMaximizer 9 z)
    (g : ℝ → ℝ) {sign : ℝ} (hs : sign = 1 ∨ sign = -1) :
    0 ≤ pairSum (fun i j => tenRadialVariation g sign (Complex.normSq (z i))
      (Complex.normSq (z j)) (realDot (z i) (z j))) := by
  let u (i : Fin 10) := Complex.normSq (z i)
  let ξ (i : Fin 10) := (g (u i) : ℂ) * z i
  let η (i : Fin 10) := Complex.I * ξ i
  have hcount : pairCount 10 = 45 := by decide +kernel
  have hpos := secondVariation_nonneg (by decide : 1 < 10) (by norm_num : (0 : ℝ) < 9) hz ξ
  have hneg := secondVariation_nonneg (by decide : 1 < 10) (by norm_num : (0 : ℝ) < 9) hz η
  norm_num [hcount] at hpos hneg
  have hmass : mass ξ = ∑ i, u i * g (u i) ^ 2 := by
    simp only [mass, ξ, Complex.normSq_mul, Complex.normSq_ofReal, sq, u]
    exact Finset.sum_congr rfl fun _ _ => by ring
  have hrot : mass η = mass ξ := by
    simp [mass, η, Complex.normSq_mul]
  rw [hrot, hmass] at hneg
  rw [hmass] at hpos
  have hquot (i j : Fin 10) :
      (((η i - η j) / (z i - z j)) ^ 2).re = -(((ξ i - ξ j) / (z i - z j)) ^ 2).re := by
    simp only [η, ← mul_sub, mul_div_assoc, mul_pow, Complex.I_sq, neg_one_mul, Complex.neg_re]
  rw [show (fun i j => (((η i - η j) / (z i - z j)) ^ 2).re) =
      (fun i j => (-1 : ℝ) * (((ξ i - ξ j) / (z i - z j)) ^ 2).re) by
        funext i j
        simpa only [neg_one_mul] using hquot i j] at hneg
  have hsplit : pairSum (fun i j => tenRadialVariation g sign (u i) (u j) (realDot (z i) (z j))) =
      5 * (∑ i, u i * g (u i) ^ 2) +
        pairSum (fun i j => sign * (((ξ i - ξ j) / (z i - z j)) ^ 2).re) := by
    simp only [tenRadialVariation, ξ, radialQuotient_sq_re, mul_div_assoc, pairSum_add, pairSum_tenVariationNorm, u]
  change 0 ≤ pairSum (fun i j => tenRadialVariation g sign (u i) (u j) (realDot (z i) (z j)))
  rw [hsplit]
  rcases hs with rfl | rfl
  · simpa only [one_mul] using hpos
  · exact hneg

def tenVariationBlock {κ : Type*} [Fintype κ] (Q : Matrix κ κ ℝ) (f : ℝ → κ → ℝ)
    (sign u v t : ℝ) : ℝ :=
  let X i := (u - t) * f u i + (v - t) * f v i
  let Y i := f u i - f v i
  (5 / 9) * (u * matrixForm Q (f u) (f u) + v * matrixForm Q (f v) (f v)) +
    sign * (matrixForm Q X X - (u * v - t ^ 2) * matrixForm Q Y Y) / (u + v - 2 * t) ^ 2

private theorem matrixForm_conjTranspose_mul_self {κ : Type*} [Fintype κ]
    (R : Matrix κ κ ℝ) (x : κ → ℝ) :
    matrixForm (Rᴴ * R) x x = ∑ k, (∑ a, R k a * x a) ^ 2 := by
  unfold matrixForm
  simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, star_trivial, Finset.sum_mul]
  calc
    _ = ∑ a, ∑ k, ∑ b, (R k a * R k b) * (x a * x b) := by
      exact Finset.sum_congr rfl fun _ _ => Finset.sum_comm
    _ = ∑ k, ∑ a, ∑ b, (R k a * R k b) * (x a * x b) := Finset.sum_comm
    _ = _ := by
      simp only [pow_two, Finset.sum_mul, Finset.mul_sum]
      exact Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => by
        ring

private theorem tenVariationBlock_conjTranspose_mul_self {κ : Type*} [Fintype κ]
    (R : Matrix κ κ ℝ) (f : ℝ → κ → ℝ) (sign u v t : ℝ) :
    tenVariationBlock (Rᴴ * R) f sign u v t =
      ∑ k, tenRadialVariation (fun r => ∑ a, R k a * f r a) sign u v t := by
  let g (k : κ) (r : ℝ) := ∑ a, R k a * f r a
  have hX (k : κ) :
      (∑ a, R k a * ((u - t) * f u a + (v - t) * f v a)) =
        (u - t) * g k u + (v - t) * g k v := by
    simp only [g, mul_add, Finset.sum_add_distrib, Finset.mul_sum]
    congr 1 <;> exact Finset.sum_congr rfl fun _ _ => by ring
  have hY (k : κ) : (∑ a, R k a * (f u a - f v a)) = g k u - g k v := by
    simp only [g, mul_sub, Finset.sum_sub_distrib]
  change tenVariationBlock (Rᴴ * R) f sign u v t = ∑ k, tenRadialVariation (g k) sign u v t
  simp only [tenVariationBlock, tenRadialVariation, matrixForm_conjTranspose_mul_self, hX, hY, g]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum, mul_add]
  rw [← Finset.sum_div]
  simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.mul_sum]

private theorem real_psd_factor {κ : Type*} [Fintype κ] [DecidableEq κ]
    {Q : Matrix κ κ ℝ} (hQ : Q.PosSemidef) :
    ∃ R : Matrix κ κ ℝ, Q = Rᴴ * R := by
  let U : Matrix κ κ ℝ := hQ.isHermitian.eigenvectorUnitary
  let D : Matrix κ κ ℝ := Matrix.diagonal fun i => Real.sqrt (hQ.isHermitian.eigenvalues i)
  refine ⟨D * Uᴴ, ?_⟩
  have hD : Dᴴ * D = Matrix.diagonal hQ.isHermitian.eigenvalues := by
    simp only [D, Matrix.diagonal_conjTranspose, star_trivial, Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    exact Real.mul_self_sqrt (hQ.eigenvalues_nonneg i)
  calc
    Q = U * Matrix.diagonal hQ.isHermitian.eigenvalues * Uᴴ := by
      simpa [U, Unitary.conjStarAlgAut_apply] using hQ.isHermitian.spectral_theorem
    _ = _ := by
      rw [← hD, Matrix.conjTranspose_mul, Matrix.conjTranspose_conjTranspose]
      simp only [mul_assoc]

/-- The archived positive matrix in either signed block gives a nonnegative
    contribution to the variation budget of every mass-nine maximizer. -/
theorem tenVariationBlock_pairSum_nonneg {κ : Type*} [Fintype κ] [DecidableEq κ]
    {z : Configuration 10} (hz : IsMaximizer 9 z) (Q : Matrix κ κ ℝ) (hQ : Q.PosSemidef)
    (f : ℝ → κ → ℝ) {sign : ℝ} (hs : sign = 1 ∨ sign = -1) :
    0 ≤ pairSum (fun i j => tenVariationBlock Q f sign (Complex.normSq (z i))
      (Complex.normSq (z j)) (realDot (z i) (z j))) := by
  obtain ⟨R, rfl⟩ := real_psd_factor hQ
  let term (k : κ) (i j : Fin 10) := tenRadialVariation (fun r => ∑ a, R k a * f r a) sign
    (Complex.normSq (z i)) (Complex.normSq (z j)) (realDot (z i) (z j))
  have he : pairSum (fun i j => ∑ k, term k i j) = ∑ k, pairSum (term k) := by
    unfold pairSum
    have h (i j : Fin 10) : (if i < j then ∑ k, term k i j else 0) =
        ∑ k, if i < j then term k i j else 0 := by split <;> simp_all
    simp only [h]
    calc
      _ = ∑ i, ∑ k, ∑ j, if i < j then term k i j else 0 :=
        Finset.sum_congr rfl fun _ _ => Finset.sum_comm
      _ = _ := Finset.sum_comm
  simp only [tenVariationBlock_conjTranspose_mul_self]
  change 0 ≤ pairSum (fun i j => ∑ k, term k i j)
  rw [he]
  exact Finset.sum_nonneg fun _ _ => tenRadialVariation_pairSum_nonneg hz _ hs

/-- The sum of the two archived signed blocks supplies the variation budget. -/
theorem tenVariation_pairSum_nonneg {κ : Type*} [Fintype κ] [DecidableEq κ]
    {z : Configuration 10} (hz : IsMaximizer 9 z) (Qpos Qneg : Matrix κ κ ℝ)
    (hpos : Qpos.PosSemidef) (hneg : Qneg.PosSemidef) (f : ℝ → κ → ℝ) :
    0 ≤ pairSum (fun i j =>
      tenVariationBlock Qpos f 1 (Complex.normSq (z i)) (Complex.normSq (z j)) (realDot (z i) (z j)) +
      tenVariationBlock Qneg f (-1) (Complex.normSq (z i)) (Complex.normSq (z j)) (realDot (z i) (z j))) := by
  rw [pairSum_add]
  exact add_nonneg (tenVariationBlock_pairSum_nonneg hz Qpos hpos f (Or.inl rfl))
    (tenVariationBlock_pairSum_nonneg hz Qneg hneg f (Or.inr rfl))

end Mordell
