import Mordell.Normalization
import Mathlib.Data.ZMod.Basic

/-! Summation of the two-point majorant, its zero-sum corrections and Gram penalties. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def baseKernel {n : ℕ} (B : ℝ) (u : Fin n → ℝ) (K : Fin n → Fin n → ℝ)
    (i j : Fin n) : ℝ :=
  B / ((n : ℝ) * (n - 1) / 2) - 1 +
    ((n : ℝ) / (2 * (n - 1))) * (u i + u j) - K i j - (K i i + K j j) / (2 * (n - 1))

theorem pairSum_baseKernel {n : ℕ} (hn : 1 < n) (B : ℝ) (u : Fin n → ℝ)
    (K : Fin n → Fin n → ℝ) (hK : ∀ i j, K i j = K j i) (hu : ∑ i, u i = (n : ℝ) - 1) :
    pairSum (baseKernel B u K) = B - (∑ i, ∑ j, K i j) / 2 := by
  classical
  have hsym (i j : Fin n) : baseKernel B u K i j = baseKernel B u K j i := by
    simp only [baseKernel, hK i j, add_comm (u i) (u j), add_comm (K i i) (K j j)]
  have he := twice_pairSum (baseKernel B u K) hsym
  have hnR : (1 : ℝ) < n := by exact_mod_cast hn
  have hne : (n : ℝ) ≠ 0 := by linarith
  have hne' : (n : ℝ) - 1 ≠ 0 := by linarith
  calc
    pairSum (baseKernel B u K) = (2 * pairSum (baseKernel B u K)) / 2 := by ring
    _ = _ := by
      rw [he]
      simp only [baseKernel, mul_add, add_mul, div_eq_mul_inv, Finset.sum_add_distrib,
        Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.sum_mul,
        Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, hu]
      field_simp [hne, hne']
      ring

theorem pairSum_symmetrized_rows {n : ℕ} (r : Fin n → Fin n → ℝ)
    (hdiag : ∀ i, r i i = 0) (hrow : ∀ i, ∑ j, r i j = 0) :
    pairSum (fun i j => r i j + r j i) = 0 := by
  have he := twice_pairSum (fun i j => r i j + r j i) (fun i j => add_comm _ _)
  have hs : (∑ i, ∑ j, r j i) = 0 := by
    rw [Finset.sum_comm]
    simp [hrow]
  simp only [Finset.sum_add_distrib, hdiag, hrow, hs, Finset.sum_const_zero, add_zero,
    sub_zero] at he
  linarith

theorem mass_row_zero (u : Fin 10 → ℝ) (hu : ∑ i, u i = 9) (i : Fin 10) :
    (∑ j, if j = i then 0 else u j - 1 + u i / 9) = 0 := by
  classical
  have he (j : Fin 10) :
      (if j = i then 0 else u j - 1 + u i / 9) =
        (u j - 1 + u i / 9) - (if j = i then u i - 1 + u i / 9 else 0) := by
    split_ifs with h
    · subst j
      ring
    · simp
  simp only [he, Finset.sum_sub_distrib, Finset.sum_add_distrib, Finset.sum_ite_eq',
    Finset.mem_univ, if_true, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul, hu]
  ring

theorem population_row_zero (outer : Fin 10 → ℝ)
    (houter : ∑ i, outer i = 8) (i : Fin 10) :
    (∑ j, if j = i then 0 else (8 - outer i) / 9 - outer j) = 0 := by
  classical
  have he (j : Fin 10) :
      (if j = i then 0 else (8 - outer i) / 9 - outer j) =
        ((8 - outer i) / 9 - outer j) -
          (if j = i then (8 - outer i) / 9 - outer i else 0) := by
    split_ifs with h
    · subst j
      ring
    · simp
  simp only [he, Finset.sum_sub_distrib, Finset.sum_ite_eq', Finset.mem_univ, if_true,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, houter]
  ring

def cyclicEndpointSum (k : ZMod 8) (f : ZMod 8 → ℝ) : ℝ :=
  ∑ i, (f i + f (i + k))

theorem cyclicEndpointSum_eq (k : ZMod 8) (f : ZMod 8 → ℝ) :
    cyclicEndpointSum k f = 2 * ∑ i, f i := by
  have hp : (∑ i : ZMod 8, f (i + k)) = ∑ i, f i :=
    Equiv.sum_comp (Equiv.addRight k) f
  simp only [cyclicEndpointSum, Finset.sum_add_distrib, hp, two_mul]

/-- The fourth cyclic class has half the directed edge count. -/
theorem cyclic_row_zero (k : ZMod 8) (f : ZMod 8 → ℝ) :
    cyclicEndpointSum k f - 2 * (cyclicEndpointSum 4 f / 2) = 0 := by
  rw [cyclicEndpointSum_eq, cyclicEndpointSum_eq]
  ring

def gram {n : ℕ} {κ : Type*} [Fintype κ]
    (f : κ → Fin n → ℝ) (i j : Fin n) : ℝ :=
  ∑ k, f k i * f k j

theorem gram_symm {n : ℕ} {κ : Type*} [Fintype κ]
    (f : κ → Fin n → ℝ) (i j : Fin n) :
    gram f i j = gram f j i := by
  simp [gram, mul_comm]

theorem sum_gram {n : ℕ} {κ : Type*} [Fintype κ] (f : κ → Fin n → ℝ) :
    (∑ i, ∑ j, gram f i j) = ∑ k, (∑ i, f k i) ^ 2 := by
  classical
  unfold gram
  calc
    (∑ i, ∑ j, ∑ k, f k i * f k j) = ∑ i, ∑ k, ∑ j, f k i * f k j := by
      exact Finset.sum_congr rfl (fun i _ => Finset.sum_comm)
    _ = ∑ k, ∑ i, ∑ j, f k i * f k j := Finset.sum_comm
    _ = ∑ k, (∑ i, f k i) ^ 2 := by
      refine Finset.sum_congr rfl fun k _ => ?_
      simp only [← Finset.mul_sum, ← Finset.sum_mul, pow_two]

theorem sum_gram_nonneg {n : ℕ} {κ : Type*} [Fintype κ] (f : κ → Fin n → ℝ) :
    0 ≤ ∑ i, ∑ j, gram f i j := by
  rw [sum_gram]
  exact Finset.sum_nonneg fun _ _ => sq_nonneg _

theorem moments_zero_of_sum_gram_zero {n : ℕ} {κ : Type*} [Fintype κ]
    (f : κ → Fin n → ℝ) (hf : ∑ i, ∑ j, gram f i j = 0) (k : κ) :
    ∑ i, f k i = 0 := by
  classical
  rw [sum_gram] at hf
  have hk : (∑ i, f k i) ^ 2 ≤ ∑ l, (∑ i, f l i) ^ 2 :=
    Finset.single_le_sum (fun l _ => sq_nonneg (∑ i, f l i)) (Finset.mem_univ k)
  rw [hf] at hk
  nlinarith [sq_nonneg (∑ i, f k i)]

def assembledKernel {n : ℕ} {κ : Type*} [Fintype κ] (B : ℝ)
    (u : Fin n → ℝ) (f : κ → Fin n → ℝ) (W V : Fin n → Fin n → ℝ)
    (i j : Fin n) : ℝ :=
  baseKernel B u (gram f) i j + W i j - V i j

/-- The pair majorant bounds the logarithmic product. Taking equality kills
    every Gram moment. The finite certificate must supply the pairwise inequalities. -/
theorem gram_majorant {n : ℕ} {κ : Type*} [Fintype κ] {z : Configuration n}
    (hn : 1 < n) (hz : Function.Injective z) (hq : mass z = (n : ℝ) - 1)
    (B : ℝ) (f : κ → Fin n → ℝ) (W V : Fin n → Fin n → ℝ)
    (hW : pairSum W = 0) (hV : 0 ≤ pairSum V)
    (hmajor : ∀ i j, i < j →
      Real.log (distanceSq z i j) ≤ assembledKernel B (fun k => Complex.normSq (z k)) f W V i j) :
    Real.log (distanceProduct z) ≤ B ∧
      (Real.log (distanceProduct z) = B → ∀ k, ∑ i, f k i = 0) := by
  have hsum := pairSum_mono hmajor
  rw [← log_distanceProduct hz] at hsum
  change Real.log (distanceProduct z) ≤
    pairSum (fun i j => baseKernel B (fun k => Complex.normSq (z k)) (gram f) i j +
      W i j - V i j) at hsum
  rw [pairSum_sub, pairSum_add, pairSum_baseKernel hn B _ _ (gram_symm f) hq, hW, add_zero] at hsum
  have hgn := sum_gram_nonneg f
  refine ⟨by linarith, fun he k => moments_zero_of_sum_gram_zero f ?_ k⟩
  rw [he] at hsum
  linarith

/-- A certificate for a fixed configuration, with the cover inequality
    and equality reconstruction exposed as proof fields. -/
structure PairMajorant {n : ℕ} (z w : Configuration n) where
  featureCount : ℕ
  feature : Fin featureCount → Fin n → ℝ
  ward : Fin n → Fin n → ℝ
  variation : Fin n → Fin n → ℝ
  ward_zero : pairSum ward = 0
  variation_nonneg : 0 ≤ pairSum variation
  majorizes : ∀ i j, i < j →
    Real.log (distanceSq z i j) ≤
      assembledKernel (Real.log (distanceProduct w)) (fun k => Complex.normSq (z k))
        feature ward variation i j
  rigid : (∀ k, ∑ i, feature k i = 0) → SameOrbit z w

theorem log_certificate_of_pair_majorants {n : ℕ} (hn : 1 < n) (w : Configuration n)
    (h : ∀ z, IsMaximizer ((n : ℝ) - 1) z → Function.Injective z → Nonempty (PairMajorant z w)) :
    LogCertificate ((n : ℝ) - 1) w := by
  intro z hz hi
  obtain ⟨cert⟩ := h z hz hi
  obtain ⟨hb, heq⟩ := gram_majorant hn hi hz.1 _ cert.feature cert.ward cert.variation
    cert.ward_zero cert.variation_nonneg cert.majorizes
  exact ⟨hb, fun he => cert.rigid (heq he)⟩

end Mordell
