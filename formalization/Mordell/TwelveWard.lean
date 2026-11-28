import Mordell.HarmonicKernel
import Mordell.Variational

/-! Exact cancellation of the four Ward blocks in the twelve-point certificate. -/

open scoped BigOperators Matrix

noncomputable section

namespace Mordell

/-- The five unprojected radial powers used by all four Ward blocks. -/
def twelveWardFeature (z : Configuration 12) (k : ℕ) (i : Fin 12) (a : Fin 5) : ℂ :=
  ((Complex.normSq (z i) - 1 : ℝ) : ℂ) ^ a.val * z i ^ k

private theorem twelveWard_linear_moment (z : Configuration 12) (c : ℝ) (k : ℕ) :
    (c : ℂ) * (∑ i, twelveWardFeature z k i 0) + ∑ i, twelveWardFeature z k i 1 =
      (∑ i, (Complex.normSq (z i) : ℂ) * z i ^ k) + ((c - 1 : ℝ) : ℂ) * powerMoment z k := by
  simp only [twelveWardFeature, Fin.val_zero, Fin.val_one, pow_zero, pow_one, one_mul,
    powerMoment, Finset.mul_sum, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun _ _ => by push_cast; ring

theorem twelveWard_moments {z : Configuration 12} (hz : Function.Injective z)
    (hs : Stationary (6 : ℂ) z) (hmass : mass z = 11) :
    (∑ i, twelveWardFeature z 1 i 0) = 0 ∧ (∑ i, twelveWardFeature z 1 i 1) = 0 ∧
    (1 / 12 : ℂ) * (∑ i, twelveWardFeature z 0 i 0) + (∑ i, twelveWardFeature z 0 i 1) = 0 ∧
    (-3 / 4 : ℂ) * (∑ i, twelveWardFeature z 2 i 0) + (∑ i, twelveWardFeature z 2 i 1) = 0 ∧
    (-2 / 3 : ℂ) * (∑ i, twelveWardFeature z 3 i 0) + (∑ i, twelveWardFeature z 3 i 1) = 0 := by
  have hp0 : powerMoment z 0 = 12 := by simp [powerMoment]
  have hp1 : powerMoment z 1 = 0 := by
    simpa [powerMoment] using centered_of_stationary (by norm_num : (6 : ℂ) ≠ 0) hs
  have w1 := ward_identity hz hs 1
  have w2 := ward_identity hz hs 2
  have w3 := ward_identity hz hs 3
  norm_num [Finset.sum_range_succ, hp0, hp1] at w1 w2 w3
  have hq : (∑ i, (Complex.normSq (z i) : ℂ) * z i ^ 0) = 11 := by
    have h := congrArg (fun r : ℝ => (r : ℂ)) hmass
    simpa [mass, Complex.ofReal_sum] using h
  refine ⟨by simpa [twelveWardFeature, powerMoment] using hp1, ?_, ?_, ?_, ?_⟩
  · have h := twelveWard_linear_moment z 0 1
    norm_num [hp1] at h
    exact h.trans w1
  · have h := twelveWard_linear_moment z (1 / 12) 0
    simp only [pow_zero, mul_one] at hq
    norm_num [hp0, hq] at h
    exact h
  · have h := twelveWard_linear_moment z (-3 / 4) 2
    push_cast at h
    linear_combination h + w2 / 12
  · have h := twelveWard_linear_moment z (-2 / 3) 3
    push_cast at h
    linear_combination h + w3 / 12

/-- The rank-two Ward matrix uses the linear moment with coefficients `c, 1`. -/
def twelveWardRankTwo (c : ℝ) (a : Fin 5 → ℝ) : Matrix (Fin 5) (Fin 5) ℝ :=
  fun r s => (![c, 1, 0, 0, 0] : Fin 5 → ℝ) r * a s + a r * (![c, 1, 0, 0, 0] : Fin 5 → ℝ) s

private theorem matrixForm_twelveWardRankTwo (c : ℝ) (a x : Fin 5 → ℝ)
    (hx : c * x 0 + x 1 = 0) : matrixForm (twelveWardRankTwo c a) x x = 0 := by
  have h : (∑ r : Fin 5, (![c, 1, 0, 0, 0] : Fin 5 → ℝ) r * x r) = 0 := by
    simpa [Fin.sum_univ_succ] using hx
  unfold matrixForm twelveWardRankTwo
  calc
    _ = (∑ r, (![c, 1, 0, 0, 0] : Fin 5 → ℝ) r * x r) * (∑ r, a r * x r) +
      (∑ r, a r * x r) * (∑ r, (![c, 1, 0, 0, 0] : Fin 5 → ℝ) r * x r) := by
      simp only [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => by ring
    _ = 0 := by simp [h]

private theorem twelveWard_rank_two_sum (z : Configuration 12) (k : ℕ) (c : ℝ) (a : Fin 5 → ℝ)
    (h : (c : ℂ) * (∑ i, twelveWardFeature z k i 0) + ∑ i, twelveWardFeature z k i 1 = 0) :
    (∑ i, ∑ j, matrixGram (twelveWardRankTwo c a) (twelveWardFeature z k) i j) = 0 := by
  have hre := congrArg Complex.re h
  have him := congrArg Complex.im h
  simp only [Complex.add_re, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, zero_mul,
    sub_zero, Complex.zero_re] at hre
  simp only [Complex.add_im, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul,
    add_zero, Complex.zero_im] at him
  rw [sum_matrixGram, matrixForm_twelveWardRankTwo _ _ _ hre, matrixForm_twelveWardRankTwo _ _ _ him, add_zero]

private theorem twelveWard_first_sum (z : Configuration 12) (A : Matrix (Fin 5) (Fin 5) ℝ)
    (hA : ∀ r s, 2 ≤ r.val → 2 ≤ s.val → A r s = 0)
    (h0 : (∑ i, twelveWardFeature z 1 i 0) = 0) (h1 : (∑ i, twelveWardFeature z 1 i 1) = 0) :
    (∑ i, ∑ j, matrixGram A (twelveWardFeature z 1) i j) = 0 := by
  have hform (x : Fin 5 → ℝ) (hx0 : x 0 = 0) (hx1 : x 1 = 0) : matrixForm A x x = 0 := by
    apply Finset.sum_eq_zero
    intro r _
    apply Finset.sum_eq_zero
    intro s _
    by_cases hr : 2 ≤ r.val
    · by_cases hs : 2 ≤ s.val
      · simp [hA r s hr hs]
      · have h : s = 0 ∨ s = 1 := by
          fin_cases s <;> simp_all
        rcases h with (rfl | rfl) <;> simp [hx0, hx1]
    · have h : r = 0 ∨ r = 1 := by
        fin_cases r <;> simp_all
      rcases h with (rfl | rfl) <;> simp [hx0, hx1]
  rw [sum_matrixGram, hform _ (by rw [h0, Complex.zero_re]) (by rw [h1, Complex.zero_re]),
    hform _ (by rw [h0, Complex.zero_im]) (by rw [h1, Complex.zero_im]), add_zero]

/-- Each archived Ward block has zero ordered sum for a mass-eleven maximizer. -/
theorem twelveWard_block_sums {z : Configuration 12} (hz : IsMaximizer 11 z)
    (A : Matrix (Fin 5) (Fin 5) ℝ) (hA : ∀ r s, 2 ≤ r.val → 2 ≤ s.val → A r s = 0)
    (a₀ a₂ a₃ : Fin 5 → ℝ) :
    (∑ i, ∑ j, matrixGram A (twelveWardFeature z 1) i j) = 0 ∧
    (∑ i, ∑ j, matrixGram (twelveWardRankTwo (1 / 12) a₀) (twelveWardFeature z 0) i j) = 0 ∧
    (∑ i, ∑ j, matrixGram (twelveWardRankTwo (-3 / 4) a₂) (twelveWardFeature z 2) i j) = 0 ∧
    (∑ i, ∑ j, matrixGram (twelveWardRankTwo (-2 / 3) a₃) (twelveWardFeature z 3) i j) = 0 := by
  have hi := injective_of_maximizer (by decide : 1 < 12) (by norm_num : (0 : ℝ) < 11) hz
  have hm : IsMaximizer ((12 : ℝ) - 1) z := by norm_num; exact hz
  have hs : Stationary (6 : ℂ) z := by
    convert stationary_at_mass_n_sub_one (by decide : 1 < 12) hm using 1
    norm_num
  obtain ⟨h0, h1, hRadial, hSecond, hThird⟩ := twelveWard_moments hi hs hz.1
  exact ⟨twelveWard_first_sum z A hA h0 h1,
    twelveWard_rank_two_sum z 0 _ a₀ (by simpa using hRadial),
    twelveWard_rank_two_sum z 2 _ a₂ (by simpa using hSecond),
    twelveWard_rank_two_sum z 3 _ a₃ (by simpa using hThird)⟩

/-- The ordered Ward kernel before distributing the diagonal among edges. -/
def twelveWardKernel (z : Configuration 12) (A : Matrix (Fin 5) (Fin 5) ℝ)
    (a₀ a₂ a₃ : Fin 5 → ℝ) (i j : Fin 12) : ℝ :=
  matrixGram A (twelveWardFeature z 1) i j +
    matrixGram (twelveWardRankTwo (1 / 12) a₀) (twelveWardFeature z 0) i j +
    matrixGram (twelveWardRankTwo (-3 / 4) a₂) (twelveWardFeature z 2) i j +
    matrixGram (twelveWardRankTwo (-2 / 3) a₃) (twelveWardFeature z 3) i j

private theorem twelveWardRankTwo_isHermitian (c : ℝ) (a : Fin 5 → ℝ) :
    (twelveWardRankTwo c a).IsHermitian := by
  ext r s
  simp only [Matrix.conjTranspose_apply, star_trivial, twelveWardRankTwo]
  ring

/-- The diagonal-compensated archived Ward correction sums to zero on unordered edges. -/
theorem twelveWard_pairSum_zero {z : Configuration 12} (hz : IsMaximizer 11 z)
    (A : Matrix (Fin 5) (Fin 5) ℝ) (hSym : A.IsHermitian)
    (hA : ∀ r s, 2 ≤ r.val → 2 ≤ s.val → A r s = 0) (a₀ a₂ a₃ : Fin 5 → ℝ) :
    pairSum (fun i j => -twelveWardKernel z A a₀ a₂ a₃ i j -
      (twelveWardKernel z A a₀ a₂ a₃ i i + twelveWardKernel z A a₀ a₂ a₃ j j) / 22) = 0 := by
  let K := twelveWardKernel z A a₀ a₂ a₃
  have hs (i j : Fin 12) : K i j = K j i := by
    unfold K twelveWardKernel
    rw [matrixGram_symm hSym _ i j, matrixGram_symm (twelveWardRankTwo_isHermitian _ _) _ i j,
      matrixGram_symm (twelveWardRankTwo_isHermitian _ _) _ i j,
      matrixGram_symm (twelveWardRankTwo_isHermitian _ _) _ i j]
  have hsum : (∑ i, ∑ j, K i j) = 0 := by
    obtain ⟨h1, h0, h2, h3⟩ := twelveWard_block_sums hz A hA a₀ a₂ a₃
    simp only [K, twelveWardKernel, Finset.sum_add_distrib, h1, h0, h2, h3, add_zero]
  let u (i : Fin 12) := Complex.normSq (z i)
  have hu : (∑ i, u i) = (12 : ℝ) - 1 := by norm_num; exact hz.1
  have hc := pairSum_diagonal_correction (by decide : 1 < 12) u hu K hs
  norm_num [hsum] at hc
  exact hc

end Mordell
