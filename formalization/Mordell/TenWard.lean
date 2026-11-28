import Mordell.HarmonicKernel
import Mordell.Variational

/-! The Ward, mass, population and radial-force corrections of the ten-point certificate. -/

open scoped BigOperators Matrix ComplexConjugate

noncomputable section

namespace Mordell

def tenWardHarmonic (row : Fin 5) : ℕ := ![0, 1, 1, 2, 3] row

def tenWardMultiplier (row : Fin 5) (u : ℝ) : ℝ :=
  ![u - 9 / 10, 1, u - 1, u - 17 / 10, u - 8 / 5] row

private theorem ten_stationary {z : Configuration 10} (hz : IsMaximizer 9 z) :
    Stationary (5 : ℂ) z := by
  have hm : IsMaximizer ((10 : ℝ) - 1) z := by norm_num; exact hz
  convert stationary_at_mass_n_sub_one (by decide : 1 < 10) hm using 1
  norm_num

private theorem shifted_powerMoment (z : Configuration 10) (c : ℝ) (k : ℕ) :
    (∑ i, ((Complex.normSq (z i) - c : ℝ) : ℂ) * z i ^ k) =
      (∑ i, (Complex.normSq (z i) : ℂ) * z i ^ k) - (c : ℂ) * powerMoment z k := by
  simp only [Complex.ofReal_sub, sub_mul, Finset.sum_sub_distrib, ← Finset.mul_sum, powerMoment]

/-- The five linear moments used by the archived Ward columns vanish at every maximizer. -/
theorem tenWard_moments {z : Configuration 10} (hz : IsMaximizer 9 z) (row : Fin 5) :
    (∑ i, (tenWardMultiplier row (Complex.normSq (z i)) : ℂ) * z i ^ tenWardHarmonic row) = 0 := by
  have hi := injective_of_maximizer (by decide : 1 < 10) (by norm_num : (0 : ℝ) < 9) hz
  have hs := ten_stationary hz
  have hp0 : powerMoment z 0 = 10 := by simp [powerMoment]
  have hp1 : powerMoment z 1 = 0 := by
    simpa [powerMoment] using centered_of_stationary (by norm_num : (5 : ℂ) ≠ 0) hs
  have w1 := ward_identity hi hs 1
  have w2 := ward_identity hi hs 2
  have w3 := ward_identity hi hs 3
  norm_num [Finset.sum_range_succ, hp0, hp1] at w1 w2 w3
  have hmass : (∑ i, (Complex.normSq (z i) : ℂ)) = 9 := by
    have h := congrArg (fun r : ℝ => (r : ℂ)) hz.1
    simpa [mass, Complex.ofReal_sum] using h
  fin_cases row
  · change (∑ i, ((Complex.normSq (z i) - 9 / 10 : ℝ) : ℂ) * z i ^ 0) = 0
    rw [shifted_powerMoment]
    norm_num [hp0, hmass]
  · change (∑ i, ((1 : ℝ) : ℂ) * z i ^ 1) = 0
    simpa [powerMoment] using hp1
  · change (∑ i, ((Complex.normSq (z i) - 1 : ℝ) : ℂ) * z i ^ 1) = 0
    rw [shifted_powerMoment]
    simpa [hp1] using w1
  · change (∑ i, ((Complex.normSq (z i) - 17 / 10 : ℝ) : ℂ) * z i ^ 2) = 0
    rw [shifted_powerMoment]
    push_cast
    linear_combination w2 / 10
  · change (∑ i, ((Complex.normSq (z i) - 8 / 5 : ℝ) : ℂ) * z i ^ 3) = 0
    rw [shifted_powerMoment]
    push_cast
    linear_combination w3 / 10

/-- The archived signed Ward column, including its distributed diagonal. -/
def tenWardCorrection (row : Fin 5) (g : ℝ → ℝ) (u v t : ℝ) : ℝ :=
  let w := tenWardMultiplier row
  let k := tenWardHarmonic row;
  -(harmonicPolynomial (u * v) t k * (w u * g v + w v * g u)) -
    (u ^ k * w u * g u + v ^ k * w v * g v) / 9

theorem tenWardCorrection_pairSum_zero {z : Configuration 10} (hz : IsMaximizer 9 z)
    (row : Fin 5) (g : ℝ → ℝ) :
    pairSum (fun i j => tenWardCorrection row g (Complex.normSq (z i))
      (Complex.normSq (z j)) (realDot (z i) (z j))) = 0 := by
  let u (i : Fin 10) := Complex.normSq (z i)
  let w := tenWardMultiplier row
  let k := tenWardHarmonic row
  let r (i : Fin 10) : Fin 2 → ℝ := ![w (u i), g (u i)]
  let Q : Matrix (Fin 2) (Fin 2) ℝ := !![0, 1; 1, 0]
  let f (i : Fin 10) (a : Fin 2) : ℂ := (r i a : ℂ) * z i ^ k
  let K := matrixGram Q f
  have hQ : Q.IsHermitian := by
    ext a b
    fin_cases a <;> fin_cases b <;> norm_num [Q]
  have hm : (∑ i, f i 0) = 0 := tenWard_moments hz row
  have hsum : (∑ i, ∑ j, K i j) = 0 := by
    have hr : (∑ i, (f i 0).re) = 0 := by
      have h : (∑ i, f i 0).re = ∑ i, (f i 0).re :=
        map_sum Complex.reAddGroupHom (fun i => f i 0) Finset.univ
      simp only [hm, Complex.zero_re] at h
      exact h.symm
    have hi : (∑ i, (f i 0).im) = 0 := by
      have h : (∑ i, f i 0).im = ∑ i, (f i 0).im :=
        map_sum Complex.imAddGroupHom (fun i => f i 0) Finset.univ
      simp only [hm, Complex.zero_im] at h
      exact h.symm
    unfold K
    rw [sum_matrixGram]
    simp [matrixForm, Fin.sum_univ_two, Q, hr, hi]
  have hK (i j : Fin 10) : K i j = harmonicPolynomial (u i * u j) (realDot (z i) (z j)) k *
      (w (u i) * g (u j) + w (u j) * g (u i)) := by
    unfold K f
    rw [matrixGram_harmonic]
    simp [matrixForm, Fin.sum_univ_two, Q, r, u]
    simp [mul_comm (g (Complex.normSq (z i)))]
  have hdiag (i : Fin 10) : K i i = 2 * u i ^ k * w (u i) * g (u i) := by
    rw [hK]
    change harmonicPolynomial (Complex.normSq (z i) * Complex.normSq (z i)) _ _ * _ = _
    rw [harmonicPolynomial_pair, realDot, Complex.mul_conj, Complex.ofReal_re, map_pow]
    simp only [u]
    ring
  have hu : (∑ i, u i) = (10 : ℝ) - 1 := by norm_num; exact hz.1
  have hc := pairSum_diagonal_correction (by decide : 1 < 10) u hu K (matrixGram_symm hQ f)
  have he (i j : Fin 10) : -K i j - (K i i + K j j) / 18 =
      tenWardCorrection row g (u i) (u j) (realDot (z i) (z j)) := by
    rw [hdiag, hdiag, hK]
    simp only [tenWardCorrection, k, w]
    ring
  norm_num [hsum] at hc
  simpa only [he] using hc

/-- The real radial-force row with its diagonal removed. -/
def tenRadialForceRow (u t d : ℝ) : ℝ := (u - t) / d - (5 / 9) * u

theorem tenRadialForceCorrection_pairSum_zero {z : Configuration 10} (hz : IsMaximizer 9 z)
    (a : Fin 10 → ℝ) :
    pairSum (fun i j =>
      a i * tenRadialForceRow (Complex.normSq (z i)) (realDot (z i) (z j)) (distanceSq z i j) +
      a j * tenRadialForceRow (Complex.normSq (z j)) (realDot (z j) (z i)) (distanceSq z j i)) = 0 := by
  apply pairSum_weighted_rows
  intro i
  have hs : Stationary ((5 : ℝ) : ℂ) z := by simpa using ten_stationary hz
  let C := (5 / 9 : ℝ) * Complex.normSq (z i)
  have hc (j : Fin 10) : (if i = j then 0 else C) = C - (if i = j then C else 0) := by
    split_ifs <;> simp
  have hsum : (∑ j, if i = j then 0 else C) = 5 * Complex.normSq (z i) := by
    simp only [hc, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
      nsmul_eq_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true, C]
    ring
  have he (j : Fin 10) :
      (if i = j then 0 else tenRadialForceRow (Complex.normSq (z i)) (realDot (z i) (z j))
        (distanceSq z i j)) =
      (if i = j then 0 else (Complex.normSq (z i) - realDot (z i) (z j)) / distanceSq z i j) -
      (if i = j then 0 else C) := by
    split_ifs <;> simp [tenRadialForceRow, C]
  simp only [he, Finset.sum_sub_distrib, hsum, radial_force_sum hs, sub_self]

theorem tenMassCorrection_pairSum_zero {z : Configuration 10} (hz : IsMaximizer 9 z)
    (a : Fin 10 → ℝ) :
    pairSum (fun i j => a i * (Complex.normSq (z j) - 1 + Complex.normSq (z i) / 9) +
      a j * (Complex.normSq (z i) - 1 + Complex.normSq (z j) / 9)) = 0 := by
  apply pairSum_weighted_rows
  intro i
  simpa only [eq_comm] using mass_row_zero (fun j => Complex.normSq (z j)) hz.1 i

theorem tenPopulationCorrection_pairSum_zero (outer a : Fin 10 → ℝ)
    (houter : ∑ i, outer i = 8) :
    pairSum (fun i j => a i * ((8 - outer i) / 9 - outer j) +
      a j * ((8 - outer j) / 9 - outer i)) = 0 := by
  apply pairSum_weighted_rows
  intro i
  simpa only [eq_comm] using population_row_zero outer houter i

/-- A directed step among the eight outer points in their cyclic order. -/
def tenOuterStep (e : ZMod 8 ↪ Fin 10) (s : ZMod 8) (i j : Fin 10) : ℝ :=
  ∑ l : ZMod 8, if i = e l ∧ j = e (l + s) then 1 else 0

theorem tenOuterStep_on_cycle (e : ZMod 8 ↪ Fin 10) (s x y : ZMod 8) :
    tenOuterStep e s (e x) (e y) = if y = x + s then 1 else 0 := by
  unfold tenOuterStep
  have h (l : ZMod 8) : (if e x = e l ∧ e y = e (l + s) then (1 : ℝ) else 0) =
      if l = x then (if y = x + s then 1 else 0) else 0 := by
    by_cases hl : l = x
    · subst l
      simp [e.injective.eq_iff]
    · simp [e.injective.eq_iff, hl, Ne.symm hl]
  simp only [h, Finset.sum_ite_eq', Finset.mem_univ, if_true]

private theorem tenOuterStep_flip (e : ZMod 8 ↪ Fin 10) (s : ZMod 8) (i j : Fin 10) :
    tenOuterStep e s j i = tenOuterStep e (-s) i j := by
  unfold tenOuterStep
  rw [← Equiv.sum_comp (Equiv.addRight s)
    (fun l => if i = e l ∧ j = e (l + -s) then (1 : ℝ) else 0)]
  apply Finset.sum_congr rfl
  intro l _
  change (if j = e l ∧ i = e (l + s) then (1 : ℝ) else 0) =
    if i = e (l + s) ∧ j = e (l + s + -s) then 1 else 0
  simp [and_comm]

private theorem tenOuterStep_sum (e : ZMod 8 ↪ Fin 10) (s : ZMod 8) (i : Fin 10) :
    (∑ j, tenOuterStep e s i j) = ∑ l, if i = e l then 1 else 0 := by
  unfold tenOuterStep
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro l _
  by_cases h : i = e l <;>
    simp [h, Finset.sum_ite_eq']

private theorem tenOuterStep_diagonal (e : ZMod 8 ↪ Fin 10) {s : ZMod 8}
    (hs : s ≠ 0) (i : Fin 10) : tenOuterStep e s i i = 0 := by
  apply Finset.sum_eq_zero
  intro l _
  apply if_neg
  rintro ⟨h, h'⟩
  have he := e.injective (h.symm.trans h')
  exact hs (by simpa using congrArg (· - l) he.symm)

/-- One of the two-way cyclic classes minus twice the opposite class has
    zero weighted edge sum. This applies to steps one, two and three. -/
theorem tenOuterAngleCorrection_pairSum_zero (e : ZMod 8 ↪ Fin 10) (s : ZMod 8)
    (hs : s ≠ 0) (a : Fin 10 → ℝ) :
    pairSum (fun i j => (a i + a j) *
      (tenOuterStep e s i j + tenOuterStep e (-s) i j - 2 * tenOuterStep e 4 i j)) = 0 := by
  let R (i j : Fin 10) := tenOuterStep e s i j + tenOuterStep e (-s) i j - 2 * tenOuterStep e 4 i j
  have hSym (i j : Fin 10) : R j i = R i j := by
    simp only [R, tenOuterStep_flip e _ i j, neg_neg, show (-4 : ZMod 8) = 4 by decide]
    rw [add_comm]
  have hDiag (i : Fin 10) : R i i = 0 := by
    simp only [R, tenOuterStep_diagonal e hs, tenOuterStep_diagonal e (neg_ne_zero.mpr hs),
      tenOuterStep_diagonal e (by decide : (4 : ZMod 8) ≠ 0), add_zero, mul_zero, sub_zero]
  have hRow (i : Fin 10) : (∑ j, if i = j then 0 else R i j) = 0 := by
    have h (j : Fin 10) : (if i = j then 0 else R i j) = R i j := by
      by_cases h : i = j
      · subst j
        simp [hDiag]
      · simp [h]
    simp only [h, R, Finset.sum_sub_distrib, Finset.sum_add_distrib, ← Finset.mul_sum, tenOuterStep_sum]
    ring
  have hEntry (i j : Fin 10) : (a i + a j) * R i j = a i * R i j + a j * R j i := by
    rw [hSym i j, add_mul]
  simpa only [← hEntry] using pairSum_weighted_rows a R hRow

end Mordell
