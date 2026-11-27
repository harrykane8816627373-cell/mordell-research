import Mordell.Basic

/-! The concrete supporting lines for the twelve-point radial population cut. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def twelveRadialCubic (u : ℝ) : ℝ :=
  u ^ 3 - 4 * u ^ 2 + 13118335840468817 / 4503599627370496

private def gap (s intercept slope x : ℝ) : ℝ :=
  s * twelveRadialCubic x - intercept - slope * x

private def gapDerivative (s slope x : ℝ) : ℝ := s * (3 * x ^ 2 - 8 * x) - slope

private theorem cubic_gap_nonneg {s p m a b x : ℝ}
    (ha : a ≤ x) (hb : x ≤ b) (hab : a < b)
    (h0 : 0 ≤ gap s p m a)
    (h1 : 0 ≤ gap s p m a + (b - a) * gapDerivative s m a / 3)
    (h2 : 0 ≤ gap s p m b - (b - a) * gapDerivative s m b / 3)
    (h3 : 0 ≤ gap s p m b) :
    0 ≤ gap s p m x := by
  have hxa := sub_nonneg.mpr ha
  have hbx := sub_nonneg.mpr hb
  have he :
      (b - a) ^ 3 * gap s p m x =
        gap s p m a * (b - x) ^ 3 +
        3 * (gap s p m a + (b - a) * gapDerivative s m a / 3) * (b - x) ^ 2 * (x - a) +
        3 * (gap s p m b - (b - a) * gapDerivative s m b / 3) * (b - x) * (x - a) ^ 2 +
        gap s p m b * (x - a) ^ 3 := by
    unfold gap gapDerivative twelveRadialCubic
    ring
  have hnonneg : 0 ≤ (b - a) ^ 3 * gap s p m x := by
    rw [he]
    positivity
  exact nonneg_of_mul_nonneg_right hnonneg (pow_pos (sub_pos.mpr hab) _)

theorem twelveRadialCubic_inner {u : ℝ} (hu : u ∈ Set.Icc (3 / 25) (6 / 25)) :
    347 / 100 - (26 / 5) * u ≤ twelveRadialCubic u ∧
      -(199 / 50) + (531 / 100) * u ≤ -twelveRadialCubic u := by
  obtain ⟨hlo, hhi⟩ := hu
  have hp := cubic_gap_nonneg (s := 1) (p := 347 / 100) (m := -(26 / 5)) hlo hhi
    (by norm_num) (by norm_num [gap, twelveRadialCubic])
    (by norm_num [gap, gapDerivative, twelveRadialCubic])
    (by norm_num [gap, gapDerivative, twelveRadialCubic])
    (by norm_num [gap, twelveRadialCubic])
  have hm := cubic_gap_nonneg (s := -1) (p := -(199 / 50)) (m := 531 / 100) hlo hhi
    (by norm_num) (by norm_num [gap, twelveRadialCubic])
    (by norm_num [gap, gapDerivative, twelveRadialCubic])
    (by norm_num [gap, gapDerivative, twelveRadialCubic])
    (by norm_num [gap, twelveRadialCubic])
  unfold gap at hp hm
  exact ⟨by linarith, by linarith⟩

theorem twelveRadialCubic_outer {u : ℝ} (hu : u ∈ Set.Icc (47 / 50) (71 / 50)) :
    509 / 100 - (26 / 5) * u ≤ twelveRadialCubic u ∧
      -(263 / 50) + (531 / 100) * u ≤ -twelveRadialCubic u := by
  obtain ⟨hlo, hhi⟩ := hu
  have hp := cubic_gap_nonneg (s := 1) (p := 509 / 100) (m := -(26 / 5)) hlo hhi
    (by norm_num) (by norm_num [gap, twelveRadialCubic])
    (by norm_num [gap, gapDerivative, twelveRadialCubic])
    (by norm_num [gap, gapDerivative, twelveRadialCubic])
    (by norm_num [gap, twelveRadialCubic])
  have hm : 0 ≤ gap (-1) (-(263 / 50)) (531 / 100) u := by
    by_cases hmid : u ≤ 59 / 50
    · exact cubic_gap_nonneg hlo hmid (by norm_num) (by norm_num [gap, twelveRadialCubic])
        (by norm_num [gap, gapDerivative, twelveRadialCubic])
        (by norm_num [gap, gapDerivative, twelveRadialCubic])
        (by norm_num [gap, twelveRadialCubic])
    · exact cubic_gap_nonneg (le_of_not_ge hmid) hhi (by norm_num)
        (by norm_num [gap, twelveRadialCubic])
        (by norm_num [gap, gapDerivative, twelveRadialCubic])
        (by norm_num [gap, gapDerivative, twelveRadialCubic])
        (by norm_num [gap, twelveRadialCubic])
  unfold gap at hp hm
  exact ⟨by linarith, by linarith⟩

/-- On the two necessary radial bands, the mass and the radial moment
    budget force three points in the inner band. -/
theorem twelve_inner_population_of_radial_budget (u : Fin 12 → ℝ)
    (hu : ∀ i, u i ∈ Set.Icc (3 / 25) (6 / 25) ∨ u i ∈ Set.Icc (47 / 50) (71 / 50))
    (hmass : ∑ i, u i = 11)
    (hbudget : |∑ i, twelveRadialCubic (u i)| < 41 / 100) :
    (Finset.univ.filter (fun i => u i ≤ 6 / 25)).card = 3 := by
  classical
  let k := (Finset.univ.filter (fun i => u i ≤ 6 / 25)).card
  let indicator (i : Fin 12) : ℝ := if u i ≤ 6 / 25 then 1 else 0
  have hcounts : (∑ i, indicator i) = (k : ℝ) := by
    simp [indicator, k]
  have h₁ : (388 - 162 * (k : ℝ)) / 100 ≤ ∑ i, twelveRadialCubic (u i) := by
    have hsum := Finset.sum_le_sum (s := (Finset.univ : Finset (Fin 12))) (fun i _ => show
        509 / 100 - (162 / 100) * indicator i - (26 / 5) * u i ≤ twelveRadialCubic (u i) from by
      rcases hu i with hinner | houter
      · simp only [indicator, if_pos hinner.2]
        linarith [(twelveRadialCubic_inner hinner).1]
      · have hn : ¬ u i ≤ 6 / 25 := by linarith [houter.1]
        simp only [indicator, if_neg hn]
        linarith [(twelveRadialCubic_outer houter).1])
    simp only [Finset.sum_sub_distrib, ← Finset.mul_sum, hcounts, hmass,
      Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
    linarith
  have h₂ : (128 * (k : ℝ) - 471) / 100 ≤ -(∑ i, twelveRadialCubic (u i)) := by
    have hsum := Finset.sum_le_sum (s := (Finset.univ : Finset (Fin 12))) (fun i _ => show
        -(263 / 50) + (64 / 50) * indicator i + (531 / 100) * u i ≤ -twelveRadialCubic (u i) from by
      rcases hu i with hinner | houter
      · simp only [indicator, if_pos hinner.2]
        linarith [(twelveRadialCubic_inner hinner).2]
      · have hn : ¬ u i ≤ 6 / 25 := by linarith [houter.1]
        simp only [indicator, if_neg hn]
        linarith [(twelveRadialCubic_outer houter).2])
    simp only [Finset.sum_add_distrib, Finset.sum_neg_distrib, ← Finset.mul_sum, hcounts, hmass,
      Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
    linarith
  change k = 3
  by_contra hne
  rcases (show k ≤ 2 ∨ 4 ≤ k by omega) with hs | hs
  · have hs' : (k : ℝ) ≤ 2 := by exact_mod_cast hs
    linarith [(abs_lt.mp hbudget).2]
  · have hs' : (4 : ℝ) ≤ k := by exact_mod_cast hs
    linarith [(abs_lt.mp hbudget).1]

end Mordell
