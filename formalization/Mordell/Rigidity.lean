import Mordell.Candidate

/-! Equality support and reconstruction of the antipodal layers. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

def radialPenalty (p : Parameters) (u : ℝ) : ℝ :=
  ((u - p.a) * (u - p.b) * (u - p.c)) ^ 2

theorem support_of_radialPenalty_zero {n : ℕ} (p : Parameters) (u : Fin n → ℝ)
    (hz : ∑ i, radialPenalty p (u i) = 0) (i : Fin n) :
    u i = p.a ∨ u i = p.b ∨ u i = p.c := by
  have hi : radialPenalty p (u i) ≤ ∑ j, radialPenalty p (u j) :=
    Finset.single_le_sum (fun j _ => show 0 ≤ radialPenalty p (u j) from sq_nonneg _) (Finset.mem_univ i)
  rw [hz] at hi
  unfold radialPenalty at hi
  have he : (u i - p.a) * (u i - p.b) * (u i - p.c) = 0 := by nlinarith
  rcases mul_eq_zero.mp he with hab | hc
  · rcases mul_eq_zero.mp hab with ha | hb
    · exact Or.inl (sub_eq_zero.mp ha)
    · exact Or.inr (Or.inl (sub_eq_zero.mp hb))
  · exact Or.inr (Or.inr (sub_eq_zero.mp hc))

private theorem mul_conj_add_eq {r : ℝ} {a b : ℂ}
    (ha : Complex.normSq a = r) (hb : Complex.normSq b = r) :
    a * b * conj (a + b) = (r : ℂ) * (a + b) := by
  calc
    a * b * conj (a + b) = (a * conj a) * b + (b * conj b) * a := by
      rw [map_add]
      ring
    _ = _ := by rw [Complex.mul_conj, Complex.mul_conj, ha, hb]; ring

/-- An unordered pair on a circle with a prescribed nonzero sum is unique. -/
theorem same_radius_pair_eq {r : ℝ} {a b c d : ℂ}
    (ha : Complex.normSq a = r) (hb : Complex.normSq b = r)
    (hc : Complex.normSq c = r) (hd : Complex.normSq d = r)
    (hs : a + b = c + d) (hne : a + b ≠ 0) :
    (a = c ∧ b = d) ∨ (a = d ∧ b = c) := by
  have hconj : conj (a + b) ≠ 0 := by
    intro h
    have hh := congrArg (starRingEnd ℂ) h
    exact hne (by simpa using hh)
  have hab := mul_conj_add_eq ha hb
  have hcd := mul_conj_add_eq hc hd
  rw [← hs] at hcd
  have he : a * b = c * d := mul_right_cancel₀ hconj (hab.trans hcd.symm)
  have hz : (a - c) * (a - d) = 0 := by
    calc
      (a - c) * (a - d) = a ^ 2 - a * (c + d) + c * d := by ring
      _ = 0 := by rw [← hs, ← he]; ring
  rcases mul_eq_zero.mp hz with h | h
  · have hac := sub_eq_zero.mp h
    exact Or.inl ⟨hac, add_left_cancel (by simpa only [hac] using hs)⟩
  · have had := sub_eq_zero.mp h
    exact Or.inr ⟨had, add_left_cancel (by simpa only [had, add_comm c d] using hs)⟩

/-- Four points on a circle whose first moment vanishes split into two
    antipodal pairs, with the three alternatives specifying the pairing. -/
theorem four_points_antipodal {r : ℝ} {a b c d : ℂ}
    (ha : Complex.normSq a = r) (hb : Complex.normSq b = r)
    (hc : Complex.normSq c = r) (hd : Complex.normSq d = r)
    (h : a + b + c + d = 0) :
    (b = -a ∧ d = -c) ∨ (c = -a ∧ d = -b) ∨ (d = -a ∧ c = -b) := by
  have hs : a + b = -c + -d := by linear_combination h
  by_cases hzero : a + b = 0
  · exact Or.inl ⟨by linear_combination hzero, by linear_combination h - hzero⟩
  · have hc' : Complex.normSq (-c) = r := by simpa using hc
    have hd' : Complex.normSq (-d) = r := by simpa using hd
    rcases same_radius_pair_eq ha hb hc' hd' hs hzero with ⟨hac, hbd⟩ | ⟨had, hbc⟩
    · exact Or.inr (Or.inl ⟨by linear_combination hac, by linear_combination hbd⟩)
    · exact Or.inr (Or.inr ⟨by linear_combination had, by linear_combination hbc⟩)

private theorem eq_or_eq_neg_of_sq {a b : ℂ} (h : a ^ 2 = b ^ 2) :
    a = b ∨ a = -b := by
  have hz : (a - b) * (a + b) = 0 := by linear_combination h
  rcases mul_eq_zero.mp hz with hz | hz
  · exact Or.inl (sub_eq_zero.mp hz)
  · exact Or.inr (by linear_combination hz)

/-- The second moment of a four-point antipodal layer fixes both pairs
    up to sign and interchange once the common phase is fixed. -/
theorem antipodal_layer_eq {r : ℝ} {a b c d : ℂ}
    (ha : Complex.normSq a = r) (hb : Complex.normSq b = r)
    (hc : Complex.normSq c = r) (hd : Complex.normSq d = r)
    (hs : a ^ 2 + b ^ 2 = c ^ 2 + d ^ 2) (hne : a ^ 2 + b ^ 2 ≠ 0) :
    ((a = c ∨ a = -c) ∧ (b = d ∨ b = -d)) ∨
      ((a = d ∨ a = -d) ∧ (b = c ∨ b = -c)) := by
  have hsqr (z : ℂ) (hz : Complex.normSq z = r) : Complex.normSq (z ^ 2) = r ^ 2 := by
    simp [pow_two, Complex.normSq_mul, hz]
  rcases same_radius_pair_eq (hsqr a ha) (hsqr b hb) (hsqr c hc) (hsqr d hd) hs hne with
    ⟨h₁, h₂⟩ | ⟨h₁, h₂⟩
  · exact Or.inl ⟨eq_or_eq_neg_of_sq h₁, eq_or_eq_neg_of_sq h₂⟩
  · exact Or.inr ⟨eq_or_eq_neg_of_sq h₁, eq_or_eq_neg_of_sq h₂⟩

theorem antipodal_cross_product (z w : ℂ) :
    Complex.normSq (z - w) * Complex.normSq (z - -w) *
      Complex.normSq (-z - w) * Complex.normSq (-z - -w) =
        Complex.normSq (z ^ 2 - w ^ 2) ^ 2 := by
  have h₁ : -z - w = -(z + w) := by ring
  have h₂ : -z - -w = -(z - w) := by ring
  have h₃ : z ^ 2 - w ^ 2 = (z - w) * (z + w) := by ring
  rw [h₁, h₂, h₃, Complex.normSq_neg, Complex.normSq_neg, Complex.normSq_mul, sub_neg_eq_add]
  ring

/-- Three distinct radial nodes are determined separately by the first
    three radial moments. Coefficients can also be complex layer moments. -/
theorem three_node_vandermonde {a b c A B C : ℂ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (h0 : A + B + C = 0) (h1 : a * A + b * B + c * C = 0)
    (h2 : a ^ 2 * A + b ^ 2 * B + c ^ 2 * C = 0) :
    A = 0 ∧ B = 0 ∧ C = 0 := by
  have ha : (a - b) * (a - c) * A = 0 := by
    linear_combination h2 - (b + c) * h1 + b * c * h0
  have hb : (b - a) * (b - c) * B = 0 := by
    linear_combination h2 - (a + c) * h1 + a * c * h0
  have ha0 := (mul_eq_zero.mp ha).resolve_left
    (mul_ne_zero (sub_ne_zero.mpr hab) (sub_ne_zero.mpr hac))
  have hb0 := (mul_eq_zero.mp hb).resolve_left
    (mul_ne_zero (sub_ne_zero.mpr hab.symm) (sub_ne_zero.mpr hbc))
  exact ⟨ha0, hb0, by simpa [ha0, hb0] using h0⟩

theorem three_points_cubes_eq {a b c : ℂ}
    (h1 : a + b + c = 0) (h2 : a ^ 2 + b ^ 2 + c ^ 2 = 0) :
    a ^ 3 = a * b * c ∧ b ^ 3 = a * b * c ∧ c ^ 3 = a * b * c := by
  have he : a * b + a * c + b * c = 0 := by
    linear_combination ((a + b + c) / 2) * h1 - h2 / 2
  constructor
  · linear_combination a ^ 2 * h1 - a * he
  · constructor
    · linear_combination b ^ 2 * h1 - b * he
    · linear_combination c ^ 2 * h1 - c * he

/-- The 3,3,6 population follows from the first three radial moments. -/
theorem populations_twelve_of_moments {a b c : ℝ} {A B C : ℕ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c) (hcount : A + B + C = 12)
    (hmass : (A : ℝ) * a + B * b + C * c = 3 * a + 3 * b + 6 * c)
    (hsecond : (A : ℝ) * a ^ 2 + B * b ^ 2 + C * c ^ 2 = 3 * a ^ 2 + 3 * b ^ 2 + 6 * c ^ 2) :
    A = 3 ∧ B = 3 ∧ C = 6 := by
  have hc0 : (A : ℂ) + B + C = 12 := by exact_mod_cast hcount
  have hc1 := congrArg Complex.ofReal hmass
  have hc2 := congrArg Complex.ofReal hsecond
  push_cast at hc1 hc2
  obtain ⟨hA, hB, hC⟩ := three_node_vandermonde
    (a := (a : ℂ)) (b := (b : ℂ)) (c := (c : ℂ))
    (A := (A : ℂ) - 3) (B := (B : ℂ) - 3) (C := (C : ℂ) - 6)
    (by exact_mod_cast hab) (by exact_mod_cast hac) (by exact_mod_cast hbc)
    (by linear_combination hc0) (by linear_combination hc1) (by linear_combination hc2)
  exact ⟨by exact_mod_cast sub_eq_zero.mp hA, by exact_mod_cast sub_eq_zero.mp hB,
    by exact_mod_cast sub_eq_zero.mp hC⟩

end Mordell
