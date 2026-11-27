import Mordell.Layers
import Mordell.Rigidity

/-! Ten-point rigidity from the equality moments of the sharp Gram kernel. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

structure TenLayerMoments (p : Parameters) (z : Configuration 10) : Prop where
  radial : ∑ i, radialPenalty p (Complex.normSq (z i)) = 0
  first_a : layerMoment z p.a 1 = 0
  first_b : layerMoment z p.b 1 = 0
  first_c : layerMoment z p.c 1 = 0
  second : ∃ η : ℂ,
    layerMoment z p.a 2 = η * (2 * (p.a : ℂ)) ∧
    layerMoment z p.b 2 = η * (4 * (p.x : ℂ)) ∧
    layerMoment z p.c 2 = η * (4 * (p.y : ℂ))

private def outerBandAffine (u : ℝ) : ℝ :=
  (2 * u - (89 / 100 + 1271 / 1000)) / (1271 / 1000 - 89 / 100)

/-- The six band-supported Chebyshev features, ordered inner block first. -/
def tenBandFeatures (u : ℝ) : Fin 6 → ℝ :=
  if u ≤ 159 / 1000 then
    let t := (2 * u - (21 / 250 + 159 / 1000)) / (159 / 1000 - 21 / 250)
    ![1, t, 2 * t ^ 2 - 1, 0, 0, 0]
  else
    let t := outerBandAffine u
    ![0, 0, 0, 1, t, 2 * t ^ 2 - 1]

def tenBandMoment (z : Configuration 10) (k : ℕ) (j : Fin 6) : ℂ :=
  ∑ i, (tenBandFeatures (Complex.normSq (z i)) j : ℂ) * z i ^ k

def tenSecondReference (p : Parameters) (j : Fin 6) : ℂ :=
  (tenBandFeatures p.a j : ℂ) * (2 * (p.a : ℂ)) +
    (tenBandFeatures p.b j : ℂ) * (4 * (p.x : ℂ)) +
    (tenBandFeatures p.c j : ℂ) * (4 * (p.y : ℂ))

structure TenGramMoments (p : Parameters) (z : Configuration 10) (pivot : Fin 6) : Prop where
  radial : ∑ i, radialPenalty p (Complex.normSq (z i)) = 0
  first : ∀ j, tenBandMoment z 1 j = 0
  secondProjected : ∀ j, j ≠ pivot → tenBandMoment z 2 j -
    (tenSecondReference p j / tenSecondReference p pivot) * tenBandMoment z 2 pivot = 0

private theorem tenBand_partition {p : Parameters} (hp : p.InRootBox) :
    p.a ≤ 159 / 1000 ∧ ¬ p.b ≤ 159 / 1000 ∧ ¬ p.c ≤ 159 / 1000 := by
  rcases hp with ⟨⟨_, ha₂⟩, ⟨hb₁, _⟩, ⟨hc₁, _⟩, _, _⟩
  exact ⟨by linarith, by linarith, by linarith⟩

/-- The second-harmonic pivot used by the ten-point certificate is nonzero. -/
theorem tenSecondReference_four_ne_zero {p : Parameters} (hp : p.InRootBox) :
    tenSecondReference p 4 ≠ 0 := by
  obtain ⟨hba, hbb, hbc⟩ := tenBand_partition hp
  obtain ⟨_, ⟨_, hb₂⟩, ⟨hc₁, _⟩, ⟨_, hx₂⟩, ⟨hy₁, _⟩⟩ := hp
  have hb : outerBandAffine p.b < 0 := by
    apply div_neg_of_neg_of_pos
    · linarith
    · norm_num
  have hc : 0 < outerBandAffine p.c := by
    apply div_pos
    · linarith
    · norm_num
  have hx : p.x < 0 := by linarith
  have hy : 0 < p.y := by linarith
  intro he
  have hh := congrArg Complex.re he
  simp [tenSecondReference, tenBandFeatures, hba, hbb, hbc, Complex.add_re, Complex.mul_re] at hh
  linarith only [hh, mul_pos_of_neg_of_neg hb (show 4 * p.x < 0 by linarith),
    mul_pos hc (show 0 < 4 * p.y by linarith)]

private theorem tenBand_layerMoments_eq {p : Parameters} {z : Configuration 10}
    (hp : p.InRootBox) (hr : ∑ i, radialPenalty p (Complex.normSq (z i)) = 0)
    (k : ℕ) (A B C : ℂ)
    (h : ∀ j, tenBandMoment z k j = (tenBandFeatures p.a j : ℂ) * A +
      (tenBandFeatures p.b j : ℂ) * B + (tenBandFeatures p.c j : ℂ) * C) :
    layerMoment z p.a k = A ∧ layerMoment z p.b k = B ∧ layerMoment z p.c k = C := by
  have ha := hp.admissible
  have hs := support_of_radialPenalty_zero p _ hr
  have hsum (j : Fin 6) : tenBandMoment z k j = (tenBandFeatures p.a j : ℂ) * layerMoment z p.a k +
      (tenBandFeatures p.b j : ℂ) * layerMoment z p.b k +
      (tenBandFeatures p.c j : ℂ) * layerMoment z p.c k := by
    have hlocal (r : ℝ) : (∑ i ∈ radiusLayer z r, (tenBandFeatures (Complex.normSq (z i)) j : ℂ) *
        z i ^ k) = (tenBandFeatures r j : ℂ) * layerMoment z r k := by
      rw [layerMoment, Finset.mul_sum]
      exact Finset.sum_congr rfl fun i hi => by rw [mem_radiusLayer.mp hi]
    simpa only [tenBandMoment, hlocal] using (sum_three_layers ha.a_lt_b.ne
      (ha.a_lt_b.trans ha.b_lt_c).ne ha.b_lt_c.ne hs
      (fun i => (tenBandFeatures (Complex.normSq (z i)) j : ℂ) * z i ^ k)).symm
  obtain ⟨hba, hbb, hbc⟩ := tenBand_partition hp
  have he (j : Fin 6) := (hsum j).symm.trans (h j)
  have hA : layerMoment z p.a k = A := by
    simpa [tenBandFeatures, hba, hbb, hbc] using he 0
  have hBC : layerMoment z p.b k + layerMoment z p.c k = B + C := by
    simpa [tenBandFeatures, hba, hbb, hbc] using he 3
  let s := outerBandAffine
  have hBC' : (s p.b : ℂ) * layerMoment z p.b k + (s p.c : ℂ) * layerMoment z p.c k =
      (s p.b : ℂ) * B + (s p.c : ℂ) * C := by
    simpa [tenBandFeatures, hba, hbb, hbc, s] using he 4
  have hn : (s p.b : ℂ) - (s p.c : ℂ) ≠ 0 := by
    intro he
    have hreal : s p.b = s p.c := by exact_mod_cast sub_eq_zero.mp he
    have h := (div_left_inj' (by norm_num : (1271 / 1000 - 89 / 100 : ℝ) ≠ 0)).mp hreal
    linarith [ha.b_lt_c]
  have hpB : ((s p.b : ℂ) - (s p.c : ℂ)) * (layerMoment z p.b k - B) = 0 := by
    linear_combination hBC' - (s p.c : ℂ) * hBC
  have hB : layerMoment z p.b k = B := sub_eq_zero.mp ((mul_eq_zero.mp hpB).resolve_left hn)
  exact ⟨hA, hB, by linear_combination hBC - hB⟩

theorem ten_layer_moments_of_gram {p : Parameters} {z : Configuration 10} {pivot : Fin 6}
    (hp : p.InRootBox) (href : tenSecondReference p pivot ≠ 0)
    (hm : TenGramMoments p z pivot) : TenLayerMoments p z := by
  obtain ⟨ha, hb, hc⟩ := tenBand_layerMoments_eq hp hm.radial 1 0 0 0 (by simpa using hm.first)
  obtain ⟨η, hη⟩ := proportional_of_projected_zero pivot (tenBandMoment z 2)
    (tenSecondReference p) href hm.secondProjected
  refine ⟨hm.radial, ha, hb, hc, η, ?_⟩
  apply tenBand_layerMoments_eq hp hm.radial 2
  intro j
  rw [hη j, tenSecondReference]
  ring

private def SignedPair (w v x : ℂ) : Prop := x = w ∨ x = -w ∨ x = v ∨ x = -v

private theorem signedPair_neg {w v x : ℂ} (h : SignedPair w v x) : SignedPair w v (-x) := by
  rcases h with h | h | h | h <;> rw [h] <;> simp [SignedPair]

private theorem paired_points_subset {r : ℝ} {a b c d w v : ℂ}
    (hb : b = -a) (hd : d = -c)
    (ha : Complex.normSq a = r) (hc : Complex.normSq c = r)
    (hw : Complex.normSq w = r) (hv : Complex.normSq v = r)
    (hs : a ^ 2 + b ^ 2 + c ^ 2 + d ^ 2 = 2 * (w ^ 2 + v ^ 2))
    (hn : w ^ 2 + v ^ 2 ≠ 0) :
    SignedPair w v a ∧ SignedPair w v b ∧ SignedPair w v c ∧ SignedPair w v d := by
  rw [hb, hd] at hs
  have hs' : a ^ 2 + c ^ 2 = w ^ 2 + v ^ 2 := by linear_combination hs / 2
  have hne : a ^ 2 + c ^ 2 ≠ 0 := by simpa only [hs'] using hn
  have h : SignedPair w v a ∧ SignedPair w v c := by
    rcases antipodal_layer_eq ha hc hw hv hs' hne with ⟨h₁, h₂⟩ | ⟨h₁, h₂⟩
    · exact ⟨h₁.elim Or.inl (fun h => Or.inr (Or.inl h)), Or.inr ∘ Or.inr <| h₂⟩
    · exact ⟨Or.inr (Or.inr h₁), h₂.elim Or.inl (fun h => Or.inr (Or.inl h))⟩
  exact ⟨h.1, hb.symm ▸ signedPair_neg h.1, h.2, hd.symm ▸ signedPair_neg h.2⟩

private theorem four_layer_subset {r : ℝ} (z : Configuration 4) {w v : ℂ}
    (hnorm : ∀ i, Complex.normSq (z i) = r)
    (hw : Complex.normSq w = r) (hv : Complex.normSq v = r)
    (hfirst : ∑ i, z i = 0)
    (hsecond : ∑ i, z i ^ 2 = 2 * (w ^ 2 + v ^ 2))
    (hne : w ^ 2 + v ^ 2 ≠ 0) :
    ∀ i, SignedPair w v (z i) := by
  have hfirst' : z 0 + z 1 + z 2 + z 3 = 0 := by
    simpa [Fin.sum_univ_succ, add_assoc] using hfirst
  have hsecond' : z 0 ^ 2 + z 1 ^ 2 + z 2 ^ 2 + z 3 ^ 2 = 2 * (w ^ 2 + v ^ 2) := by
    simpa [Fin.sum_univ_succ, add_assoc] using hsecond
  have h :
      SignedPair w v (z 0) ∧ SignedPair w v (z 1) ∧ SignedPair w v (z 2) ∧ SignedPair w v (z 3) := by
    rcases four_points_antipodal (hnorm 0) (hnorm 1) (hnorm 2) (hnorm 3) hfirst' with
      ⟨h₁, h₂⟩ | ⟨h₁, h₂⟩ | ⟨h₁, h₂⟩
    · exact paired_points_subset h₁ h₂ (hnorm 0) (hnorm 2) hw hv hsecond' hne
    · have he : z 0 ^ 2 + z 2 ^ 2 + z 1 ^ 2 + z 3 ^ 2 = 2 * (w ^ 2 + v ^ 2) := by
        linear_combination hsecond'
      obtain ⟨h0, h2, h1, h3⟩ := paired_points_subset h₁ h₂ (hnorm 0) (hnorm 1) hw hv he hne
      exact ⟨h0, h1, h2, h3⟩
    · have he : z 0 ^ 2 + z 3 ^ 2 + z 1 ^ 2 + z 2 ^ 2 = 2 * (w ^ 2 + v ^ 2) := by
        linear_combination hsecond'
      obtain ⟨h0, h3, h1, h2⟩ := paired_points_subset h₁ h₂ (hnorm 0) (hnorm 1) hw hv he hne
      exact ⟨h0, h1, h2, h3⟩
  intro i
  fin_cases i
  · exact h.1
  · exact h.2.1
  · exact h.2.2.1
  · exact h.2.2.2

private theorem rectangle_square_sum {r t : ℝ} (h₁ : -r ≤ t) (h₂ : t ≤ r) :
    rectanglePoint r t ^ 2 + conj (rectanglePoint r t) ^ 2 = 2 * (t : ℂ) := by
  apply Complex.ext
  · simp only [pow_two, Complex.add_re, Complex.mul_re, Complex.conj_re, Complex.conj_im,
      Complex.ofReal_re, rectanglePoint]
    nlinarith [Real.sq_sqrt (show 0 ≤ (r + t) / 2 by linarith),
      Real.sq_sqrt (show 0 ≤ (r - t) / 2 by linarith), show (2 : ℂ).re = 2 by norm_num,
      show (2 : ℂ).im = 0 by norm_num, Complex.ofReal_im t]
  · simp only [pow_two, Complex.add_im, Complex.mul_im, Complex.conj_re, Complex.conj_im,
      Complex.ofReal_re, Complex.ofReal_im, rectanglePoint]
    norm_num
    ring

private theorem scaled_rectangle_subset {r t : ℝ} {s : ℂ} (z : Configuration 4)
    (hs : Complex.normSq s = 1) (h₁ : -r ≤ t) (h₂ : t ≤ r) (ht : t ≠ 0)
    (hnorm : ∀ i, Complex.normSq (z i) = r)
    (hfirst : ∑ i, z i = 0)
    (hsecond : ∑ i, z i ^ 2 = s ^ 2 * (4 * (t : ℂ))) :
    ∀ i, SignedPair (s * rectanglePoint r t) (s * conj (rectanglePoint r t)) (z i) := by
  have hr := normSq_rectanglePoint h₁ h₂
  have hrc (w : ℂ) (hw : Complex.normSq w = r) : Complex.normSq (s * w) = r := by
    rw [Complex.normSq_mul, hs, hw, one_mul]
  have hsum :
      (s * rectanglePoint r t) ^ 2 + (s * conj (rectanglePoint r t)) ^ 2 =
        s ^ 2 * (2 * (t : ℂ)) := by
    rw [← rectangle_square_sum h₁ h₂]
    ring
  have hne : (s * rectanglePoint r t) ^ 2 + (s * conj (rectanglePoint r t)) ^ 2 ≠ 0 := by
    have hs0 : s ≠ 0 := Complex.normSq_pos.mp (by rw [hs]; norm_num)
    rw [hsum]
    exact mul_ne_zero (pow_ne_zero _ hs0) (mul_ne_zero (by norm_num) (by exact_mod_cast ht))
  apply four_layer_subset z hnorm (hrc _ hr) (hrc _ (by simpa using hr)) hfirst _ hne
  rw [hsecond, hsum]
  ring

/-- The equality moments force the entire ten-point orbit, with the
    inner antipodal pair fixing the common phase. -/
theorem ten_orbit_of_layer_moments {p : Parameters} {z : Configuration 10}
    (hp : p.InRootBox) (hz : Function.Injective z) (hq : mass z = 9)
    (hm : TenLayerMoments p z) : SameOrbit z (candidate p) := by
  have ha := hp.admissible
  have hab := ha.a_lt_b.ne
  have hac := (ha.a_lt_b.trans ha.b_lt_c).ne
  have hbc := ha.b_lt_c.ne
  have hsupport := support_of_radialPenalty_zero p _ hm.radial
  obtain ⟨hA, hB, hC⟩ := populations_of_mass hp (count_three_layers hab hac hbc hsupport)
    ((mass_of_three_layers hab hac hbc hsupport).trans hq)
  let A := listedLayer hA
  let B := listedLayer hB
  let C := listedLayer hC
  let a0 : ℂ := (Real.sqrt p.a : ℝ)
  let b0 := rectanglePoint p.b p.x
  let c0 := rectanglePoint p.c p.y
  let s := A 0 / a0
  have ha0 : Complex.normSq a0 = p.a := by
    change Complex.normSq ((Real.sqrt p.a : ℝ) : ℂ) = p.a
    rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt ha.a_pos.le]
  have ha0Sq : a0 ^ 2 = (p.a : ℂ) := by
    change ((Real.sqrt p.a : ℝ) : ℂ) ^ 2 = (p.a : ℂ)
    exact_mod_cast Real.sq_sqrt ha.a_pos.le
  have ha0ne : a0 ≠ 0 := Complex.normSq_pos.mp (by rw [ha0]; exact ha.a_pos)
  have hs : Complex.normSq s = 1 := by
    rw [show s = A 0 / a0 from rfl, Complex.normSq_div, normSq_listedLayer hA 0,
      ha0, div_self ha.a_pos.ne']
  have hA0 : A 0 = s * a0 := (div_mul_cancel₀ _ ha0ne).symm
  have hA1 : A 1 = -(A 0) := by
    have h : A 0 + A 1 = 0 := by
      simpa [A, Fin.sum_univ_succ] using (moment_listedLayer hA 1).trans hm.first_a
    linear_combination h
  obtain ⟨η, hηa, hηb, hηc⟩ := hm.second
  have heta : η = s ^ 2 := by
    have h : A 0 ^ 2 + A 1 ^ 2 = η * (2 * (p.a : ℂ)) := by
      simpa [A, Fin.sum_univ_succ] using (moment_listedLayer hA 2).trans hηa
    rw [hA1, neg_sq, hA0, mul_pow, ha0Sq] at h
    have he : η * (2 * (p.a : ℂ)) = s ^ 2 * (2 * (p.a : ℂ)) := by
      linear_combination -h
    exact mul_right_cancel₀ (mul_ne_zero (by norm_num) (by exact_mod_cast ha.a_pos.ne')) he
  have hB' : ∀ i, SignedPair (s * b0) (s * conj b0) (B i) := by
    apply scaled_rectangle_subset B hs ha.x_lower.le ha.x_upper.le ha.x_ne (normSq_listedLayer hB)
    · simpa [B, pow_one] using (moment_listedLayer hB 1).trans hm.first_b
    · exact (moment_listedLayer hB 2).trans (hηb.trans (by rw [heta]))
  have hC' : ∀ i, SignedPair (s * c0) (s * conj c0) (C i) := by
    apply scaled_rectangle_subset C hs ha.y_lower.le ha.y_upper.le ha.y_ne (normSq_listedLayer hC)
    · simpa [C, pow_one] using (moment_listedLayer hC 1).trans hm.first_c
    · exact (moment_listedLayer hC 2).trans (hηc.trans (by rw [heta]))
  apply sameOrbit_of_range_subset hz hs
  rintro _ ⟨i, rfl⟩
  rcases hsupport i with hi | hi | hi
  · obtain ⟨j, hj⟩ := mem_range_listedLayer hA hi
    change A j = z i at hj
    rw [← hj]
    fin_cases j
    · exact ⟨0, hA0.symm⟩
    · exact ⟨1, by simp [scale, candidate, ← hA0, hA1, a0]⟩
  · obtain ⟨j, hj⟩ := mem_range_listedLayer hB hi
    change B j = z i at hj
    rw [← hj]
    rcases hB' j with h | h | h | h
    · exact ⟨2, by simpa [scale, candidate, b0] using h.symm⟩
    · exact ⟨3, by simpa [scale, candidate, b0] using h.symm⟩
    · exact ⟨4, by simpa [scale, candidate, b0, Complex.star_def] using h.symm⟩
    · exact ⟨5, by simpa [scale, candidate, b0, Complex.star_def] using h.symm⟩
  · obtain ⟨j, hj⟩ := mem_range_listedLayer hC hi
    change C j = z i at hj
    rw [← hj]
    rcases hC' j with h | h | h | h
    · exact ⟨6, by simpa [scale, candidate, c0] using h.symm⟩
    · exact ⟨7, by simpa [scale, candidate, c0] using h.symm⟩
    · exact ⟨8, by simpa [scale, candidate, c0, Complex.star_def] using h.symm⟩
    · exact ⟨9, by simpa [scale, candidate, c0, Complex.star_def] using h.symm⟩

end Mordell
