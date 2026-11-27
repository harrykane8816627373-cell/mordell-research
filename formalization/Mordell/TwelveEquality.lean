import Mordell.Layers
import Mordell.SixPointLayer
import Mordell.HarmonicKernel

/-! Twelve-point equality reconstruction from radial and angular moments. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

namespace TriangleParameters

def a (p : TriangleParameters) : ℝ := p.inner ^ 2

def b (p : TriangleParameters) : ℝ := p.middle ^ 2

def c (p : TriangleParameters) : ℝ := Complex.normSq p.outer

def outerCubeSum (p : TriangleParameters) : ℂ := p.outer ^ 3 + (conj p.outer) ^ 3

def thirdReference (p : TriangleParameters) (d : Fin 3) : ℂ :=
  ((p.a - 1 : ℝ) : ℂ) ^ d.val * (3 * (p.inner : ℂ) ^ 3) +
  ((p.b - 1 : ℝ) : ℂ) ^ d.val * (3 * (-(p.middle : ℂ)) ^ 3) +
  ((p.c - 1 : ℝ) : ℂ) ^ d.val * (3 * p.outerCubeSum)

end TriangleParameters

structure TwelveRadialMoments (p : TriangleParameters) (z : Configuration 12) : Prop where
  radial : ∑ i, ((Complex.normSq (z i) - p.a) * (Complex.normSq (z i) - p.b) *
    (Complex.normSq (z i) - p.c)) ^ 2 = 0
  radialSecond : ∑ i, Complex.normSq (z i) ^ 2 = 3 * p.a ^ 2 + 3 * p.b ^ 2 + 6 * p.c ^ 2

structure TwelveLayerMoments (p : TriangleParameters) (z : Configuration 12) : Prop
    extends TwelveRadialMoments p z where
  innerVanishing : layerMoment z p.a 1 = 0 ∧ layerMoment z p.a 2 = 0
  middleVanishing : layerMoment z p.b 1 = 0 ∧ layerMoment z p.b 2 = 0
  outerVanishing : layerMoment z p.c 1 = 0 ∧ layerMoment z p.c 2 = 0 ∧
    layerMoment z p.c 4 = 0 ∧ layerMoment z p.c 5 = 0
  third : ∃ η : ℂ,
    layerMoment z p.a 3 = η * (3 * (p.inner : ℂ) ^ 3) ∧
    layerMoment z p.b 3 = η * (3 * (-(p.middle : ℂ)) ^ 3) ∧
    layerMoment z p.c 3 = η * (3 * p.outerCubeSum)

/-- The radial zero moments and the first three radial coefficients in the relevant harmonics. -/
structure TwelveGramMoments (p : TriangleParameters) (z : Configuration 12) : Prop
    extends TwelveRadialMoments p z where
  harmonic : ∀ k : ℕ, k = 1 ∨ k = 2 ∨ k = 4 ∨ k = 5 →
    ∀ d : Fin 3, weightedMoment z 1 d.val k = 0
  thirdProjected : ∀ d : Fin 3, d ≠ 0 → weightedMoment z 1 d.val 3 -
    (p.thirdReference d / p.thirdReference 0) * weightedMoment z 1 0 3 = 0

private theorem support_of_twelve_moments {p : TriangleParameters} {z : Configuration 12}
    (hm : TwelveRadialMoments p z) (i : Fin 12) :
    Complex.normSq (z i) = p.a ∨ Complex.normSq (z i) = p.b ∨ Complex.normSq (z i) = p.c := by
  let q (j : Fin 12) := (Complex.normSq (z j) - p.a) * (Complex.normSq (z j) - p.b) *
    (Complex.normSq (z j) - p.c)
  have h : q i ^ 2 ≤ ∑ j, q j ^ 2 :=
    Finset.single_le_sum (fun j _ => sq_nonneg (q j)) (Finset.mem_univ i)
  have hs : ∑ j, q j ^ 2 = 0 := hm.radial
  rw [hs] at h
  have hz : q i = 0 := sq_eq_zero_iff.mp (le_antisymm h (sq_nonneg _))
  change (Complex.normSq (z i) - p.a) * (Complex.normSq (z i) - p.b) *
    (Complex.normSq (z i) - p.c) = 0 at hz
  rcases mul_eq_zero.mp hz with hab | hc
  · exact (mul_eq_zero.mp hab).imp sub_eq_zero.mp (Or.inl ∘ sub_eq_zero.mp)
  · exact Or.inr (Or.inr (sub_eq_zero.mp hc))

theorem twelve_layer_moments_of_gram {p : TriangleParameters} {z : Configuration 12}
    (hab : p.a ≠ p.b) (hac : p.a ≠ p.c) (hbc : p.b ≠ p.c)
    (href : p.thirdReference 0 ≠ 0) (hm : TwelveGramMoments p z) : TwelveLayerMoments p z := by
  have hs := support_of_twelve_moments hm.toTwelveRadialMoments
  have hvanish (k : ℕ) (hk : k = 1 ∨ k = 2 ∨ k = 4 ∨ k = 5) :
      layerMoment z p.a k = 0 ∧ layerMoment z p.b k = 0 ∧ layerMoment z p.c k = 0 :=
    layerMoments_eq_of_weightedMoments hab hac hbc hs 1 k 0 0 0 (by simpa using hm.harmonic k hk)
  have h1 := hvanish 1 (by decide)
  have h2 := hvanish 2 (by decide)
  have h4 := hvanish 4 (by decide)
  have h5 := hvanish 5 (by decide)
  obtain ⟨η, hη⟩ := proportional_of_projected_zero 0 (fun d : Fin 3 => weightedMoment z 1 d.val 3)
    p.thirdReference href hm.thirdProjected
  refine ⟨hm.toTwelveRadialMoments, ⟨h1.1, h2.1⟩, ⟨h1.2.1, h2.2.1⟩,
    ⟨h1.2.2, h2.2.2, h4.2.2, h5.2.2⟩, η, ?_⟩
  apply layerMoments_eq_of_weightedMoments hab hac hbc hs 1 3
  intro d
  rw [hη d, TriangleParameters.thirdReference]
  ring

private theorem three_layer_cubes {z : Configuration 3}
    (h1 : ∑ i, z i ^ 1 = 0) (h2 : ∑ i, z i ^ 2 = 0) :
    (∀ i, z i ^ 3 = z 0 ^ 3) ∧ (∑ i, z i ^ 3) = 3 * z 0 ^ 3 := by
  have he1 : z 0 + z 1 + z 2 = 0 := by simpa [Fin.sum_univ_succ, add_assoc] using h1
  have he2 : z 0 ^ 2 + z 1 ^ 2 + z 2 ^ 2 = 0 := by
    simpa [Fin.sum_univ_succ, add_assoc] using h2
  obtain ⟨ha, hb, hc⟩ := three_points_cubes_eq he1 he2
  have he (i : Fin 3) : z i ^ 3 = z 0 ^ 3 := by
    fin_cases i
    · rfl
    · exact hb.trans ha.symm
    · exact hc.trans ha.symm
  refine ⟨he, ?_⟩
  simp [he, nsmul_eq_mul]

private theorem mem_scaled_candidate_of_seed_cube {p : TriangleParameters} {s w : ℂ}
    (j : Fin 4) (hw : w ^ 3 = (s * p.seeds j) ^ 3) :
    w ∈ Set.range (scale s (candidateTwelve p)) := by
  have h (k : Fin 3) : scale s (candidateTwelve p) (finProdFinEquiv (j, k)) =
      omega ^ k.val * (s * p.seeds j) := by
    simp only [scale, candidateTwelve, triangleOrbits, Equiv.symm_apply_apply]
    ring
  rcases cube_roots_same_orbit hw with he | he | he
  · exact ⟨finProdFinEquiv (j, (0 : Fin 3)), by simpa [h 0] using he.symm⟩
  · exact ⟨finProdFinEquiv (j, (1 : Fin 3)), by simpa [h 1] using he.symm⟩
  · exact ⟨finProdFinEquiv (j, (2 : Fin 3)), by simpa [h 2] using he.symm⟩

/-- The three radial layers and their equality moments determine the twelve-point orbit. -/
theorem twelve_orbit_of_layer_moments {p : TriangleParameters} {z : Configuration 12}
    (hp : p.RadiusMass) (ha : p.inner ≠ 0) (hab : p.a ≠ p.b) (hac : p.a ≠ p.c) (hbc : p.b ≠ p.c)
    (houter : p.outerCubeSum ≠ 0) (hz : Function.Injective z) (hq : mass z = 11)
    (hm : TwelveLayerMoments p z) : SameOrbit z (candidateTwelve p) := by
  have hsupport := support_of_twelve_moments hm.toTwelveRadialMoments
  have hq' : 3 * p.a + 3 * p.b + 6 * p.c = 11 := by
    simpa [TriangleParameters.a, TriangleParameters.b, TriangleParameters.c, TriangleParameters.RadiusMass,
      TriangleParameters.outer, Complex.normSq_mk, pow_two] using hp
  have hr (r : ℝ) : (∑ i ∈ radiusLayer z r, Complex.normSq (z i) ^ 2) =
      ((radiusLayer z r).card : ℝ) * r ^ 2 := by
    calc
      _ = ∑ _i ∈ radiusLayer z r, r ^ 2 :=
        Finset.sum_congr rfl fun _ hi => congrArg (· ^ 2) (mem_radiusLayer.mp hi)
      _ = _ := by simp
  have hr2 := (sum_three_layers hab hac hbc hsupport (fun i => Complex.normSq (z i) ^ 2)).trans
    hm.radialSecond
  simp only [hr] at hr2
  obtain ⟨hA, hB, hC⟩ := populations_twelve_of_moments hab hac hbc
    (count_three_layers hab hac hbc hsupport)
    ((mass_of_three_layers hab hac hbc hsupport).trans (hq.trans hq'.symm)) hr2
  let A := listedLayer hA
  let B := listedLayer hB
  let C := listedLayer hC
  let a0 : ℂ := p.inner
  let b0 : ℂ := -(p.middle : ℂ)
  let s := A 0 / a0
  have ha0 : Complex.normSq a0 = p.a := by simp [a0, TriangleParameters.a, pow_two]
  have ha0ne : a0 ≠ 0 := by
    change (p.inner : ℂ) ≠ 0
    exact_mod_cast ha
  have hapos : 0 < p.a := ha0 ▸ Complex.normSq_pos.mpr ha0ne
  have hs : Complex.normSq s = 1 := by
    rw [show s = A 0 / a0 from rfl, Complex.normSq_div, normSq_listedLayer hA 0, ha0,
      div_self hapos.ne']
  have hs0 : s ≠ 0 := Complex.normSq_pos.mp (by rw [hs]; norm_num)
  have hA0 : A 0 = s * a0 := (div_mul_cancel₀ _ ha0ne).symm
  obtain ⟨hAi, hA3⟩ := three_layer_cubes ((moment_listedLayer hA 1).trans hm.innerVanishing.1)
    ((moment_listedLayer hA 2).trans hm.innerVanishing.2)
  obtain ⟨hBi, hB3⟩ := three_layer_cubes ((moment_listedLayer hB 1).trans hm.middleVanishing.1)
    ((moment_listedLayer hB 2).trans hm.middleVanishing.2)
  obtain ⟨j, _, hCi, hC3⟩ := six_layer_cubic_images (injective_listedLayer hC hz)
    ((moment_listedLayer hC 1).trans hm.outerVanishing.1)
    ((moment_listedLayer hC 2).trans hm.outerVanishing.2.1)
    ((moment_listedLayer hC 4).trans hm.outerVanishing.2.2.1)
    ((moment_listedLayer hC 5).trans hm.outerVanishing.2.2.2)
  obtain ⟨η, hηa, hηb, hηc⟩ := hm.third
  have hη : η = s ^ 3 := by
    have he := hA3.symm.trans ((moment_listedLayer hA 3).trans hηa)
    change 3 * A 0 ^ 3 = η * (3 * a0 ^ 3) at he
    rw [hA0, mul_pow] at he
    apply mul_right_cancel₀ (mul_ne_zero (by norm_num : (3 : ℂ) ≠ 0) (pow_ne_zero 3 ha0ne))
    linear_combination -he
  have hB0 : B 0 ^ 3 = (s * b0) ^ 3 := by
    have he := hB3.symm.trans ((moment_listedLayer hB 3).trans hηb)
    change 3 * B 0 ^ 3 = η * (3 * b0 ^ 3) at he
    rw [hη] at he
    rw [mul_pow]
    linear_combination he / 3
  have hCsum : C 0 ^ 3 + C j ^ 3 = (s * p.outer) ^ 3 + (s * conj p.outer) ^ 3 := by
    have he := hC3.symm.trans ((moment_listedLayer hC 3).trans hηc)
    change 3 * (C 0 ^ 3 + C j ^ 3) = η * (3 * p.outerCubeSum) at he
    rw [hη, TriangleParameters.outerCubeSum] at he
    linear_combination he / 3
  have hne : C 0 ^ 3 + C j ^ 3 ≠ 0 := by
    rw [hCsum, mul_pow, mul_pow, ← mul_add]
    exact mul_ne_zero (pow_ne_zero 3 hs0) houter
  have hO : Complex.normSq (s * p.outer) = p.c := by
    simp [Complex.normSq_mul, hs, TriangleParameters.c]
  have hObar : Complex.normSq (s * conj p.outer) = p.c := by
    simpa only [Complex.normSq_mul, Complex.normSq_conj] using hO
  have hout := cubic_layer_images_eq (normSq_listedLayer hC 0) (normSq_listedLayer hC j)
    hO hObar hCsum hne
  have hrange (i : Fin 6) : C i ∈ Set.range (scale s (candidateTwelve p)) := by
    rcases hout with ⟨h0, hj⟩ | ⟨h0, hj⟩ <;> rcases hCi i with hi | hi
    · exact mem_scaled_candidate_of_seed_cube 2 (hi.trans h0)
    · exact mem_scaled_candidate_of_seed_cube 3 (hi.trans hj)
    · exact mem_scaled_candidate_of_seed_cube 3 (hi.trans h0)
    · exact mem_scaled_candidate_of_seed_cube 2 (hi.trans hj)
  apply sameOrbit_of_range_subset hz hs
  rintro _ ⟨i, rfl⟩
  rcases hsupport i with hi | hi | hi
  · obtain ⟨k, hk⟩ := mem_range_listedLayer hA hi
    rw [← hk]
    exact mem_scaled_candidate_of_seed_cube 0 ((hAi k).trans (congrArg (· ^ 3) hA0))
  · obtain ⟨k, hk⟩ := mem_range_listedLayer hB hi
    rw [← hk]
    exact mem_scaled_candidate_of_seed_cube 1 ((hBi k).trans hB0)
  · obtain ⟨k, hk⟩ := mem_range_listedLayer hC hi
    rw [← hk]
    exact hrange k

theorem sharp_twelve_of_gram_bounds (p : TriangleParameters)
    (hp : p.RadiusMass) (hab : p.a ≠ p.b) (hac : p.a ≠ p.c) (hbc : p.b ≠ p.c)
    (href : p.thirdReference 0 ≠ 0) (houter : p.outerCubeSum ≠ 0) (hnz : ∀ j, p.seeds j ≠ 0)
    (hc : Function.Injective (fun j => p.seeds j ^ 3))
    (h : ∀ z : Configuration 12, IsMaximizer 11 z → Function.Injective z →
      Nonempty {cert : BlockBound z (candidateTwelve p) //
        (∀ k a, ∑ i, cert.feature k i a = 0) → TwelveGramMoments p z}) :
    IsMaximizer 12 (scale upScaleTwelve (candidateTwelve p)) ∧
      (∀ z : Configuration 12, IsMaximizer 12 z →
        SameOrbit z (scale upScaleTwelve (candidateTwelve p))) := by
  apply sharp_twelve_of_log_certificate p hp hnz hc
  intro z hz hi
  obtain ⟨⟨cert, he⟩⟩ := h z hz hi
  have ha : p.inner ≠ 0 := by
    intro ha
    exact hnz 0 (by simp [TriangleParameters.seeds, ha])
  have hmass : mass z = (↑(12 : ℕ) : ℝ) - 1 := hz.1.trans (by norm_num)
  obtain ⟨hb, hzero⟩ := cert.log_le_and_moments_zero (by decide) hi hmass
  refine ⟨hb, fun hlog => twelve_orbit_of_layer_moments hp ha hab hac hbc houter hi hz.1 ?_⟩
  exact twelve_layer_moments_of_gram hab hac hbc href (he (hzero hlog))

end Mordell
