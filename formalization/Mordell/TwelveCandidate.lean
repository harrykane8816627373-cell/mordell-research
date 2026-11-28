import Mordell.TwelveShape
import Mordell.TwelveEquality

/-! Algebraic nondegeneracy checks at the isolated twelve-point candidate. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell.TwelveRootBox

private abbrev p := isolatedTriangleParameters
private abbrev z := exactTwelveConfiguration

private theorem ball_upper {b : RationalBall} {x : ℝ} {u : ℚ} (hx : b.Contains x)
    (hu : b.mid + b.rad < u) : x < (u : ℝ) := by
  have h : |x - (b.mid : ℝ)| ≤ b.rad := hx
  have he := (abs_le.mp h).2
  have hu' : (b.mid : ℝ) + b.rad < (u : ℝ) := by exact_mod_cast hu
  linarith

private theorem ball_lower {b : RationalBall} {x : ℝ} {l : ℚ} (hx : b.Contains x)
    (hl : l < b.mid - b.rad) : (l : ℝ) < x := by
  have h : |x - (b.mid : ℝ)| ≤ b.rad := hx
  have he := (abs_le.mp h).1
  have hl' : (l : ℝ) < (b.mid : ℝ) - b.rad := by exact_mod_cast hl
  linarith

private def radiusBall (i : Fin 12) : RationalBall :=
  ((rootComponentBall i 0).mul (rootComponentBall i 0)).add
    ((rootComponentBall i 1).mul (rootComponentBall i 1))

private theorem radiusBall_contains (i : Fin 12) : (radiusBall i).Contains (Complex.normSq (z i)) := by
  have x := rootComponentBall_contains i 0
  have y := rootComponentBall_contains i 1
  simpa [radiusBall, complexComponent, Complex.normSq_apply] using
    RationalBall.add_sound (RationalBall.mul_sound x x) (RationalBall.mul_sound y y)

private theorem aligned_radii (i : Fin 12) : Complex.normSq (axisAligned i) = (11 / 12 : ℝ) * Complex.normSq (z i) := by
  simp [axisAligned, massElevenTwelve, scale, Complex.normSq_mul, axisPhase_unit, normSq_downScaleTwelve]

theorem isolated_radii_order : p.a < p.b ∧ p.b < p.c := by
  have a : Complex.normSq (z 7) < 1 := by
    simpa using ball_upper (radiusBall_contains 7) (by decide +kernel : (radiusBall 7).mid + (radiusBall 7).rad < 1)
  have b₁ : 1 < Complex.normSq (z 3) := by
    simpa using ball_lower (radiusBall_contains 3) (by decide +kernel : (1 : ℚ) < (radiusBall 3).mid - (radiusBall 3).rad)
  have b₂ : Complex.normSq (z 3) < (5 / 4 : ℝ) := by
    simpa using ball_upper (radiusBall_contains 3)
      (by decide +kernel : (radiusBall 3).mid + (radiusBall 3).rad < 5 / 4)
  have c : (5 / 4 : ℝ) < Complex.normSq (z 0) := by
    simpa using ball_lower (radiusBall_contains 0)
      (by decide +kernel : (5 / 4 : ℚ) < (radiusBall 0).mid - (radiusBall 0).rad)
  have ha := (congrArg Complex.normSq (aligned_seeds 0)).trans (aligned_radii 7)
  have hb := (congrArg Complex.normSq (aligned_seeds 1)).trans (aligned_radii 3)
  have hc := (congrArg Complex.normSq (aligned_seeds 2)).trans (aligned_radii 0)
  change Complex.normSq (p.inner : ℂ) = _ at ha
  change Complex.normSq (-(p.middle : ℂ)) = _ at hb
  change p.c = _ at hc
  simp only [Complex.normSq_neg, Complex.normSq_ofReal, ← sq] at ha hb
  change p.a = _ at ha
  change p.b = _ at hb
  constructor <;> linarith

private def cubicReBall (i : Fin 12) : RationalBall :=
  let x := rootComponentBall i 0
  let y := rootComponentBall i 1
  ((x.mul x).mul x).sub (((RationalBall.point 3).mul x).mul (y.mul y))

private theorem cubicReBall_contains (i : Fin 12) : (cubicReBall i).Contains ((z i) ^ 3).re := by
  have x := rootComponentBall_contains i 0
  have y := rootComponentBall_contains i 1
  have he : ((z i) ^ 3).re = (z i).re * (z i).re * (z i).re - 3 * (z i).re * ((z i).im * (z i).im) := by
    simp only [show (3 : ℕ) = 2 + 1 by decide, pow_add, pow_two, pow_one, Complex.mul_re, Complex.mul_im]
    ring
  rw [he]
  simpa [cubicReBall, complexComponent] using
    RationalBall.sub_sound (RationalBall.mul_sound (RationalBall.mul_sound x x) x)
      (RationalBall.mul_sound (RationalBall.mul_sound (RationalBall.point_sound 3) x) (RationalBall.mul_sound y y))

private def outerCubicBall : RationalBall := (cubicReBall 0).add (cubicReBall 11)

private theorem outerCubes_ne : (z 0) ^ 3 + (z 11) ^ 3 ≠ 0 := by
  have h : outerCubicBall.Contains ((z 0) ^ 3 + (z 11) ^ 3).re := by
    simpa only [outerCubicBall, Complex.add_re] using
      RationalBall.add_sound (cubicReBall_contains 0) (cubicReBall_contains 11)
  have hp : 0 < ((z 0) ^ 3 + (z 11) ^ 3).re := by
    simpa using ball_lower h (by decide +kernel : (0 : ℚ) < outerCubicBall.mid - outerCubicBall.rad)
  intro he
  simp [he] at hp

private def cubicSumBall : RationalBall := ((cubicReBall 7).add (cubicReBall 3)).add outerCubicBall

private theorem allCubes_ne : (z 7) ^ 3 + (z 3) ^ 3 + ((z 0) ^ 3 + (z 11) ^ 3) ≠ 0 := by
  have h : cubicSumBall.Contains ((z 7) ^ 3 + (z 3) ^ 3 + ((z 0) ^ 3 + (z 11) ^ 3)).re := by
    simpa only [cubicSumBall, outerCubicBall, Complex.add_re] using
      RationalBall.add_sound (RationalBall.add_sound (cubicReBall_contains 7) (cubicReBall_contains 3))
        (RationalBall.add_sound (cubicReBall_contains 0) (cubicReBall_contains 11))
  have hp : 0 < ((z 7) ^ 3 + (z 3) ^ 3 + ((z 0) ^ 3 + (z 11) ^ 3)).re := by
    simpa using ball_lower h (by decide +kernel : (0 : ℚ) < cubicSumBall.mid - cubicSumBall.rad)
  intro he
  simp [he] at hp

private theorem aligned_factor_ne : (axisPhase * downScaleTwelve) ^ 3 ≠ 0 := by
  have h₁ : axisPhase ≠ 0 := Complex.normSq_pos.mp (by rw [axisPhase_unit]; norm_num)
  have h₂ : downScaleTwelve ≠ 0 := Complex.normSq_pos.mp (by rw [normSq_downScaleTwelve]; norm_num)
  exact pow_ne_zero _ (mul_ne_zero h₁ h₂)

theorem isolated_outerCubeSum_ne : p.outerCubeSum ≠ 0 := by
  have h : p.outerCubeSum = (axisPhase * downScaleTwelve) ^ 3 * ((z 0) ^ 3 + (z 11) ^ 3) := by
    change (p.seeds 2) ^ 3 + (p.seeds 3) ^ 3 = _
    simp only [aligned_seeds, axisAligned, massElevenTwelve, scale, mul_pow, seedIndex]
    change (axisPhase^3 * (downScaleTwelve^3 * (z 0)^3)) + (axisPhase^3 * (downScaleTwelve^3 * (z 11)^3)) = _
    ring
  rw [h]
  exact mul_ne_zero aligned_factor_ne outerCubes_ne

theorem isolated_thirdReference_ne : p.thirdReference 0 ≠ 0 := by
  have h₁ : p.thirdReference 0 = 3 * ((p.seeds 0) ^ 3 + (p.seeds 1) ^ 3 + ((p.seeds 2) ^ 3 + (p.seeds 3) ^ 3)) := by
    simp [TriangleParameters.thirdReference, TriangleParameters.seeds, TriangleParameters.outerCubeSum]
    ring
  rw [h₁]
  have h₂ : (p.seeds 0) ^ 3 + (p.seeds 1) ^ 3 + ((p.seeds 2) ^ 3 + (p.seeds 3) ^ 3) =
      (axisPhase * downScaleTwelve) ^ 3 * ((z 7) ^ 3 + (z 3) ^ 3 + ((z 0) ^ 3 + (z 11) ^ 3)) := by
    simp only [aligned_seeds, axisAligned, massElevenTwelve, scale, seedIndex]
    change (axisPhase * (downScaleTwelve * z 7))^3 + (axisPhase * (downScaleTwelve * z 3))^3 +
      ((axisPhase * (downScaleTwelve * z 0))^3 + (axisPhase * (downScaleTwelve * z 11))^3) = _
    ring
  rw [h₂]
  exact mul_ne_zero (by norm_num) (mul_ne_zero aligned_factor_ne allCubes_ne)

theorem isolated_candidate_injective : Function.Injective (candidateTwelve p) := by
  have hs : downScaleTwelve ≠ 0 := Complex.normSq_pos.mp (by rw [normSq_downScaleTwelve]; norm_num)
  have hv : Function.Injective massElevenTwelve := fun _ _ h =>
    exactTwelveConfiguration_injective (mul_left_cancel₀ hs h)
  obtain ⟨s, e, _, h⟩ := isolated_sameOrbit_candidate
  intro i j hij
  apply e.symm.injective
  apply hv
  rw [h (e.symm i), h (e.symm j), Equiv.apply_symm_apply, Equiv.apply_symm_apply, hij]

private theorem candidate_at (j : Fin 4) (k : Fin 3) :
    candidateTwelve p (finProdFinEquiv (j, k)) = omega ^ k.val * p.seeds j := by
  simp only [candidateTwelve, triangleOrbits, Equiv.symm_apply_apply]
  exact mul_comm _ _

theorem isolated_seeds_nonzero (j : Fin 4) : p.seeds j ≠ 0 := by
  intro h
  have he : candidateTwelve p (finProdFinEquiv (j, (0 : Fin 3))) =
      candidateTwelve p (finProdFinEquiv (j, (1 : Fin 3))) := by
    rw [candidate_at j 0, candidate_at j 1, h, mul_zero, mul_zero]
  have hi := finProdFinEquiv.injective (isolated_candidate_injective he)
  exact (by decide : (0 : Fin 3) ≠ 1) (congrArg Prod.snd hi)

theorem isolated_cube_injective : Function.Injective (fun j => p.seeds j ^ 3) := by
  intro i j hij
  have hx : ∃ k : Fin 3, p.seeds i = omega ^ k.val * p.seeds j := by
    rcases cube_roots_same_orbit hij with h | h | h
    · exact ⟨0, by simpa using h⟩
    · exact ⟨1, by simpa using h⟩
    · exact ⟨2, by simpa using h⟩
  obtain ⟨k, hk⟩ := hx
  have he : candidateTwelve p (finProdFinEquiv (i, (0 : Fin 3))) =
      candidateTwelve p (finProdFinEquiv (j, k)) := by
    rw [candidate_at i 0, candidate_at j k]
    simpa using hk
  exact congrArg Prod.fst (finProdFinEquiv.injective (isolated_candidate_injective he))

theorem isolated_sameOrbit_scaled_candidate :
    SameOrbit exactTwelveConfiguration (scale upScaleTwelve (candidateTwelve p)) := by
  have hmul : upScaleTwelve * downScaleTwelve = 1 := by
    change ((Real.sqrt (12 / 11 : ℝ) : ℝ) : ℂ) * ((Real.sqrt (11 / 12 : ℝ) : ℝ) : ℂ) = 1
    rw [← Complex.ofReal_mul, ← Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 12 / 11)]
    norm_num
  have hscale : scale upScaleTwelve massElevenTwelve = exactTwelveConfiguration := by
    funext i
    simp [massElevenTwelve, scale, ← mul_assoc, hmul]
  simpa only [hscale] using isolated_sameOrbit_candidate.scale upScaleTwelve

theorem isolated_sharp_of_gram_bounds
    (h : ∀ z : Configuration 12, IsMaximizer 11 z → Function.Injective z →
      Nonempty {cert : BlockBound z (candidateTwelve p) //
        (∀ k a, ∑ i, cert.feature k i a = 0) → TwelveGramMoments p z}) :
    IsMaximizer 12 exactTwelveConfiguration ∧
      (∀ z : Configuration 12, IsMaximizer 12 z → SameOrbit z exactTwelveConfiguration) := by
  obtain ⟨hab, hbc⟩ := isolated_radii_order
  obtain ⟨hm, hu⟩ := sharp_twelve_of_gram_bounds p isolatedTriangleParameters_mass hab.ne (hab.trans hbc).ne hbc.ne
    isolated_thirdReference_ne isolated_outerCubeSum_ne isolated_seeds_nonzero isolated_cube_injective h
  have hOrbit := isolated_sameOrbit_scaled_candidate
  refine ⟨⟨hOrbit.mass_eq.trans hm.1, fun w hw => (hm.2 w hw).trans_eq hOrbit.distanceProduct_eq.symm⟩, ?_⟩
  intro w hw
  obtain ⟨s, e, hs, he⟩ := hu w hw
  obtain ⟨r, f, hr, hf⟩ := hOrbit
  have hr0 : r ≠ 0 := Complex.normSq_pos.mp (by rw [hr]; norm_num)
  refine ⟨s / r, e.trans f.symm, by rw [Complex.normSq_div, hs, hr, div_one], fun i => ?_⟩
  rw [Equiv.trans_apply, hf, Equiv.apply_symm_apply, ← mul_assoc, div_mul_cancel₀ _ hr0, he]

end Mordell.TwelveRootBox
