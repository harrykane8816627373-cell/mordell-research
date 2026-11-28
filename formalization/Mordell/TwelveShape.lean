import Mordell.TwelveSymmetry
import Mordell.TriangleRigidity
import Mordell.Layers

/-! Identification of the isolated stationary configuration with the four-triangle model. -/

open scoped ComplexConjugate

noncomputable section

namespace Mordell.TwelveRootBox

private abbrev z := exactTwelveConfiguration

private def rotationPhase : ℂ := conj (symmetryImage false 0) / (‖symmetryImage false 0‖ : ℝ)

private def reflectionPhase : ℂ := conj (symmetryImage true 0) / (‖symmetryImage true 0‖ : ℝ)

private theorem rotation_identity (i : Fin 12) : rotationPhase * z (rotationPerm i) = z i :=
  congrFun (regauged_symmetry_fixed false) i

private theorem reflection_identity (i : Fin 12) : reflectionPhase * conj (z (reflectionPerm i)) = z i :=
  congrFun (regauged_symmetry_fixed true) i

private theorem rotation_order (i : Fin 12) : rotationPerm (rotationPerm (rotationPerm i)) = i := by
  fin_cases i <;> decide +kernel

private theorem pivot_nonzero : z 0 ≠ 0 := by
  have h : 0 < (z 0).re := pivot_pos_on_box isolatedCoordinates_spec.1
  intro he
  simp [he] at h

private theorem rotationPhase_cube : rotationPhase ^ 3 = 1 := by
  have h : rotationPhase ^ 3 * z 0 = 1 * z 0 := by
    calc
      _ = rotationPhase * (rotationPhase * (rotationPhase * z (rotationPerm (rotationPerm (rotationPerm 0))))) := by
        rw [rotation_order]
        ring
      _ = _ := by rw [rotation_identity, rotation_identity, rotation_identity, one_mul]
  exact mul_right_cancel₀ pivot_nonzero h

private theorem rotation_cubes (i : Fin 12) : z (rotationPerm i) ^ 3 = z i ^ 3 := by
  have h := congrArg (fun v : ℂ => v ^ 3) (rotation_identity i)
  simpa only [mul_pow, rotationPhase_cube, one_mul] using h

private theorem axis_nonzero : z 7 ≠ 0 := by
  let j : Fin 11 := 6
  let k : Fin 23 := (finProdFinEquiv (j, (0 : Fin 2))).succ
  have h : |isolatedCoordinates k - center k| ≤ radius := by
    simpa only [Real.dist_eq] using
      (dist_pi_le_iff radius_pos.le).1 (Metric.mem_closedBall.mp isolatedCoordinates_spec.1) k
  have hc : radius < center k := (Rat.cast_lt (K := ℝ)).2 (by decide +kernel : radiusRat < centerRat k)
  have hr : 0 < (z 7).re := by
    change 0 < isolatedCoordinates k
    have h' := (abs_le.mp h).1
    linarith
  intro he
  simp [he] at hr

def axisPhase : ℂ := conj (z 7) / (‖z 7‖ : ℝ)

theorem axisPhase_unit : Complex.normSq axisPhase = 1 := by
  rw [axisPhase, Complex.normSq_div, Complex.normSq_conj, Complex.normSq_eq_norm_sq,
    Complex.normSq_ofReal, ← sq, div_self (pow_ne_zero 2 (norm_ne_zero_iff.mpr axis_nonzero))]

private theorem axisPhase_reflection : axisPhase * reflectionPhase = conj axisPhase := by
  have h := reflection_identity 7
  change reflectionPhase * conj (z 7) = z 7 at h
  rw [axisPhase, div_mul_eq_mul_div, mul_comm (conj (z 7)), h, map_div₀, Complex.conj_conj,
    Complex.conj_ofReal]

def axisAligned : Configuration 12 := scale axisPhase massElevenTwelve

private theorem aligned_reflection (i : Fin 12) : axisAligned (reflectionPerm i) = conj (axisAligned i) := by
  have hd : conj downScaleTwelve = downScaleTwelve := Complex.conj_ofReal _
  have he : axisAligned i = conj (axisAligned (reflectionPerm i)) := by
    simp only [axisAligned, massElevenTwelve, scale, map_mul, hd]
    change axisPhase * (downScaleTwelve * z i) = conj axisPhase * (downScaleTwelve * conj (z (reflectionPerm i)))
    rw [← axisPhase_reflection, ← reflection_identity i]
    ring
  simpa only [Complex.conj_conj] using congrArg (starRingEnd ℂ) he.symm

private theorem aligned_real (i : Fin 12) (hi : reflectionPerm i = i) : (axisAligned i).im = 0 := by
  have h := congrArg Complex.im (aligned_reflection i)
  rw [hi, Complex.conj_im] at h
  linarith

def isolatedTriangleParameters : TriangleParameters :=
  ⟨(axisAligned 7).re, -(axisAligned 3).re, (axisAligned 0).re, (axisAligned 0).im⟩

def seedIndex : Fin 4 → Fin 12 := ![7, 3, 0, 11]

theorem aligned_seeds (j : Fin 4) : isolatedTriangleParameters.seeds j = axisAligned (seedIndex j) := by
  fin_cases j
  · apply Complex.ext <;> simp [TriangleParameters.seeds, isolatedTriangleParameters, seedIndex,
      aligned_real 7 (by decide +kernel)]
  · apply Complex.ext <;> simp [TriangleParameters.seeds, isolatedTriangleParameters, seedIndex,
      aligned_real 3 (by decide +kernel)]
  · apply Complex.ext <;> rfl
  · exact (aligned_reflection 0).symm

private theorem cycles_cover (i : Fin 12) : ∃ j : Fin 4, i = seedIndex j ∨
    i = rotationPerm (seedIndex j) ∨ i = rotationPerm (rotationPerm (seedIndex j)) := by
  fin_cases i <;> decide +kernel

private theorem aligned_cubes (i : Fin 12) : (axisAligned (rotationPerm i)) ^ 3 = (axisAligned i) ^ 3 := by
  simp only [axisAligned, massElevenTwelve, scale, mul_pow, rotation_cubes]

private theorem aligned_mem_candidate (i : Fin 12) : axisAligned i ∈ Set.range (candidateTwelve isolatedTriangleParameters) := by
  obtain ⟨j, hj⟩ := cycles_cover i
  have hc : (axisAligned i) ^ 3 = (isolatedTriangleParameters.seeds j) ^ 3 := by
    rw [aligned_seeds]
    rcases hj with h | h | h <;> rw [h] <;> simp only [aligned_cubes]
  have h (k : Fin 3) : candidateTwelve isolatedTriangleParameters (finProdFinEquiv (j, k)) =
      omega ^ k.val * isolatedTriangleParameters.seeds j := by
    simp only [candidateTwelve, triangleOrbits, Equiv.symm_apply_apply]
    exact mul_comm _ _
  rcases cube_roots_same_orbit hc with hh | hh | hh
  · exact ⟨finProdFinEquiv (j, (0 : Fin 3)), by simpa [h 0] using hh.symm⟩
  · exact ⟨finProdFinEquiv (j, (1 : Fin 3)), by simpa [h 1] using hh.symm⟩
  · exact ⟨finProdFinEquiv (j, (2 : Fin 3)), by simpa [h 2] using hh.symm⟩

/-- The archived solution, normalized to mass eleven, lies in the orbit of a four-triangle candidate. -/
theorem isolated_sameOrbit_candidate : SameOrbit massElevenTwelve (candidateTwelve isolatedTriangleParameters) := by
  have haxis : axisPhase ≠ 0 := Complex.normSq_pos.mp (by rw [axisPhase_unit]; norm_num)
  have hdown : downScaleTwelve ≠ 0 := by
    unfold downScaleTwelve
    exact Complex.ofReal_ne_zero.mpr (Real.sqrt_ne_zero'.mpr (by norm_num))
  have hinj : Function.Injective massElevenTwelve := fun _ _ h =>
    exactTwelveConfiguration_injective (mul_left_cancel₀ hdown h)
  have hu : Complex.normSq axisPhase⁻¹ = 1 := by rw [map_inv₀, axisPhase_unit, inv_one]
  apply sameOrbit_of_range_subset hinj hu
  rintro _ ⟨i, rfl⟩
  obtain ⟨j, hj⟩ := aligned_mem_candidate i
  refine ⟨j, ?_⟩
  change axisPhase⁻¹ * candidateTwelve isolatedTriangleParameters j = massElevenTwelve i
  rw [hj, axisAligned, scale, ← mul_assoc, inv_mul_cancel₀ haxis, one_mul]

theorem isolatedTriangleParameters_mass : isolatedTriangleParameters.RadiusMass := by
  rw [TriangleParameters.RadiusMass, ← mass_candidateTwelve, ← isolated_sameOrbit_candidate.mass_eq]
  exact massElevenTwelve_mass

end Mordell.TwelveRootBox
