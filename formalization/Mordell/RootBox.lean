import Mordell.Candidate
import Mordell.RootData
import Mathlib.Topology.MetricSpace.Contracting

/-! Exact existence and uniqueness of the root used in the preprint. -/

open scoped BigOperators Matrix

noncomputable section

namespace Mordell

private abbrev Coord := Fin 5 → ℝ

private def fromCentered (v : Coord) : Parameters :=
  ⟨(118037 + v 0) / 1048576, (1084060 + v 1) / 1048576,
    (1216217 + v 2) / 1048576, (-797427 + v 3) / 1048576,
    (810046 + v 4) / 1048576⟩

private def toCentered (p : Parameters) : Coord :=
  ![1048576 * p.a - 118037, 1048576 * p.b - 1084060,
    1048576 * p.c - 1216217, 1048576 * p.x + 797427,
    1048576 * p.y - 810046]

private theorem fromCentered_toCentered (p : Parameters) : fromCentered (toCentered p) = p := by
  rcases p with ⟨a, b, c, x, y⟩
  simp [fromCentered, toCentered]

private def unitBox : Set Coord := Metric.closedBall 0 1

private theorem unitBox_iff {v : Coord} : v ∈ unitBox ↔ ∀ i, |v i| ≤ 1 := by
  rw [unitBox, Metric.mem_closedBall, dist_pi_le_iff (by norm_num : (0 : ℝ) ≤ 1)]
  simp

private theorem fromCentered_in_box (v : Coord) : (fromCentered v).InRootBox ↔ v ∈ unitBox := by
  rw [unitBox_iff]
  simp only [Parameters.InRootBox, fromCentered, Set.mem_Icc]
  constructor
  · rintro ⟨⟨ha₁, ha₂⟩, ⟨hb₁, hb₂⟩, ⟨hc₁, hc₂⟩, ⟨hx₁, hx₂⟩, ⟨hy₁, hy₂⟩⟩ i
    have ha : |v 0| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
    have hb : |v 1| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
    have hc : |v 2| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
    have hx : |v 3| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
    have hy : |v 4| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
    fin_cases i <;> assumption
  · intro h
    have ha := abs_le.mp (h 0)
    have hb := abs_le.mp (h 1)
    have hc := abs_le.mp (h 2)
    have hx := abs_le.mp (h 3)
    have hy := abs_le.mp (h 4)
    obtain ⟨ha₁, ha₂⟩ := ha
    obtain ⟨hb₁, hb₂⟩ := hb
    obtain ⟨hc₁, hc₂⟩ := hc
    obtain ⟨hx₁, hx₂⟩ := hx
    obtain ⟨hy₁, hy₂⟩ := hy
    exact ⟨⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩,
      ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩⟩

private def equationVector (v : Coord) : Coord := fun k => (fromCentered v).equation k

private def rootMap (v : Coord) : Coord :=
  fun i => v i - 1048576 * (rootPreconditioner *ᵥ equationVector v) i

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 12000 in
private theorem rootMap_evaluate (v : Coord) (i : Fin 5) :
    rootMap v i = evaluate v (rootCoefficients i) := by
  fin_cases i
  all_goals first
    | change rootMap v 0 = evaluate v (rootCoefficients 0)
    | change rootMap v 1 = evaluate v (rootCoefficients 1)
    | change rootMap v 2 = evaluate v (rootCoefficients 2)
    | change rootMap v 3 = evaluate v (rootCoefficients 3)
    | change rootMap v 4 = evaluate v (rootCoefficients 4)
  all_goals
    norm_num [rootMap, rootPreconditioner, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
      equationVector, fromCentered, Parameters.equation, Parameters.powerSum,
      Parameters.weightedPowerSum, Parameters.quadPower, Finset.sum_range_succ,
      rootCoefficients, evaluate, monomialValue, Matrix.cons_val_two, Matrix.cons_val_three,
      Matrix.cons_val_four, Matrix.vecHead, Matrix.vecTail]
    ring

set_option maxRecDepth 12000 in
set_option maxHeartbeats 0 in
private theorem rootCoefficients_bounds :
    ∀ i : Fin 5, magnitude (rootCoefficients i) ≤ 3 / 4 ∧
      lipschitzBound (rootCoefficients i) ≤ 1 / 4 := by
  decide +kernel

private theorem rootMap_into_box {v : Coord} (hv : v ∈ unitBox) : rootMap v ∈ unitBox := by
  rw [unitBox_iff] at hv ⊢
  intro i
  rw [rootMap_evaluate]
  have hb : (magnitude (rootCoefficients i) : ℝ) ≤ 3 / 4 := by
    simpa only [Rat.cast_div, Rat.cast_ofNat] using
      (Rat.cast_le (K := ℝ)).2 (rootCoefficients_bounds i).1
  exact (abs_evaluate_le hv _).trans (hb.trans (by norm_num))

private theorem rootMap_dist {v w : Coord} (hv : v ∈ unitBox) (hw : w ∈ unitBox) :
    dist (rootMap v) (rootMap w) ≤ (1 / 4) * dist v w := by
  rw [unitBox_iff] at hv hw
  apply (dist_pi_le_iff (by positivity)).2
  intro i
  rw [Real.dist_eq, rootMap_evaluate, rootMap_evaluate]
  have hdist : ∀ i, |v i - w i| ≤ dist v w := by
    intro k
    simpa only [← Real.dist_eq] using (dist_pi_le_iff dist_nonneg).mp (le_rfl (a := dist v w)) k
  have hb : (lipschitzBound (rootCoefficients i) : ℝ) ≤ 1 / 4 := by
    simpa using (Rat.cast_le (K := ℝ)).2 (rootCoefficients_bounds i).2
  exact (abs_evaluate_sub_le hv hw dist_nonneg hdist _).trans
    (mul_le_mul_of_nonneg_right hb dist_nonneg)

private theorem rootPreconditioner_left_inv : rootLeftInverse * rootPreconditioner = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [rootLeftInverse, rootPreconditioner, Matrix.mul_apply, Matrix.one_apply, Fin.sum_univ_succ]

private theorem rootMap_fixed_iff (v : Coord) :
    rootMap v = v ↔ ∀ k : Fin 5, (fromCentered v).equation k = 0 := by
  constructor
  · intro h
    have hp : rootPreconditioner *ᵥ equationVector v = 0 := by
      funext i
      have hi := congrFun h i
      change v i - 1048576 * (rootPreconditioner *ᵥ equationVector v) i = v i at hi
      change (rootPreconditioner *ᵥ equationVector v) i = 0
      linarith
    have hq : rootLeftInverse *ᵥ (rootPreconditioner *ᵥ equationVector v) =
        rootLeftInverse *ᵥ (0 : Coord) := congrArg (fun w => rootLeftInverse *ᵥ w) hp
    rw [Matrix.mulVec_mulVec, rootPreconditioner_left_inv, Matrix.one_mulVec,
      Matrix.mulVec_zero] at hq
    exact fun k => congrFun hq k
  · intro h
    have hz : equationVector v = 0 := funext h
    funext i
    change v i - 1048576 * (rootPreconditioner *ᵥ equationVector v) i = v i
    rw [hz, Matrix.mulVec_zero]
    simp

/-- The preprint's five equations have a unique real solution in the
    stated dyadic box. Both rational bounds on the correction map are
    computed by kernel reduction. -/
theorem existsUnique_candidate_root : ∃! p : Parameters, p.IsRoot := by
  haveI : CompleteSpace unitBox :=
    (Metric.isClosed_closedBall : IsClosed unitBox).isComplete.completeSpace_coe
  haveI : Nonempty unitBox :=
    ⟨⟨0, unitBox_iff.2 (by simp)⟩⟩
  let g : unitBox → unitBox := fun v => ⟨rootMap v.val, rootMap_into_box v.property⟩
  have hg : ContractingWith (1 / 4) g := by
    refine ⟨by norm_num, LipschitzWith.of_dist_le_mul ?_⟩
    intro v w
    exact rootMap_dist v.property w.property
  let m := ContractingWith.fixedPoint g hg
  have hmf : rootMap m.val = m.val := congrArg Subtype.val hg.fixedPoint_isFixedPt
  refine ⟨fromCentered m.val, ⟨(fromCentered_in_box _).2 m.property, (rootMap_fixed_iff _).1 hmf⟩, ?_⟩
  intro p hp
  let v : unitBox := ⟨toCentered p, (fromCentered_in_box _).1 (by
    simpa only [fromCentered_toCentered] using hp.1)⟩
  have hv : Function.IsFixedPt g v := by
    apply Subtype.ext
    apply (rootMap_fixed_iff _).2
    simpa only [v, fromCentered_toCentered] using hp.2
  have he := congrArg (fun v : unitBox => fromCentered v.val) (hg.fixedPoint_unique hv)
  simpa only [v, fromCentered_toCentered, m] using he

noncomputable def isolatedParameters : Parameters :=
  Classical.choose existsUnique_candidate_root.exists

theorem isolatedParameters_isRoot : isolatedParameters.IsRoot :=
  Classical.choose_spec existsUnique_candidate_root.exists

end Mordell
