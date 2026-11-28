import Mordell.Basic
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.Topology.MetricSpace.Contracting

/-! A preconditioned root theorem using bounded real divided differences. -/

open scoped BigOperators Matrix NNReal

noncomputable section

namespace Mordell

theorem dist_mulVec_zero_le {n : ℕ} {A : Matrix (Fin n) (Fin n) ℝ} {q : ℝ}
    (hq : 0 ≤ q) (hA : ∀ i, ∑ j, |A i j| ≤ q) (x : Fin n → ℝ) :
    dist (A *ᵥ x) 0 ≤ q * dist x 0 := by
  apply (dist_pi_le_iff (mul_nonneg hq dist_nonneg)).2
  intro i
  have hx (j : Fin n) : |x j| ≤ dist x 0 := by
    simpa only [Pi.zero_apply, Real.dist_eq, sub_zero] using
      (dist_pi_le_iff (f := x) (g := 0) dist_nonneg).1 (le_refl (dist x 0)) j
  simp only [Pi.zero_apply, Real.dist_eq, sub_zero, Matrix.mulVec, dotProduct]
  calc
    _ ≤ ∑ j, |A i j * x j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ j, |A i j| * dist x 0 := Finset.sum_le_sum fun j _ => by
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (hx j) (abs_nonneg _)
    _ ≤ _ := by rw [← Finset.sum_mul]; exact mul_le_mul_of_nonneg_right (hA i) dist_nonneg

private theorem mulVec_injective_of_rows {n : ℕ} {D : Matrix (Fin n) (Fin n) ℝ} {q : ℝ}
    (hq₀ : 0 ≤ q) (hq₁ : q < 1) (hD : ∀ i, ∑ j, |(1 - D) i j| ≤ q) :
    Function.Injective D.mulVec := by
  intro x y h
  have he : (1 - D) *ᵥ (x - y) = x - y := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.mulVec_sub, h, sub_self, sub_zero]
  have hl := dist_mulVec_zero_le hq₀ hD (x - y)
  rw [he] at hl
  have hd : dist (x - y) 0 = 0 := by nlinarith [dist_nonneg (x := x - y) (y := 0)]
  exact sub_eq_zero.mp (dist_eq_zero.mp hd)

/-- If each secant of `f` is represented by a matrix whose preconditioned
    error has row sums below one, the residual test gives a unique zero. -/
theorem existsUnique_root_of_box_slope {n : ℕ} (f : (Fin n → ℝ) → Fin n → ℝ)
    (C : Matrix (Fin n) (Fin n) ℝ) (m : Fin n → ℝ) {r : ℝ} (hr : 0 ≤ r)
    {q : ℝ≥0} (hq : q < 1) (hres : dist (C *ᵥ f m) 0 ≤ (1 - (q : ℝ)) * r)
    (hs : ∀ x ∈ Metric.closedBall m r, ∀ y ∈ Metric.closedBall m r,
      ∃ D : Matrix (Fin n) (Fin n) ℝ, f x - f y = D *ᵥ (x - y) ∧
        ∀ i, ∑ j, |(1 - C * D) i j| ≤ (q : ℝ)) :
    ∃! x, x ∈ Metric.closedBall m r ∧ f x = 0 := by
  let box := Metric.closedBall m r
  let g (x : Fin n → ℝ) := x - C *ᵥ f x
  have hm : m ∈ box := Metric.mem_closedBall_self hr
  have hqr : (q : ℝ) < 1 := by exact_mod_cast hq
  obtain ⟨D, _, hD⟩ := hs m hm m hm
  have hCD := Matrix.mulVec_injective_iff_isUnit.mp (mulVec_injective_of_rows q.coe_nonneg hqr hD)
  have hC : Function.Injective C.mulVec := by
    apply Matrix.mulVec_injective_iff_isUnit.mpr
    apply (Matrix.isUnit_iff_isUnit_det C).mpr
    have hdet := (Matrix.isUnit_iff_isUnit_det (C * D)).mp hCD
    rw [Matrix.det_mul] at hdet
    exact isUnit_of_mul_isUnit_left hdet
  have hdist (x y : Fin n → ℝ) (hx : x ∈ box) (hy : y ∈ box) :
      dist (g x) (g y) ≤ q * dist x y := by
    obtain ⟨E, hE, hrow⟩ := hs x hx y hy
    have hge : g x - g y = (1 - C * E) *ᵥ (x - y) := by
      rw [Matrix.sub_mulVec, Matrix.one_mulVec, ← Matrix.mulVec_mulVec, ← hE, Matrix.mulVec_sub]
      dsimp [g]
      abel
    rw [dist_eq_norm, hge]
    simpa only [dist_eq_norm, sub_zero] using dist_mulVec_zero_le q.coe_nonneg hrow (x - y)
  have hg (x : box) : g x.val ∈ box := by
    apply Metric.mem_closedBall.mpr
    have hm' : dist (g m) m ≤ (1 - (q : ℝ)) * r := by
      simpa [g, dist_eq_norm] using hres
    calc
      _ ≤ dist (g x.val) (g m) + dist (g m) m := dist_triangle _ _ _
      _ ≤ (q : ℝ) * r + (1 - (q : ℝ)) * r := add_le_add
        ((hdist _ _ x.property hm).trans (mul_le_mul_of_nonneg_left x.property q.coe_nonneg)) hm'
      _ = r := by ring
  haveI : CompleteSpace box := Metric.isClosed_closedBall.isComplete.completeSpace_coe
  haveI : Nonempty box := ⟨⟨m, hm⟩⟩
  let G (x : box) : box := ⟨g x.val, hg x⟩
  have hG : ContractingWith q G :=
    ⟨hq, LipschitzWith.of_dist_le_mul (fun x y => hdist _ _ x.property y.property)⟩
  have hfixed (x : box) : Function.IsFixedPt G x ↔ f x.val = 0 := by
    constructor
    · intro h
      have hh := congrArg Subtype.val h
      change x.val - C *ᵥ f x.val = x.val at hh
      apply hC
      rw [Matrix.mulVec_zero]
      simpa using congrArg (fun y => x.val - y) hh
    · intro h
      apply Subtype.ext
      change x.val - C *ᵥ f x.val = x.val
      rw [h, Matrix.mulVec_zero, sub_zero]
  let root := ContractingWith.fixedPoint G hG
  refine ⟨root.val, ⟨root.property, (hfixed root).mp hG.fixedPoint_isFixedPt⟩, ?_⟩
  rintro x ⟨hx, hf⟩
  exact congrArg Subtype.val (hG.fixedPoint_unique ((hfixed ⟨x, hx⟩).mpr hf))

end Mordell
