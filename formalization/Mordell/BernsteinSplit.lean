import Mordell.BernsteinTransform
import Mathlib.Data.Fin.Rev

/-! Exact de Casteljau midpoint subdivision, including a coordinatewise lower-control rule. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def midpointWeight (R : Type*) [Field R] {d : ℕ} (right : Bool) (i k : Fin (d + 1)) : R :=
  let a := if right then i.rev else i
  let b := if right then k.rev else k
  (a.val.choose b.val : R) / 2 ^ a.val

private theorem bernsteinBasis_half {d : ℕ} (k : Fin (d + 1)) (x : ℝ) :
    bernsteinBasis d k (x / 2) = ∑ i : Fin (d + 1), midpointWeight ℝ false i k * bernsteinBasis d i x := by
  let f (i : ℕ) : ℝ := (i.choose k.val : ℝ) / 2 ^ i *
    ((d.choose i : ℝ) * x ^ i * (1 - x) ^ (d - i))
  change (d.choose k.val : ℝ) * (x / 2) ^ k.val * (1 - x / 2) ^ (d - k.val) = ∑ i : Fin (d + 1), f i.val
  rw [Fin.sum_univ_eq_sum_range f]
  have hk : k.val ≤ d := Nat.le_of_lt_succ k.isLt
  have hzero : (∑ i ∈ Finset.range k.val, f i) = 0 := by
    exact Finset.sum_eq_zero fun i hi => by simp [f, Nat.choose_eq_zero_of_lt (Finset.mem_range.mp hi)]
  conv_rhs => rw [show d + 1 = k.val + (d - k.val + 1) by omega, Finset.sum_range_add, hzero, zero_add]
  rw [show 1 - x / 2 = x / 2 + (1 - x) by ring, add_pow, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i hi => ?_
  have hc : ((k.val + i).choose k.val : ℝ) * (d.choose (k.val + i) : ℝ) =
      (d.choose k.val : ℝ) * ((d - k.val).choose i : ℝ) := by
    exact_mod_cast (by
      rw [Nat.mul_comm, Nat.choose_mul (Nat.le_add_right _ _), Nat.add_sub_cancel_left] :
        (k.val + i).choose k.val * d.choose (k.val + i) = d.choose k.val * (d - k.val).choose i)
  have hs : d - (k.val + i) = d - k.val - i := by omega
  have hp : (x / 2) ^ k.val * (x / 2) ^ i = x ^ (k.val + i) / 2 ^ (k.val + i) := by rw [← pow_add, div_pow]
  simp only [f, hs]
  calc
    _ = (d.choose k.val : ℝ) * ((d - k.val).choose i : ℝ) * ((x / 2) ^ k.val * (x / 2) ^ i) *
        (1 - x) ^ (d - k.val - i) := by ring
    _ = _ := by rw [hp, ← hc]; ring

private theorem bernsteinBasis_reflect {d : ℕ} (i : Fin (d + 1)) (x : ℝ) :
    bernsteinBasis d i (1 - x) = bernsteinBasis d i.rev x := by
  have hi : i.val ≤ d := Nat.le_of_lt_succ i.isLt
  have hr : i.rev.val = d - i.val := by
    rw [Fin.val_rev]
    omega
  rw [bernsteinBasis, bernsteinBasis, hr, Nat.choose_symm hi,
    show d - (d - i.val) = i.val by omega, show 1 - (1 - x) = x by ring]
  ring

private theorem bernsteinBasis_right_half {d : ℕ} (k : Fin (d + 1)) (x : ℝ) :
    bernsteinBasis d k ((1 + x) / 2) = ∑ i : Fin (d + 1), midpointWeight ℝ true i k * bernsteinBasis d i x := by
  have h : bernsteinBasis d k ((1 + x) / 2) = bernsteinBasis d k.rev ((1 - x) / 2) := by
    rw [← bernsteinBasis_reflect, show 1 - (1 - x) / 2 = (1 + x) / 2 by ring]
  rw [h, bernsteinBasis_half]
  let f (i : Fin (d + 1)) := midpointWeight ℝ false i k.rev * bernsteinBasis d i (1 - x)
  change (∑ i, f i) = _
  rw [← Equiv.sum_comp Fin.revPerm f]
  refine Finset.sum_congr rfl fun i _ => ?_
  change midpointWeight ℝ false i.rev k.rev * bernsteinBasis d i.rev (1 - x) = _
  rw [bernsteinBasis_reflect, Fin.rev_rev]
  rfl

def midpointParameter (right : Bool) (x : ℝ) : ℝ := if right then (1 + x) / 2 else x / 2

theorem midpointWeight_sound {d : ℕ} (right : Bool) (k : Fin (d + 1)) (x : ℝ) :
    bernsteinBasis d k (midpointParameter right x) =
      ∑ i : Fin (d + 1), midpointWeight ℝ right i k * bernsteinBasis d i x := by
  cases right
  · exact bernsteinBasis_half k x
  · exact bernsteinBasis_right_half k x

def splitControls {n : ℕ} {R : Type*} [Field R] (d : Fin n → ℕ)
    (c : (∀ j, Fin (d j + 1)) → R) (axis : Fin n) (right : Bool) (i : ∀ j, Fin (d j + 1)) : R :=
  ∑ k : Fin (d axis + 1), midpointWeight R right (i axis) k * c (Function.update i axis k)

private def swapAxis {n : ℕ} (d : Fin n → ℕ) (axis : Fin n) : Equiv.Perm ((∀ j, Fin (d j + 1)) × Fin (d axis + 1)) :=
  Function.Involutive.toPerm (fun p => (Function.update p.1 axis p.2, p.1 axis)) (by
    intro p
    apply Prod.ext
    · funext j
      by_cases h : j = axis <;> simp [h]
    · simp)

theorem splitControls_sound {n : ℕ} (d : Fin n → ℕ) (c : (∀ j, Fin (d j + 1)) → ℝ)
    (axis : Fin n) (right : Bool) (x : Fin n → ℝ) :
    tensorBernstein d (splitControls d c axis right) x =
      tensorBernstein d c (Function.update x axis (midpointParameter right (x axis))) := by
  let p (i : ∀ j, Fin (d j + 1)) (j : Fin n) := bernsteinBasis (d j) (i j) (x j)
  let other (i : ∀ j, Fin (d j + 1)) := ∏ j ∈ Finset.univ.erase axis, p i j
  have hprod (i : ∀ j, Fin (d j + 1)) : other i * p i axis = ∏ j, p i j :=
    Finset.prod_erase_mul _ _ (Finset.mem_univ _)
  have hs (i : ∀ j, Fin (d j + 1)) (j : Fin (d axis + 1)) :
      other (Function.update i axis j) = other i := by
    apply Finset.prod_congr rfl
    intro k hk
    simp only [p, Function.update_of_ne (Finset.ne_of_mem_erase hk)]
  let f (ik : (∀ j, Fin (d j + 1)) × Fin (d axis + 1)) :=
    c ik.1 * midpointWeight ℝ right ik.2 (ik.1 axis) * (other ik.1 * bernsteinBasis (d axis) ik.2 (x axis))
  calc
    _ = ∑ i, ∑ k, f (swapAxis d axis (i, k)) := by
      unfold tensorBernstein splitControls
      simp only [Finset.sum_mul]
      refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun k _ => ?_
      change (midpointWeight ℝ right (i axis) k * c (Function.update i axis k)) * (∏ j, p i j) = _
      rw [← hprod i]
      have he : swapAxis d axis (i, k) = (Function.update i axis k, i axis) := rfl
      simp only [f, he, Function.update_self, hs, p]
      ring
    _ = ∑ i, ∑ k, f (i, k) := by
      rw [← Fintype.sum_prod_type (fun p => f (swapAxis d axis p)), ← Fintype.sum_prod_type f]
      exact Equiv.sum_comp (swapAxis d axis) f
    _ = _ := by
      unfold tensorBernstein
      refine Finset.sum_congr rfl fun i _ => ?_
      have ht : (∏ j, bernsteinBasis (d j) (i j)
          (Function.update x axis (midpointParameter right (x axis)) j)) =
          other i * bernsteinBasis (d axis) (i axis) (midpointParameter right (x axis)) := by
        rw [← Finset.prod_erase_mul _ _ (Finset.mem_univ axis), Function.update_self]
        congr 1
        exact Finset.prod_congr rfl fun j hj => by simp only [p, Function.update_of_ne (Finset.ne_of_mem_erase hj)]
      rw [ht, midpointWeight_sound]
      simp only [Finset.mul_sum]
      exact Finset.sum_congr rfl fun k _ => by simp only [f]; ring

theorem splitControls_lower {n : ℕ} (d : Fin n → ℕ) (c lower : (∀ j, Fin (d j + 1)) → ℝ)
    (axis : Fin n) (right : Bool) {x : Fin n → ℝ} (hx : ∀ j, 0 ≤ x j ∧ x j ≤ 1)
    (h : ∀ i, lower i ≤ splitControls d c axis right i) :
    tensorBernstein d lower x ≤
      tensorBernstein d c (Function.update x axis (midpointParameter right (x axis))) := by
  rw [← splitControls_sound]
  exact tensorBernstein_mono hx h

end Mordell
