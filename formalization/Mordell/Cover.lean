import Mordell.LogBounds

/-! Soundness of a box subdivision and of the tensor Bernstein sign test. -/

open scoped BigOperators

noncomputable section

namespace Mordell

structure Box (n : ℕ) where
  lower : Fin n → ℝ
  upper : Fin n → ℝ

namespace Box

def Contains {n : ℕ} (b : Box n) (x : Fin n → ℝ) : Prop :=
  ∀ i, b.lower i ≤ x i ∧ x i ≤ b.upper i

def midpoint {n : ℕ} (b : Box n) (i : Fin n) : ℝ :=
  (b.lower i + b.upper i) / 2

def left {n : ℕ} (b : Box n) (i : Fin n) : Box n :=
  ⟨b.lower, Function.update b.upper i (b.midpoint i)⟩

def right {n : ℕ} (b : Box n) (i : Fin n) : Box n :=
  ⟨Function.update b.lower i (b.midpoint i), b.upper⟩

theorem contained_in_children {n : ℕ} {b : Box n} {x : Fin n → ℝ}
    (h : b.Contains x) (i : Fin n) : (b.left i).Contains x ∨ (b.right i).Contains x := by
  by_cases hx : x i ≤ b.midpoint i
  · left
    intro j
    refine ⟨(h j).1, ?_⟩
    by_cases hji : j = i
    · subst j
      simpa [left] using hx
    · simpa [left, Function.update_of_ne hji] using (h j).2
  · right
    intro j
    refine ⟨?_, (h j).2⟩
    by_cases hji : j = i
    · subst j
      simpa [right] using (le_of_not_ge hx : b.midpoint i ≤ x i)
    · simpa [right, Function.update_of_ne hji] using (h j).1

end Box

/-- The analytic facts at leaves and the midpoint split rule compose to a
    whole-domain inequality. An external transcript needs a proof of every leaf. -/
inductive BoxCover {n : ℕ} (domain goal : (Fin n → ℝ) → Prop) : Box n → Prop
  | excluded {b} : (∀ x, b.Contains x → ¬ domain x) → BoxCover domain goal b
  | checked {b} : (∀ x, b.Contains x → domain x → goal x) → BoxCover domain goal b
  | split {b} (i : Fin n) :
      BoxCover domain goal (b.left i) → BoxCover domain goal (b.right i) →
        BoxCover domain goal b

theorem BoxCover.sound {n : ℕ} {domain goal : (Fin n → ℝ) → Prop}
    {root : Box n} (cover : BoxCover domain goal root) :
    ∀ x, root.Contains x → domain x → goal x := by
  induction cover with
  | excluded h => exact fun x hx hd => (h x hx hd).elim
  | checked h => exact h
  | split i _ _ hl hr =>
    intro x hx hd
    rcases Box.contained_in_children hx i with h | h
    · exact hl x h hd
    · exact hr x h hd

def bernsteinBasis (n : ℕ) (i : Fin (n + 1)) (x : ℝ) : ℝ :=
  (n.choose i.val : ℝ) * x ^ i.val * (1 - x) ^ (n - i.val)

def tensorBernstein {n : ℕ} (degree : Fin n → ℕ)
    (c : (∀ j, Fin (degree j + 1)) → ℝ) (x : Fin n → ℝ) : ℝ :=
  ∑ i, c i * ∏ j, bernsteinBasis (degree j) (i j) (x j)

theorem tensorBernstein_mono {n : ℕ} {degree : Fin n → ℕ}
    {a b : (∀ j, Fin (degree j + 1)) → ℝ} {x : Fin n → ℝ}
    (hx : ∀ j, 0 ≤ x j ∧ x j ≤ 1) (hc : ∀ i, a i ≤ b i) :
    tensorBernstein degree a x ≤ tensorBernstein degree b x := by
  classical
  unfold tensorBernstein
  apply Finset.sum_le_sum
  intro i _
  apply mul_le_mul_of_nonneg_right (hc i)
  apply Finset.prod_nonneg
  intro j _
  unfold bernsteinBasis
  exact mul_nonneg (mul_nonneg (Nat.cast_nonneg _) (pow_nonneg (hx j).1 _))
    (pow_nonneg (sub_nonneg.mpr (hx j).2) _)

end Mordell
