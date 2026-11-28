import Mordell.Cover

/-! Exact conversion from affine-rescaled monomials to tensor Bernstein controls. -/

open scoped BigOperators

noncomputable section

namespace Mordell

private def basisNat (d i : ℕ) (x : ℝ) : ℝ :=
  (d.choose i : ℝ) * x ^ i * (1 - x) ^ (d - i)

private theorem sum_basisNat (d : ℕ) (x : ℝ) :
    (∑ i ∈ Finset.range (d + 1), basisNat d i x) = 1 := by
  have h := (add_pow x (1 - x) d).symm
  rw [show x + (1 - x) = 1 by ring, one_pow] at h
  simpa only [basisNat, mul_left_comm, mul_comm, mul_assoc] using h

private theorem sum_choose_basisNat {d k : ℕ} (hk : k ≤ d) (x : ℝ) :
    (∑ i ∈ Finset.range (d + 1), (i.choose k : ℝ) * basisNat d i x) = (d.choose k : ℝ) * x ^ k := by
  let f (i : ℕ) := (i.choose k : ℝ) * basisNat d i x
  have hzero : (∑ i ∈ Finset.range k, f i) = 0 := by
    apply Finset.sum_eq_zero
    intro i hi
    simp [f, Nat.choose_eq_zero_of_lt (Finset.mem_range.mp hi)]
  have hpoint (i : ℕ) : f (k + i) = (d.choose k : ℝ) * x ^ k * basisNat (d - k) i x := by
    have hchoose : ((k + i).choose k : ℝ) * (d.choose (k + i) : ℝ) =
        (d.choose k : ℝ) * ((d - k).choose i : ℝ) := by
      exact_mod_cast (by
        rw [Nat.mul_comm, Nat.choose_mul (Nat.le_add_right k i), Nat.add_sub_cancel_left] :
          (k + i).choose k * d.choose (k + i) = d.choose k * (d - k).choose i)
    have hsub : d - (k + i) = d - k - i := by omega
    simp only [f, basisNat, hsub, pow_add]
    linear_combination x ^ k * x ^ i * (1 - x) ^ (d - k - i) * hchoose
  change (∑ i ∈ Finset.range (d + 1), f i) = _
  rw [show d + 1 = k + (d - k + 1) by omega, Finset.sum_range_add, hzero, zero_add]
  simp only [hpoint, ← Finset.mul_sum, sum_basisNat, mul_one]

theorem bernstein_monomial {d k : ℕ} (hk : k ≤ d) (x : ℝ) :
    (∑ i : Fin (d + 1), ((i.val.choose k : ℝ) / (d.choose k : ℝ)) * bernsteinBasis d i x) = x ^ k := by
  change (∑ i : Fin (d + 1), ((i.val.choose k : ℝ) / (d.choose k : ℝ)) * basisNat d i.val x) = _
  rw [Fin.sum_univ_eq_sum_range (fun j => ((j.choose k : ℝ) / (d.choose k : ℝ)) * basisNat d j x)]
  simp only [div_eq_mul_inv]
  have h (i : ℕ) : (i.choose k : ℝ) * (d.choose k : ℝ)⁻¹ * basisNat d i x =
      ((i.choose k : ℝ) * basisNat d i x) * (d.choose k : ℝ)⁻¹ := by ring
  simp only [h, ← Finset.sum_mul, sum_choose_basisNat hk]
  apply (mul_inv_eq_iff_eq_mul₀ (by exact_mod_cast (Nat.choose_pos hk).ne')).2
  ring

theorem bernstein_weights_sum {n : ℕ} (degree : Fin n → ℕ) (x : Fin n → ℝ) :
    (∑ i : ∀ j, Fin (degree j + 1), ∏ j, bernsteinBasis (degree j) (i j) (x j)) = 1 := by
  rw [← Fintype.prod_sum (fun j (i : Fin (degree j + 1)) => bernsteinBasis (degree j) i (x j))]
  have h (j : Fin n) : (∑ i : Fin (degree j + 1), bernsteinBasis (degree j) i (x j)) = 1 := by
    change (∑ i : Fin (degree j + 1), basisNat (degree j) i.val (x j)) = _
    rw [Fin.sum_univ_eq_sum_range (fun i => basisNat (degree j) i (x j))]
    exact sum_basisNat _ _
  simp only [h, Finset.prod_const_one]

/-- The matrix entry for the monomial `x^k` in a degree-`d` Bernstein basis on `[lo, hi]`. -/
def axisPowerControl (R : Type*) [Field R] (d : ℕ) (lo hi : R) (i k : ℕ) : R :=
  ∑ b ∈ Finset.range (k + 1),
    ((i.choose b : R) / (d.choose b : R)) * (k.choose b : R) * lo ^ (k - b) * (hi - lo) ^ b

theorem axisPowerControl_sound {d k : ℕ} (hk : k ≤ d) (lo hi x : ℝ) :
    (∑ i : Fin (d + 1), axisPowerControl ℝ d lo hi i.val k * bernsteinBasis d i x) =
      (lo + (hi - lo) * x) ^ k := by
  unfold axisPowerControl
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  conv_rhs => rw [add_comm lo ((hi - lo) * x), add_pow]
  refine Finset.sum_congr rfl fun j hj => ?_
  have hj' : j ≤ d := (Nat.le_of_lt_succ (Finset.mem_range.mp hj)).trans hk
  calc
    _ = (∑ i : Fin (d + 1), ((i.val.choose j : ℝ) / (d.choose j : ℝ)) * bernsteinBasis d i x) *
        ((k.choose j : ℝ) * lo ^ (k - j) * (hi - lo) ^ j) := by
      rw [Finset.sum_mul]
      exact Finset.sum_congr rfl fun _ _ => by ring
    _ = _ := by rw [bernstein_monomial hj', mul_pow]; ring

/-- Exponents and a scalar coefficient of a sparse multivariate monomial. -/
abbrev PowerTerms (n : ℕ) (R : Type*) := List ((Fin n → ℕ) × R)

def powerEval {n : ℕ} {R : Type*} [CommSemiring R] (p : PowerTerms n R) (x : Fin n → R) : R :=
  (p.map fun (e, a) => a * ∏ j, x j ^ e j).sum

def powerControls {n : ℕ} {R : Type*} [Field R] (d : Fin n → ℕ)
    (lo hi : Fin n → R) (p : PowerTerms n R) (i : ∀ j, Fin (d j + 1)) : R :=
  (p.map fun (e, a) => a * ∏ j, axisPowerControl R (d j) (lo j) (hi j) (i j).val (e j)).sum

theorem powerControls_sound {n : ℕ} (d : Fin n → ℕ) (lo hi x : Fin n → ℝ) (p : PowerTerms n ℝ)
    (hp : ∀ t ∈ p, ∀ j, t.1 j ≤ d j) :
    tensorBernstein d (powerControls d lo hi p) x = powerEval p (fun j => lo j + (hi j - lo j) * x j) := by
  induction p with
  | nil => simp [tensorBernstein, powerControls, powerEval]
  | cons t p ih =>
    have ht := hp t (List.mem_cons_self)
    have hsplit (i : ∀ j, Fin (d j + 1)) : powerControls d lo hi (t :: p) i =
        t.2 * (∏ j, axisPowerControl ℝ (d j) (lo j) (hi j) (i j).val (t.1 j)) + powerControls d lo hi p i := rfl
    have hind := ih (fun u hu => hp u (List.mem_cons_of_mem _ hu))
    unfold tensorBernstein at hind ⊢
    simp only [hsplit, add_mul, Finset.sum_add_distrib]
    change _ = t.2 * (∏ j, (lo j + (hi j - lo j) * x j) ^ t.1 j) + powerEval p _
    rw [hind]
    congr 1
    let f (j : Fin n) (i : Fin (d j + 1)) :=
      axisPowerControl ℝ (d j) (lo j) (hi j) i.val (t.1 j) * bernsteinBasis (d j) i (x j)
    calc
      _ = t.2 * ∑ i : ∀ j, Fin (d j + 1), ∏ j, f j (i j) := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun _ _ => by simp only [f, Finset.prod_mul_distrib, mul_assoc]
      _ = _ := by
        rw [← Fintype.prod_sum f]
        congr 1
        exact Finset.prod_congr rfl fun j _ => axisPowerControl_sound (ht j) _ _ _

end Mordell
