import Mordell.GridBernsteinReplay

/-! Separable initialization of a three-dimensional Bernstein tensor from integer witnesses. -/

open scoped BigOperators

namespace Mordell

def tensorOf3 {d : Fin 3 → ℕ} {R : Type*}
    (c : Fin (d 0 + 1) → Fin (d 1 + 1) → Fin (d 2 + 1) → R) (i : GridIndex d) : R :=
  c (i 0) (i 1) (i 2)

private def gridOf3 {d : Fin 3 → ℕ}
    (a : Fin (d 0 + 1)) (b : Fin (d 1 + 1)) (c : Fin (d 2 + 1)) : GridIndex d :=
  Function.update (Function.update (Function.update (fun _ => 0) 0 a) 1 b) 2 c

@[simp] private theorem gridOf3_zero {d : Fin 3 → ℕ} (a : Fin (d 0 + 1)) (b : Fin (d 1 + 1))
    (c : Fin (d 2 + 1)) : gridOf3 a b c 0 = a := by simp [gridOf3]

@[simp] private theorem gridOf3_one {d : Fin 3 → ℕ} (a : Fin (d 0 + 1)) (b : Fin (d 1 + 1))
    (c : Fin (d 2 + 1)) : gridOf3 a b c 1 = b := by simp [gridOf3]

@[simp] private theorem gridOf3_two {d : Fin 3 → ℕ} (a : Fin (d 0 + 1)) (b : Fin (d 1 + 1))
    (c : Fin (d 2 + 1)) : gridOf3 a b c 2 = c := by simp [gridOf3]

private def gridEquiv {d : Fin 3 → ℕ} : Fin (d 0 + 1) × Fin (d 1 + 1) × Fin (d 2 + 1) ≃ GridIndex d where
  toFun p := gridOf3 p.1 p.2.1 p.2.2
  invFun i := (i 0, i 1, i 2)
  left_inv p := by simp
  right_inv i := by funext j; fin_cases j <;> simp

private theorem sum_gridIndex {d : Fin 3 → ℕ} {R : Type*} [AddCommMonoid R] (f : GridIndex d → R) :
    (∑ i, f i) = ∑ a, ∑ b, ∑ c, f (gridOf3 a b c) := by
  rw [← Equiv.sum_comp gridEquiv f]
  simp only [Fintype.sum_prod_type]
  rfl

/-- Apply a matrix along one coordinate, leaving the other indices fixed. -/
def tensorAxisMap {d : Fin 3 → ℕ} {R : Type*} [Semiring R] (axis : Fin 3)
    (m : Fin (d axis + 1) → Fin (d axis + 1) → R) (c : GridIndex d → R) (i : GridIndex d) : R :=
  (List.ofFn fun k => m (i axis) k * c (Function.update i axis k)).sum

theorem tensorAxisMap_three {d : Fin 3 → ℕ} {R : Type*} [CommSemiring R]
    (m : ∀ j, Fin (d j + 1) → Fin (d j + 1) → R) (c : GridIndex d → R) (i : GridIndex d) :
    tensorAxisMap 2 (m 2) (tensorAxisMap 1 (m 1) (tensorAxisMap 0 (m 0) c)) i =
      ∑ k, c k * ∏ j, m j (i j) (k j) := by
  rw [sum_gridIndex]
  simp only [tensorAxisMap, List.sum_ofFn, Function.update_of_ne (by decide : (1 : Fin 3) ≠ 2),
    Function.update_of_ne (by decide : (0 : Fin 3) ≠ 1), Function.update_of_ne (by decide : (0 : Fin 3) ≠ 2),
    Finset.mul_sum, Fin.prod_univ_three, gridOf3_zero, gridOf3_one, gridOf3_two]
  have h (a : Fin (d 0 + 1)) (b : Fin (d 1 + 1)) (c' : Fin (d 2 + 1)) :
      Function.update (Function.update (Function.update i 2 c') 1 b) 0 a = gridOf3 a b c' := by
    funext j
    fin_cases j <;> simp
  simp only [h]
  let f (a : Fin (d 0 + 1)) (b : Fin (d 1 + 1)) (k : Fin (d 2 + 1)) :=
    c (gridOf3 a b k) * (m 0 (i 0) a * m 1 (i 1) b * m 2 (i 2) k)
  calc
    _ = ∑ k, ∑ b, ∑ a, f a b k := by
      refine Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => ?_
      dsimp [f]
      ac_rfl
    _ = ∑ a, ∑ b, ∑ k, f a b k := by
      calc
        _ = ∑ k, ∑ a, ∑ b, f a b k := Finset.sum_congr rfl fun _ _ => Finset.sum_comm
        _ = ∑ a, ∑ k, ∑ b, f a b k := Finset.sum_comm
        _ = _ := Finset.sum_congr rfl fun _ _ => Finset.sum_comm

noncomputable def realTensor {d : Fin 3 → ℕ} (den : ℕ) (c : GridIndex d → ℤ) (i : GridIndex d) : ℝ :=
  (c i : ℝ) / den

theorem realTensor_axisMap {d : Fin 3 → ℕ} (axis : Fin 3) (s t : ℕ)
    (m : Fin (d axis + 1) → Fin (d axis + 1) → ℤ) (c : GridIndex d → ℤ) (i : GridIndex d) :
    realTensor (s * t) (tensorAxisMap axis m c) i =
      tensorAxisMap axis (fun j k => (m j k : ℝ) / t) (realTensor s c) i := by
  simp only [realTensor, tensorAxisMap, List.sum_ofFn, Int.cast_sum, Int.cast_mul, Nat.cast_mul,
    div_eq_mul_inv, Finset.sum_mul, mul_inv_rev]
  exact Finset.sum_congr rfl fun _ _ => by ring

def densePowerEval {d : Fin 3 → ℕ} {R : Type*} [CommSemiring R]
    (c : GridIndex d → R) (x : Fin 3 → R) : R := ∑ k, c k * ∏ j, x j ^ (k j).val

noncomputable def rootPowerControls (d : Fin 3 → ℕ) (b : RationalBox 3) (c : GridIndex d → ℝ) (i : GridIndex d) : ℝ :=
  ∑ k, c k * ∏ j, axisPowerControl ℝ (d j) (b.lower j) (b.upper j) (i j).val (k j).val

theorem rootPowerControls_sound (d : Fin 3 → ℕ) (b : RationalBox 3) (c : GridIndex d → ℝ) (u : Fin 3 → ℝ) :
    tensorBernstein d (rootPowerControls d b c) u =
      densePowerEval c (fun j => (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u j) := by
  unfold tensorBernstein rootPowerControls densePowerEval
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  let f (j : Fin 3) (i : Fin (d j + 1)) :=
    axisPowerControl ℝ (d j) (b.lower j) (b.upper j) i.val (k j).val * bernsteinBasis (d j) i (u j)
  calc
    _ = c k * ∑ i : GridIndex d, ∏ j, f j (i j) := by
      simp only [f, Finset.mul_sum, Finset.prod_mul_distrib, mul_assoc]
    _ = _ := by
      rw [← Fintype.prod_sum f]
      congr 1
      exact Finset.prod_congr rfl fun j _ => axisPowerControl_sound (Nat.le_of_lt_succ (k j).isLt) _ _ _

theorem gridEnvelope_of_densePower {d : Fin 3 → ℕ} {b : RationalBox 3} {packed : ℕ} {c : GridIndex d → ℝ}
    (h : ∀ i, (gridControls d packed i : ℝ) ≤ rootPowerControls d b c i) :
    GridEnvelope d b packed (densePowerEval c) := by
  intro u hu
  rw [← rootPowerControls_sound d b c]
  exact tensorBernstein_mono hu h

/-- Replace every discarded negative monomial by its value at the upper corner. -/
def discardedLower {n : ℕ} (cap : Fin n → ℚ) (p : PowerTerms n ℚ) : ℚ :=
  powerEval (p.map fun (e, a) => (e, min a 0)) cap

theorem discardedLower_sound {n : ℕ} (cap : Fin n → ℚ) {x : Fin n → ℝ}
    (h : ∀ j, 0 ≤ x j ∧ x j ≤ (cap j : ℝ)) (p : PowerTerms n ℚ) :
    (discardedLower cap p : ℝ) ≤ powerEval (castPowerTerms p) x := by
  induction p with
  | nil => simp [discardedLower, castPowerTerms, powerEval]
  | cons t p ih =>
    have hMon : (∏ j, x j ^ t.1 j) ≤ ∏ j, (cap j : ℝ) ^ t.1 j := by
      exact Finset.prod_le_prod (fun j _ => pow_nonneg (h j).1 _) (fun j _ => pow_le_pow_left₀ (h j).1 (h j).2 _)
    have he : ((discardedLower cap (t :: p)) : ℝ) =
        (min t.2 0 : ℚ) * (∏ j, (cap j : ℝ) ^ t.1 j) + (discardedLower cap p : ℝ) := by
      simp only [discardedLower, powerEval, List.map_cons, List.sum_cons, Rat.cast_add, Rat.cast_mul,
        Rat.cast_prod, Rat.cast_pow]
    have ht : ((min t.2 0 : ℚ) : ℝ) * (∏ j, (cap j : ℝ) ^ t.1 j) ≤ (t.2 : ℝ) * ∏ j, x j ^ t.1 j := by
      by_cases ht : t.2 ≤ 0
      · rw [min_eq_left ht]
        exact mul_le_mul_of_nonpos_left hMon (by exact_mod_cast ht)
      · rw [min_eq_right (le_of_lt (lt_of_not_ge ht)), Rat.cast_zero, zero_mul]
        exact mul_nonneg (by exact_mod_cast le_of_lt (lt_of_not_ge ht))
          (Finset.prod_nonneg fun j _ => pow_nonneg (h j).1 _)
    rw [he]
    exact add_le_add ht ih

theorem powerEval_mono_same_exponents {n : ℕ} (p : PowerTerms n ℚ) (a : (Fin n → ℕ) × ℚ → ℝ)
    {x : Fin n → ℝ} (hx : ∀ j, 0 ≤ x j) (ha : ∀ t ∈ p, (t.2 : ℝ) ≤ a t) :
    powerEval (castPowerTerms p) x ≤ powerEval (p.map fun t => (t.1, a t)) x := by
  induction p with
  | nil => simp [powerEval, castPowerTerms]
  | cons t p ih =>
    exact add_le_add (mul_le_mul_of_nonneg_right (ha t List.mem_cons_self)
      (Finset.prod_nonneg fun j _ => pow_nonneg (hx j) _)) (ih (fun u hu => ha u (List.mem_cons_of_mem _ hu)))

end Mordell
