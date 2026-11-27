import Mordell.Basic

/-! The algebraic specification of the ten-point candidate at mass nine. -/

open scoped BigOperators

noncomputable section

namespace Mordell

structure Parameters where
  a : ℝ
  b : ℝ
  c : ℝ
  x : ℝ
  y : ℝ

namespace Parameters

def InRootBox (p : Parameters) : Prop :=
  p.a ∈ Set.Icc (118036 / 1048576 : ℝ) (118038 / 1048576) ∧
  p.b ∈ Set.Icc (1084059 / 1048576 : ℝ) (1084061 / 1048576) ∧
  p.c ∈ Set.Icc (1216216 / 1048576 : ℝ) (1216218 / 1048576) ∧
  p.x ∈ Set.Icc (-797428 / 1048576 : ℝ) (-797426 / 1048576) ∧
  p.y ∈ Set.Icc (810045 / 1048576 : ℝ) (810047 / 1048576)

structure Admissible (p : Parameters) : Prop where
  a_pos : 0 < p.a
  a_lt_b : p.a < p.b
  b_lt_c : p.b < p.c
  x_lower : -p.b < p.x
  x_upper : p.x < p.b
  y_lower : -p.c < p.y
  y_upper : p.y < p.c
  x_ne : p.x ≠ 0
  y_ne : p.y ≠ 0

theorem InRootBox.admissible {p : Parameters} (h : p.InRootBox) : p.Admissible := by
  rcases h with ⟨⟨ha₁, ha₂⟩, ⟨hb₁, hb₂⟩, ⟨hc₁, hc₂⟩, ⟨hx₁, hx₂⟩, ⟨hy₁, hy₂⟩⟩
  constructor <;> norm_num at * <;> linarith

def quadPower (r t : ℝ) : ℕ → ℝ
  | 0 => 2
  | 1 => 2 * t
  | k + 2 => 2 * t * quadPower r t (k + 1) - r ^ 2 * quadPower r t k

def powerSum (p : Parameters) (k : ℕ) : ℝ :=
  p.a ^ k + quadPower p.b p.x k + quadPower p.c p.y k

def weightedPowerSum (p : Parameters) (k : ℕ) : ℝ :=
  p.a * p.a ^ k + p.b * quadPower p.b p.x k + p.c * quadPower p.c p.y k

def equation (p : Parameters) (k : ℕ) : ℝ :=
  2 * (∑ j ∈ Finset.range (k + 1), p.powerSum j * p.powerSum (k - j)) -
    (2 * (k : ℝ) + 1) * p.powerSum k - 10 * p.weightedPowerSum k

def IsRoot (p : Parameters) : Prop :=
  p.InRootBox ∧ ∀ k : Fin 5, p.equation k = 0

theorem equation_zero (p : Parameters) :
    p.equation 0 = 45 - 10 * (p.a + 2 * p.b + 2 * p.c) := by
  simp [equation, powerSum, weightedPowerSum, quadPower]
  ring

theorem IsRoot.weighted_mass {p : Parameters} (h : p.IsRoot) :
    2 * p.a + 4 * p.b + 4 * p.c = 9 := by
  have he : p.equation 0 = 0 := h.2 0
  rw [equation_zero] at he
  linarith

end Parameters

def rectanglePoint (r t : ℝ) : ℂ :=
  ⟨Real.sqrt ((r + t) / 2), Real.sqrt ((r - t) / 2)⟩

theorem normSq_rectanglePoint {r t : ℝ} (h₁ : -r ≤ t) (h₂ : t ≤ r) :
    Complex.normSq (rectanglePoint r t) = r := by
  simp only [Complex.normSq_apply, rectanglePoint]
  rw [← sq, ← sq, Real.sq_sqrt (by linarith : 0 ≤ (r + t) / 2),
    Real.sq_sqrt (by linarith : 0 ≤ (r - t) / 2)]
  ring

def candidate (p : Parameters) : Configuration 10 :=
  let za : ℂ := (Real.sqrt p.a : ℝ)
  let zb := rectanglePoint p.b p.x
  let zc := rectanglePoint p.c p.y
  ![za, -za, zb, -zb, star zb, -star zb, zc, -zc, star zc, -star zc]

theorem mass_candidate (p : Parameters) (hp : p.Admissible) :
    mass (candidate p) = 2 * p.a + 4 * p.b + 4 * p.c := by
  have ha : Complex.normSq ((Real.sqrt p.a : ℝ) : ℂ) = p.a := by
    rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt hp.a_pos.le]
  have hb := normSq_rectanglePoint hp.x_lower.le hp.x_upper.le
  have hc := normSq_rectanglePoint hp.y_lower.le hp.y_upper.le
  unfold mass candidate
  simp only [Fin.sum_univ_succ, Matrix.cons_val_zero, Matrix.cons_val_succ, Finset.univ_eq_empty,
    Finset.sum_empty, Complex.normSq_neg, Complex.star_def, Complex.normSq_conj, ha, hb, hc]
  ring

theorem mass_candidate_of_root {p : Parameters} (hp : p.IsRoot) :
    mass (candidate p) = 9 :=
  (mass_candidate p hp.1.admissible).trans hp.weighted_mass

set_option maxHeartbeats 1500000 in
theorem candidate_injective {p : Parameters} (hp : p.Admissible) :
    Function.Injective (candidate p) := by
  have ha : Complex.normSq ((Real.sqrt p.a : ℝ) : ℂ) = p.a := by
    rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt hp.a_pos.le]
  have hb := normSq_rectanglePoint hp.x_lower.le hp.x_upper.le
  have hc := normSq_rectanglePoint hp.y_lower.le hp.y_upper.le
  have hsa := Real.sqrt_pos.2 hp.a_pos
  have hsbr := Real.sqrt_pos.2 (show 0 < (p.b + p.x) / 2 by linarith [hp.x_lower])
  have hsbi := Real.sqrt_pos.2 (show 0 < (p.b - p.x) / 2 by linarith [hp.x_upper])
  have hscr := Real.sqrt_pos.2 (show 0 < (p.c + p.y) / 2 by linarith [hp.y_lower])
  have hsci := Real.sqrt_pos.2 (show 0 < (p.c - p.y) / 2 by linarith [hp.y_upper])
  intro i j h
  have hn := congrArg Complex.normSq h
  have hr := congrArg Complex.re h
  have hi := congrArg Complex.im h
  fin_cases i <;> fin_cases j
  all_goals first | rfl | skip
  all_goals
    simp [candidate, ha, hb, hc] at hn
    try { exfalso; linarith [hp.a_lt_b, hp.b_lt_c] }
  all_goals
    simp [candidate, rectanglePoint, -Real.sqrt_div, -Real.sqrt_div'] at hr hi
    exfalso
    linarith

private def populationCenter (a b c : ℕ) : ℕ :=
  118037 * a + 1084060 * b + 1216217 * c

set_option maxRecDepth 10000 in
set_option maxHeartbeats 0 in
private theorem checked_population :
    ∀ a b c : Fin 11, a.val + b.val + c.val = 10 →
      9 * 1048576 - 10 ≤ populationCenter a.val b.val c.val →
      populationCenter a.val b.val c.val ≤ 9 * 1048576 + 10 →
      a.val = 2 ∧ b.val = 4 ∧ c.val = 4 := by
  decide

/-- The isolating box and the exact mass constraint lock the populations. -/
theorem populations_of_mass {p : Parameters} (hp : p.InRootBox)
    {a b c : ℕ} (htotal : a + b + c = 10)
    (hmass : (a : ℝ) * p.a + (b : ℝ) * p.b + (c : ℝ) * p.c = 9) :
    a = 2 ∧ b = 4 ∧ c = 4 := by
  rcases hp with ⟨⟨ha₁, ha₂⟩, ⟨hb₁, hb₂⟩, ⟨hc₁, hc₂⟩, _⟩
  have ht : (a : ℝ) + (b : ℝ) + (c : ℝ) = 10 := by exact_mod_cast htotal
  have lower :
      (9 * 1048576 - 10 : ℝ) ≤ (a : ℝ) * 118037 + (b : ℝ) * 1084060 + (c : ℝ) * 1216217 := by
    have ha := mul_le_mul_of_nonneg_left ha₂ (Nat.cast_nonneg (α := ℝ) a)
    have hb := mul_le_mul_of_nonneg_left hb₂ (Nat.cast_nonneg (α := ℝ) b)
    have hc := mul_le_mul_of_nonneg_left hc₂ (Nat.cast_nonneg (α := ℝ) c)
    linarith
  have upper :
      (a : ℝ) * 118037 + (b : ℝ) * 1084060 + (c : ℝ) * 1216217 ≤ (9 * 1048576 + 10 : ℝ) := by
    have ha := mul_le_mul_of_nonneg_left ha₁ (Nat.cast_nonneg (α := ℝ) a)
    have hb := mul_le_mul_of_nonneg_left hb₁ (Nat.cast_nonneg (α := ℝ) b)
    have hc := mul_le_mul_of_nonneg_left hc₁ (Nat.cast_nonneg (α := ℝ) c)
    linarith
  have hlo : 9 * 1048576 - 10 ≤ populationCenter a b c := by
    have hh : (9437174 : ℝ) ≤ (118037 * a + 1084060 * b + 1216217 * c : ℕ) := by
      push_cast
      linarith
    change 9437174 ≤ 118037 * a + 1084060 * b + 1216217 * c
    exact_mod_cast hh
  have hhi : populationCenter a b c ≤ 9 * 1048576 + 10 := by
    have hh : ((118037 * a + 1084060 * b + 1216217 * c : ℕ) : ℝ) ≤ 9437194 := by
      push_cast
      linarith
    change 118037 * a + 1084060 * b + 1216217 * c ≤ 9437194
    exact_mod_cast hh
  exact checked_population ⟨a, by omega⟩ ⟨b, by omega⟩ ⟨c, by omega⟩ htotal hlo hhi

end Mordell
