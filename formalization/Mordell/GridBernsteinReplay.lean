import Mordell.BernsteinSplit
import Mordell.BernsteinLeaf

/-! Packed fixed-grid Bernstein replay. Packing is only storage; inequalities are checked on decoded lanes. -/

open scoped BigOperators

noncomputable section

namespace Mordell

abbrev GridIndex (d : Fin 3 → ℕ) := ∀ j, Fin (d j + 1)

def gridIndex (d : Fin 3 → ℕ) (i : GridIndex d) : ℕ :=
  ((i 0).val * (d 1 + 1) + (i 1).val) * (d 2 + 1) + (i 2).val

def gridLane (d : Fin 3 → ℕ) (packed : ℕ) (i : GridIndex d) : ℤ :=
  ((packed >>> (64 * gridIndex d i)) % 2 ^ (64 : ℕ) : ℕ) - (2 : ℤ) ^ (48 : ℕ)

def gridControls (d : Fin 3 → ℕ) (packed : ℕ) (i : GridIndex d) : ℚ :=
  (gridLane d packed i : ℚ) / 2 ^ (32 : ℕ)

def splitIndex {d : ℕ} (right : Bool) (k : Fin (d + 1)) : Fin (d + 1) := if right then k.rev else k

def splitInt (d : Fin 3 → ℕ) (axis : Fin 3) (right : Bool) (packed : ℕ) (i : GridIndex d) : ℤ :=
  (List.ofFn fun k : Fin (d axis + 1) =>
    ((splitIndex right (i axis)).val.choose (splitIndex right k).val : ℤ) * gridLane d packed (Function.update i axis k)).sum

def gridSplitOK (d : Fin 3 → ℕ) (axis : Fin 3) (right : Bool) (parent child : ℕ) : Prop :=
  ∀ i, 2 ^ (splitIndex right (i axis)).val * gridLane d child i ≤ splitInt d axis right parent i

instance (d : Fin 3 → ℕ) (axis : Fin 3) (right : Bool) (parent child : ℕ) :
    Decidable (gridSplitOK d axis right parent child) :=
  inferInstanceAs (Decidable (∀ i, (2 : ℤ) ^ (splitIndex right (i axis)).val * gridLane d child i ≤
    splitInt d axis right parent i))

theorem gridSplitOK_sound {d : Fin 3 → ℕ} {axis : Fin 3} {right : Bool} {parent child : ℕ}
    (h : gridSplitOK d axis right parent child) (i : GridIndex d) :
    (gridControls d child i : ℝ) ≤ splitControls d (fun j => (gridControls d parent j : ℝ)) axis right i := by
  let r := (splitIndex right (i axis)).val
  have hg : (0 : ℝ) < 2 ^ (32 : ℕ) := by positivity
  have hr : (0 : ℝ) < 2 ^ r := by positivity
  have hcast : (2 : ℝ) ^ r * (gridLane d child i : ℝ) ≤ (splitInt d axis right parent i : ℝ) := by
    exact_mod_cast h i
  simp only [gridControls, Rat.cast_div, Rat.cast_intCast, Rat.cast_pow, Rat.cast_ofNat]
  calc
    _ ≤ ((splitInt d axis right parent i : ℝ) / 2 ^ r) / 2 ^ (32 : ℕ) := by
      exact div_le_div_of_nonneg_right ((le_div_iff₀ hr).2 (by linarith only [hcast])) hg.le
    _ = _ := by
      simp only [splitInt, List.sum_ofFn, Int.cast_sum, Int.cast_mul, Int.cast_natCast, splitControls,
        midpointWeight, splitIndex, r, div_eq_mul_inv, Finset.sum_mul]
      exact Finset.sum_congr rfl fun _ _ => by ring

def boxHalf {n : ℕ} (b : RationalBox n) (axis : Fin n) (right : Bool) : RationalBox n :=
  if right then b.right axis else b.left axis

private theorem boxHalf_affine (b : RationalBox 3) (axis : Fin 3) (right : Bool) (u : Fin 3 → ℝ) :
    (fun j => ((boxHalf b axis right).lower j : ℝ) +
      (((boxHalf b axis right).upper j : ℝ) - (boxHalf b axis right).lower j) * u j) =
      fun j => (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) *
        Function.update u axis (midpointParameter right (u axis)) j := by
  funext j
  cases right <;> by_cases h : j = axis
  all_goals first
    | subst j
      simp [boxHalf, RationalBox.left, RationalBox.right, RationalBox.midpoint, midpointParameter]
      ring
    | simp [boxHalf, RationalBox.left, RationalBox.right, Function.update_of_ne h]

def GridEnvelope (d : Fin 3 → ℕ) (b : RationalBox 3) (packed : ℕ) (H : (Fin 3 → ℝ) → ℝ) : Prop :=
  ∀ u, (∀ j, 0 ≤ u j ∧ u j ≤ 1) → tensorBernstein d (fun i => (gridControls d packed i : ℝ)) u ≤
    H (fun j => (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u j)

theorem gridEnvelope_half {d : Fin 3 → ℕ} {b : RationalBox 3} {axis : Fin 3} {right : Bool}
    {parent child : ℕ} {H : (Fin 3 → ℝ) → ℝ}
    (h : GridEnvelope d b parent H) (hs : gridSplitOK d axis right parent child) :
    GridEnvelope d (boxHalf b axis right) child H := by
  intro u hu
  let v := Function.update u axis (midpointParameter right (u axis))
  have hv (j : Fin 3) : 0 ≤ v j ∧ v j ≤ 1 := by
    by_cases hj : j = axis
    · subst j
      have h := hu axis
      cases right <;> simp only [v, Function.update_self, midpointParameter, Bool.false_eq_true, if_false, if_true] <;>
        constructor <;> linarith [h.1, h.2]
    · simpa only [v, Function.update_of_ne hj] using hu j
  rw [boxHalf_affine]
  exact (splitControls_lower d _ _ axis right hu (gridSplitOK_sound hs)).trans (h v hv)

structure LogKnot where
  anchor : ℚ
  upper : ℚ

def gridLogOK (terms : ℕ) (d : Fin 3 → ℕ) (b : RationalBox 3) (packed : ℕ) (log : LogKnot) : Prop :=
  0 < log.anchor ∧ logUpperRat terms log.anchor ≤ log.upper ∧
    ∀ i, log.upper - 1 + gapDistance (affineKnot d b i) / log.anchor ≤ gridControls d packed i

instance (terms : ℕ) (d : Fin 3 → ℕ) (b : RationalBox 3) (packed : ℕ) (log : LogKnot) :
    Decidable (gridLogOK terms d b packed log) :=
  inferInstanceAs (Decidable (0 < log.anchor ∧ logUpperRat terms log.anchor ≤ log.upper ∧
    ∀ i, log.upper - 1 + gapDistance (affineKnot d b i) / log.anchor ≤ gridControls d packed i))

theorem gridLogOK_sound {terms : ℕ} {d : Fin 3 → ℕ} (hd : ∀ j, 1 ≤ d j) {b : RationalBox 3}
    {packed : ℕ} {log : LogKnot} (h : gridLogOK terms d b packed log) {H : (Fin 3 → ℝ) → ℝ}
    (henv : GridEnvelope d b packed H) {x : Fin 3 → ℝ} (hx : b.Contains x) (hx' : 0 < gapDistance x) :
    Real.log (gapDistance x) ≤ H x := by
  obtain ⟨ha, hu, hc⟩ := h
  apply log_tangent_controls_sound terms ha d hd b (gridControls d packed) H _ henv x hx hx'
  intro i
  change powerControls d b.lower b.upper (affineLogTerms log.anchor (logUpperRat terms log.anchor)) i ≤ _
  rw [affineLogControls_eq]
  linarith only [hc i, hu]

inductive GridTree where
  | pruned
  | tangent (log : LogKnot)
  | split (axis : Fin 3) (leftPacked : ℕ) (left : GridTree) (rightPacked : ℕ) (right : GridTree)

def GridTree.traceCode : GridTree → List Char
  | .pruned => ['P']
  | .tangent _ => ['A']
  | .split a _ left _ right => Char.ofNat ('0'.toNat + a.val) :: (left.traceCode ++ right.traceCode)

def GridTree.accept (terms : ℕ) (d : Fin 3 → ℕ) (minDistance : ℚ) (ordered : Bool)
    (b : RationalBox 3) (packed : ℕ) : GridTree → Bool
  | .pruned => gramPrune minDistance ordered b
  | .tangent log => decide (gridLogOK terms d b packed log)
  | .split a ml left mr right => decide (gridSplitOK d a false packed ml) &&
      left.accept terms d minDistance ordered (b.left a) ml &&
      decide (gridSplitOK d a true packed mr) && right.accept terms d minDistance ordered (b.right a) mr

theorem GridTree.accept_split {terms : ℕ} {d : Fin 3 → ℕ} {minDistance : ℚ} {ordered : Bool}
    {b : RationalBox 3} {packed ml mr : ℕ} {a : Fin 3} {left right : GridTree}
    (hl : gridSplitOK d a false packed ml) (hlt : left.accept terms d minDistance ordered (b.left a) ml = true)
    (hr : gridSplitOK d a true packed mr) (hrt : right.accept terms d minDistance ordered (b.right a) mr = true) :
    (GridTree.split a ml left mr right).accept terms d minDistance ordered b packed = true := by
  simp only [GridTree.accept, decide_eq_true_eq, Bool.and_eq_true]
  exact ⟨⟨⟨hl, hlt⟩, hr⟩, hrt⟩

theorem GridTree.sound {terms : ℕ} {d : Fin 3 → ℕ} (hd : ∀ j, 1 ≤ d j) {minDistance : ℚ}
    (hmin : 0 < minDistance) {ordered : Bool} {H : (Fin 3 → ℝ) → ℝ} (tree : GridTree)
    {b : RationalBox 3} {packed : ℕ} (h : tree.accept terms d minDistance ordered b packed = true)
    (henv : GridEnvelope d b packed H) :
    ∀ x, b.Contains x → gramPairDomain (minDistance : ℝ) ordered x → Real.log (gapDistance x) ≤ H x := by
  induction tree generalizing b packed with
  | pruned => exact fun _ hx hp => (gramPrune_sound h hx hp).elim
  | tangent log =>
    intro x hx hp
    have h₀ : (0 : ℝ) < (minDistance : ℝ) := by exact_mod_cast hmin
    exact gridLogOK_sound hd (of_decide_eq_true h) henv hx (h₀.trans_le hp.2.2.2.1)
  | split a ml left mr right lh rh =>
    simp only [accept, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨⟨⟨hl, hlt⟩, hr⟩, hrt⟩ := h
    intro x hx hp
    rcases RationalBox.contained_in_children hx a with hxl | hxr
    · exact lh hlt (gridEnvelope_half henv hl) x hxl hp
    · exact rh hrt (gridEnvelope_half henv hr) x hxr hp

end Mordell
