import Mordell.Cover

/-! Replay and soundness for the prefix-order box transcripts. -/

namespace Mordell

structure RationalBox (n : ℕ) where
  lower : Fin n → ℚ
  upper : Fin n → ℚ

namespace RationalBox

def toBox {n : ℕ} (b : RationalBox n) : Box n := ⟨fun j => b.lower j, fun j => b.upper j⟩

abbrev Contains {n : ℕ} (b : RationalBox n) (x : Fin n → ℝ) : Prop := b.toBox.Contains x

def midpoint {n : ℕ} (b : RationalBox n) (i : Fin n) : ℚ := (b.lower i + b.upper i) / 2

def left {n : ℕ} (b : RationalBox n) (i : Fin n) : RationalBox n :=
  ⟨b.lower, Function.update b.upper i (b.midpoint i)⟩

def right {n : ℕ} (b : RationalBox n) (i : Fin n) : RationalBox n :=
  ⟨Function.update b.lower i (b.midpoint i), b.upper⟩

theorem left_toBox {n : ℕ} (b : RationalBox n) (i : Fin n) : (b.left i).toBox = b.toBox.left i := by
  unfold toBox left Box.left
  apply congrArg (Box.mk _)
  funext j
  by_cases h : j = i
  · subst j
    simp only [Function.update_self, midpoint, Box.midpoint, Rat.cast_div, Rat.cast_ofNat, Rat.cast_add]
  · simp only [Function.update_of_ne h]

theorem right_toBox {n : ℕ} (b : RationalBox n) (i : Fin n) : (b.right i).toBox = b.toBox.right i := by
  unfold toBox right Box.right
  apply congrArg (fun lo => Box.mk lo _)
  funext j
  by_cases h : j = i
  · subst j
    simp only [Function.update_self, midpoint, Box.midpoint, Rat.cast_div, Rat.cast_ofNat, Rat.cast_add]
  · simp only [Function.update_of_ne h]

theorem contained_in_children {n : ℕ} {b : RationalBox n} {x : Fin n → ℝ} (hx : b.Contains x) (i : Fin n) :
    (b.left i).Contains x ∨ (b.right i).Contains x := by
  simpa only [Contains, left_toBox, right_toBox] using Box.contained_in_children hx i

theorem affineUnit_coordinates {n : ℕ} {b : RationalBox n} {x : Fin n → ℝ} (hx : b.Contains x) :
    ∃ u : Fin n → ℝ, (∀ j, 0 ≤ u j ∧ u j ≤ 1) ∧
      ∀ j, x j = (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u j := by
  have h (j : Fin n) : ∃ u : ℝ, (0 ≤ u ∧ u ≤ 1) ∧
      x j = (b.lower j : ℝ) + ((b.upper j : ℝ) - b.lower j) * u := by
    have ⟨hlo, hhi⟩ := hx j
    change (b.lower j : ℝ) ≤ x j at hlo
    change x j ≤ (b.upper j : ℝ) at hhi
    by_cases he : (b.lower j : ℝ) = (b.upper j : ℝ)
    · exact ⟨0, by norm_num, by linarith⟩
    · have hp := sub_pos.mpr (lt_of_le_of_ne (hlo.trans hhi) he)
      refine ⟨(x j - b.lower j) / ((b.upper j : ℝ) - b.lower j), ?_, ?_⟩
      · exact ⟨div_nonneg (sub_nonneg.mpr hlo) hp.le, (div_le_one hp).mpr (by linarith)⟩
      · field_simp [hp.ne']
        ring
  choose u hu using h
  exact ⟨u, fun j => (hu j).1, fun j => (hu j).2⟩

def includes {n : ℕ} (outer inner : RationalBox n) : Bool :=
  (List.ofFn fun i => decide (outer.lower i ≤ inner.lower i ∧ inner.upper i ≤ outer.upper i)).all id

theorem includes_sound {n : ℕ} {outer inner : RationalBox n} (h : outer.includes inner = true)
    {x : Fin n → ℝ} (hx : inner.Contains x) : outer.Contains x := by
  intro j
  have hj := (List.all_eq_true.mp h) _ (List.mem_ofFn.mpr ⟨j, rfl⟩)
  obtain ⟨hlo, hhi⟩ := of_decide_eq_true (show decide
    (outer.lower j ≤ inner.lower j ∧ inner.upper j ≤ outer.upper j) = true from hj)
  exact ⟨(by exact_mod_cast hlo : (outer.lower j : ℝ) ≤ inner.lower j).trans (hx j).1,
    (hx j).2.trans (by exact_mod_cast hhi : (inner.upper j : ℝ) ≤ outer.upper j)⟩

end RationalBox

inductive ReplayAction (n : ℕ) (tag : Type) where
  | split : Fin n → ReplayAction n tag
  | leaf : tag → ReplayAction n tag

/-- The left child is visited first. A successful replay exhausts the stack. -/
def replayBoxes {n : ℕ} {tag : Type} (check : RationalBox n → tag → Bool) :
    List (ReplayAction n tag) → List (RationalBox n) → Bool
  | [], stack => stack.isEmpty
  | _ :: _, [] => false
  | .split i :: actions, box :: stack => replayBoxes check actions (box.left i :: box.right i :: stack)
  | .leaf label :: actions, box :: stack => check box label && replayBoxes check actions stack

theorem replayBoxes_sound {n : ℕ} {tag : Type} {domain goal : (Fin n → ℝ) → Prop}
    (check : RationalBox n → tag → Bool)
    (sound : ∀ b label, check b label = true → ∀ x, b.Contains x → domain x → goal x)
    {actions : List (ReplayAction n tag)} {stack : List (RationalBox n)}
    (ok : replayBoxes check actions stack = true) :
    ∀ b ∈ stack, ∀ x, b.Contains x → domain x → goal x := by
  induction actions generalizing stack with
  | nil =>
    have he : stack = [] := List.isEmpty_iff.mp ok
    simp [he]
  | cons action actions ih =>
    cases stack with
    | nil => simp
    | cons box stack =>
      cases action with
      | split i =>
        have hsub := ih (stack := box.left i :: box.right i :: stack) ok
        intro b hb x hx hd
        rcases List.mem_cons.mp hb with rfl | htail
        · rcases RationalBox.contained_in_children hx i with hleft | hright
          · exact hsub _ (by simp) x hleft hd
          · exact hsub _ (by simp) x hright hd
        · exact hsub _ (by simp [htail]) x hx hd
      | leaf label =>
        have hok : check box label = true ∧ replayBoxes check actions stack = true := by
          simpa only [replayBoxes, Bool.and_eq_true] using ok
        intro b hb x hx hd
        rcases List.mem_cons.mp hb with rfl | htail
        · exact sound _ _ hok.1 x hx hd
        · exact ih hok.2 _ htail x hx hd

def gapDistance {R : Type*} [Ring R] (x : Fin 3 → R) : R := x 0 + x 1 - 2 * x 2

def gramPairDomain (minDistance : ℝ) (ordered : Bool) (x : Fin 3 → ℝ) : Prop :=
  0 ≤ x 0 ∧ 0 ≤ x 1 ∧ x 2 ^ 2 ≤ x 0 * x 1 ∧ minDistance ≤ gapDistance x ∧
    (ordered = true → x 0 ≤ x 1)

def minimumAbs (b : RationalBox 3) : ℚ := max 0 (max (b.lower 2) (-(b.upper 2)))

def maxDistance (b : RationalBox 3) : ℚ := b.upper 0 + b.upper 1 - 2 * b.lower 2

def gramPrune (minDistance : ℚ) (ordered : Bool) (b : RationalBox 3) : Bool :=
  (ordered && decide (b.upper 1 < b.lower 0)) ||
    decide (b.upper 0 * b.upper 1 < minimumAbs b ^ 2) || decide (maxDistance b < minDistance)

theorem gramPrune_sound {minDistance : ℚ} {ordered : Bool} {b : RationalBox 3}
    (h : gramPrune minDistance ordered b = true) {x : Fin 3 → ℝ} (hx : b.Contains x) :
    ¬ gramPairDomain (minDistance : ℝ) ordered x := by
  rintro ⟨hu, hv, hg, hd, ho⟩
  simp only [gramPrune, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq] at h
  rcases h with (⟨hord, hsep⟩ | hgram) | hdist
  · have hs : (b.upper 1 : ℝ) < b.lower 0 := by exact_mod_cast hsep
    exact hs.not_ge ((hx 0).1.trans ((ho hord).trans (hx 1).2))
  · have hmin : (minimumAbs b : ℝ) ≤ |x 2| := by
      simp only [minimumAbs, Rat.cast_max, Rat.cast_zero, Rat.cast_neg, max_le_iff]
      exact ⟨abs_nonneg _, (hx 2).1.trans (le_abs_self _), (neg_le_neg (hx 2).2).trans (neg_le_abs _)⟩
    have hsq : (minimumAbs b : ℝ) ^ 2 ≤ x 2 ^ 2 := by
      rw [← sq_abs (x 2)]
      exact pow_le_pow_left₀ (by simp [minimumAbs]) hmin 2
    have hop : x 0 * x 1 ≤ (b.upper 0 : ℝ) * b.upper 1 :=
      mul_le_mul (hx 0).2 (hx 1).2 hv (hu.trans (hx 0).2)
    have hsep : (b.upper 0 : ℝ) * b.upper 1 < (minimumAbs b : ℝ) ^ 2 := by exact_mod_cast hgram
    exact hsep.not_ge (hsq.trans (hg.trans hop))
  · have hsep : (maxDistance b : ℝ) < minDistance := (Rat.cast_lt (K := ℝ)).2 hdist
    have h' : gapDistance x ≤ (maxDistance b : ℝ) := by
      simp only [gapDistance, maxDistance, Rat.cast_add, Rat.cast_sub, Rat.cast_mul, Rat.cast_ofNat]
      have h0 := (hx 0).2
      have h1 := (hx 1).2
      have ht := (hx 2).1
      change x 0 ≤ (b.upper 0 : ℝ) at h0
      change x 1 ≤ (b.upper 1 : ℝ) at h1
      change (b.lower 2 : ℝ) ≤ x 2 at ht
      linarith
    exact hsep.not_ge (hd.trans h')

end Mordell
