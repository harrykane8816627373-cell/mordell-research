import Mordell.Basic
import Mathlib.Tactic.GCongr

/-! Exact rational enclosures for scalar expressions evaluated over the reals. -/

namespace Mordell

structure RationalBall where
  mid : ℚ
  rad : ℚ

namespace RationalBall

def Contains (b : RationalBall) (x : ℝ) : Prop := |x - (b.mid : ℝ)| ≤ (b.rad : ℝ)

def point (x : ℚ) : RationalBall := ⟨x, 0⟩

def add (a b : RationalBall) : RationalBall := ⟨a.mid + b.mid, a.rad + b.rad⟩

def neg (a : RationalBall) : RationalBall := ⟨-a.mid, a.rad⟩

def sub (a b : RationalBall) : RationalBall := a.add b.neg

def mul (a b : RationalBall) : RationalBall :=
  ⟨a.mid * b.mid, |a.mid| * b.rad + |b.mid| * a.rad + a.rad * b.rad⟩

def inv (a : RationalBall) : RationalBall :=
  ⟨a.mid⁻¹, a.rad / (|a.mid| * (|a.mid| - a.rad))⟩

theorem point_sound (x : ℚ) : (point x).Contains (x : ℝ) := by simp [point, Contains]

theorem add_sound {a b : RationalBall} {x y : ℝ} (hx : a.Contains x) (hy : b.Contains y) :
    (a.add b).Contains (x + y) := by
  simp only [Contains, add, Rat.cast_add] at *
  rw [show x + y - ((a.mid : ℝ) + b.mid) = (x - a.mid) + (y - b.mid) by ring]
  exact (abs_add_le _ _).trans (add_le_add hx hy)

theorem neg_sound {a : RationalBall} {x : ℝ} (hx : a.Contains x) : a.neg.Contains (-x) := by
  simp only [Contains, neg, Rat.cast_neg]
  rw [show -x - -(a.mid : ℝ) = -(x - (a.mid : ℝ)) by ring, abs_neg]
  exact hx

theorem sub_sound {a b : RationalBall} {x y : ℝ} (hx : a.Contains x) (hy : b.Contains y) :
    (a.sub b).Contains (x - y) := by
  change (a.add b.neg).Contains (x - y)
  rw [sub_eq_add_neg]
  exact add_sound hx (neg_sound hy)

theorem mul_sound {a b : RationalBall} {x y : ℝ} (hx : a.Contains x) (hy : b.Contains y) :
    (a.mul b).Contains (x * y) := by
  change |x - (a.mid : ℝ)| ≤ a.rad at hx
  change |y - (b.mid : ℝ)| ≤ b.rad at hy
  let u := x - (a.mid : ℝ)
  let v := y - (b.mid : ℝ)
  simp only [Contains, mul, Rat.cast_mul, Rat.cast_add, Rat.cast_abs]
  rw [show x * y - (a.mid : ℝ) * b.mid = (a.mid : ℝ) * v + (b.mid : ℝ) * u + u * v by
    dsimp [u, v]; ring]
  calc
    _ ≤ |(a.mid : ℝ) * v| + |(b.mid : ℝ) * u| + |u * v| :=
      (abs_add_le _ _).trans (add_le_add (abs_add_le _ _) le_rfl)
    _ ≤ _ := by
      simp only [abs_mul]
      exact add_le_add (add_le_add (mul_le_mul_of_nonneg_left hy (abs_nonneg _))
        (mul_le_mul_of_nonneg_left hx (abs_nonneg _)))
        (mul_le_mul hx hy (abs_nonneg _) ((abs_nonneg _).trans hx))

theorem inv_sound {a : RationalBall} {x : ℝ} (hx : a.Contains x) (h : a.rad < |a.mid|) :
    a.inv.Contains x⁻¹ := by
  let m : ℝ := a.mid
  let r : ℝ := a.rad
  have hh : r < |m| := by change (a.rad : ℝ) < |(a.mid : ℝ)|; exact_mod_cast h
  have hx' : |x - m| ≤ r := hx
  have hr : 0 ≤ r := (abs_nonneg _).trans hx'
  have hm : 0 < |m| := hr.trans_lt hh
  have hl : |m| - r ≤ |x| := by
    have he := abs_add_le (m - x) x
    rw [sub_add_cancel, abs_sub_comm m x] at he
    linarith
  have hxn : x ≠ 0 := abs_pos.mp ((sub_pos.mpr hh).trans_le hl)
  have hmn : m ≠ 0 := abs_pos.mp hm
  simp only [Contains, inv, Rat.cast_inv, Rat.cast_div, Rat.cast_mul, Rat.cast_abs, Rat.cast_sub]
  change |x⁻¹ - m⁻¹| ≤ r / (|m| * (|m| - r))
  rw [show x⁻¹ - m⁻¹ = (m - x) / (x * m) by field_simp [hxn, hmn],
    abs_div, abs_mul, abs_sub_comm m x, mul_comm |x| |m|]
  exact div_le_div₀ hr hx' (mul_pos hm (sub_pos.mpr hh)) (mul_le_mul_of_nonneg_left hl hm.le)

end RationalBall
end Mordell
