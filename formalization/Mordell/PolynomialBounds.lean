import Mordell.Basic

/-! Sound coefficient bounds for a sparse polynomial on a unit box. -/

noncomputable section

namespace Mordell

abbrev Coefficients (n : ℕ) := List (ℚ × List (Fin n))

def monomialValue {n : ℕ} (x : Fin n → ℝ) (indices : List (Fin n)) : ℝ :=
  (indices.map x).prod

def evaluate {n : ℕ} (x : Fin n → ℝ) : Coefficients n → ℝ
  | [] => 0
  | (c, indices) :: ts => (c : ℝ) * monomialValue x indices + evaluate x ts

def magnitude {n : ℕ} : Coefficients n → ℚ
  | [] => 0
  | (c, _) :: ts => |c| + magnitude ts

def lipschitzBound {n : ℕ} : Coefficients n → ℚ
  | [] => 0
  | (c, indices) :: ts => |c| * indices.length + lipschitzBound ts

theorem abs_monomialValue_le_one {n : ℕ} {x : Fin n → ℝ}
    (hx : ∀ i, |x i| ≤ 1) (indices : List (Fin n)) :
    |monomialValue x indices| ≤ 1 := by
  induction indices with
  | nil => simp [monomialValue]
  | cons i indices ih =>
    simp only [monomialValue, List.map_cons, List.prod_cons, abs_mul]
    simpa using mul_le_mul (hx i) ih (abs_nonneg _) (by norm_num : 0 ≤ (1 : ℝ))

theorem abs_monomialValue_sub_le {n : ℕ} {x y : Fin n → ℝ} {d : ℝ}
    (hx : ∀ i, |x i| ≤ 1) (hy : ∀ i, |y i| ≤ 1) (hd : 0 ≤ d)
    (hxy : ∀ i, |x i - y i| ≤ d) (indices : List (Fin n)) :
    |monomialValue x indices - monomialValue y indices| ≤ indices.length * d := by
  induction indices with
  | nil => simp [monomialValue]
  | cons i indices ih =>
    let u := monomialValue x indices
    let v := monomialValue y indices
    change |x i * u - y i * v| ≤ ((indices.length + 1 : ℕ) : ℝ) * d
    have he : x i * u - y i * v = x i * (u - v) + (x i - y i) * v := by ring
    rw [he]
    calc
      |x i * (u - v) + (x i - y i) * v| ≤ |x i| * |u - v| + |x i - y i| * |v| := by
        simpa only [abs_mul] using abs_add_le (x i * (u - v)) ((x i - y i) * v)
      _ ≤ (indices.length : ℝ) * d + d * 1 := add_le_add
        (by simpa using mul_le_mul (hx i) ih (abs_nonneg _) (by norm_num : 0 ≤ (1 : ℝ)))
        (mul_le_mul (hxy i) (abs_monomialValue_le_one hy indices) (abs_nonneg _) hd)
      _ = _ := by push_cast; ring

theorem abs_evaluate_le {n : ℕ} {x : Fin n → ℝ}
    (hx : ∀ i, |x i| ≤ 1) (cs : Coefficients n) :
    |evaluate x cs| ≤ (magnitude cs : ℝ) := by
  induction cs with
  | nil => simp [evaluate, magnitude]
  | cons t cs ih =>
    obtain ⟨c, ix⟩ := t
    simp only [evaluate, magnitude, Rat.cast_add, Rat.cast_abs]
    calc
      |(c : ℝ) * monomialValue x ix + evaluate x cs| ≤
          |(c : ℝ) * monomialValue x ix| + |evaluate x cs| := abs_add_le _ _
      _ ≤ |(c : ℝ)| + (magnitude cs : ℝ) := by
        rw [abs_mul]
        have h := mul_le_mul_of_nonneg_left (abs_monomialValue_le_one hx ix) (abs_nonneg (c : ℝ))
        rw [mul_one] at h
        exact add_le_add h ih

theorem abs_evaluate_sub_le {n : ℕ} {x y : Fin n → ℝ} {d : ℝ}
    (hx : ∀ i, |x i| ≤ 1) (hy : ∀ i, |y i| ≤ 1) (hd : 0 ≤ d)
    (hxy : ∀ i, |x i - y i| ≤ d) (cs : Coefficients n) :
    |evaluate x cs - evaluate y cs| ≤ (lipschitzBound cs : ℝ) * d := by
  induction cs with
  | nil => simp [evaluate, lipschitzBound]
  | cons t cs ih =>
    obtain ⟨c, ix⟩ := t
    simp only [evaluate, lipschitzBound, Rat.cast_add, Rat.cast_mul, Rat.cast_abs,
      Rat.cast_natCast]
    have he :
        (c : ℝ) * monomialValue x ix + evaluate x cs -
          ((c : ℝ) * monomialValue y ix + evaluate y cs) =
          (c : ℝ) * (monomialValue x ix - monomialValue y ix) +
            (evaluate x cs - evaluate y cs) := by ring
    rw [he]
    calc
      |(c : ℝ) * (monomialValue x ix - monomialValue y ix) +
          (evaluate x cs - evaluate y cs)| ≤
          |(c : ℝ)| * |monomialValue x ix - monomialValue y ix| +
            |evaluate x cs - evaluate y cs| := by
              simpa only [abs_mul] using
                abs_add_le ((c : ℝ) * (monomialValue x ix - monomialValue y ix))
                  (evaluate x cs - evaluate y cs)
      _ ≤ |(c : ℝ)| * (ix.length * d) + (lipschitzBound cs : ℝ) * d :=
        add_le_add (mul_le_mul_of_nonneg_left
          (abs_monomialValue_sub_le hx hy hd hxy ix) (abs_nonneg _)) ih
      _ = _ := by ring

end Mordell
