import Mordell.TwelvePoint
import Mordell.Rigidity

/-! Recovery of a triangle or a pair of cubic images from the equality moments. -/

noncomputable section

namespace Mordell

theorem cube_roots_same_orbit {z w : ℂ} (h : z ^ 3 = w ^ 3) :
    z = w ∨ z = omega * w ∨ z = omega ^ 2 * w := by
  have hsum : omega + omega ^ 2 = -1 := by linear_combination omega_triangle
  have hprod : omega * omega ^ 2 = 1 := by
    rw [← pow_succ', omega_cube]
  have hfactor : (z - omega * w) * (z - omega ^ 2 * w) = z ^ 2 + w * z + w ^ 2 := by
    calc
      _ = z ^ 2 - (omega + omega ^ 2) * (w * z) + (omega * omega ^ 2) * w ^ 2 := by ring
      _ = _ := by rw [hsum, hprod]; ring
  have hz : (z - w) * ((z - omega * w) * (z - omega ^ 2 * w)) = 0 := by
    rw [hfactor]
    linear_combination h
  rcases mul_eq_zero.mp hz with h₀ | h₁
  · exact Or.inl (sub_eq_zero.mp h₀)
  · exact Or.inr ((mul_eq_zero.mp h₁).imp sub_eq_zero.mp sub_eq_zero.mp)

theorem three_point_layer_is_triangle {a b c : ℂ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (h1 : a + b + c = 0) (h2 : a ^ 2 + b ^ 2 + c ^ 2 = 0) :
    (b = omega * a ∧ c = omega ^ 2 * a) ∨ (b = omega ^ 2 * a ∧ c = omega * a) := by
  obtain ⟨ha, hb, hc⟩ := three_points_cubes_eq h1 h2
  have hb' := (cube_roots_same_orbit (hb.trans ha.symm)).resolve_left hab.symm
  have hc' := (cube_roots_same_orbit (hc.trans ha.symm)).resolve_left hac.symm
  rcases hb' with hb' | hb' <;> rcases hc' with hc' | hc'
  · exact (hbc (hb'.trans hc'.symm)).elim
  · exact Or.inl ⟨hb', hc'⟩
  · exact Or.inr ⟨hb', hc'⟩
  · exact (hbc (hb'.trans hc'.symm)).elim

theorem cubic_layer_images_eq {r : ℝ} {a b c d : ℂ}
    (ha : Complex.normSq a = r) (hb : Complex.normSq b = r)
    (hc : Complex.normSq c = r) (hd : Complex.normSq d = r)
    (hs : a ^ 3 + b ^ 3 = c ^ 3 + d ^ 3) (hne : a ^ 3 + b ^ 3 ≠ 0) :
    (a ^ 3 = c ^ 3 ∧ b ^ 3 = d ^ 3) ∨ (a ^ 3 = d ^ 3 ∧ b ^ 3 = c ^ 3) := by
  have hp (x : ℂ) (h : Complex.normSq x = r) : Complex.normSq (x ^ 3) = r ^ 3 := by
    rw [map_pow, h]
  exact same_radius_pair_eq (hp a ha) (hp b hb) (hp c hc) (hp d hd) hs hne

end Mordell
