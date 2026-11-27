import Mordell.TriangleRigidity
import Mathlib.RingTheory.MvPolynomial.Symmetric.NewtonIdentities
import Mathlib.RingTheory.Polynomial.Vieta

/-! Newton identities and cubic fibers for the six-point outer layer. -/

open scoped BigOperators

noncomputable section

namespace Mordell

private def layerElements {n : ℕ} (z : Configuration n) : Multiset ℂ :=
  (Finset.univ : Finset (Fin n)).val.map z

private def elementary {n : ℕ} (z : Configuration n) (k : ℕ) : ℂ :=
  (layerElements z).esymm k

private theorem scalar_newton {n : ℕ} (z : Configuration n) {k : ℕ} (hk : 0 < k) :
    (∑ i, z i ^ k) = (-1 : ℂ) ^ (k + 1) * (k : ℂ) * elementary z k -
      ∑ j ∈ Finset.range (k + 1),
        if 0 < j ∧ j < k then (-1 : ℂ) ^ j * elementary z j * (∑ i, z i ^ (k - j)) else 0 := by
  have h := congrArg (MvPolynomial.aeval z : MvPolynomial (Fin n) ℂ →ₐ[ℂ] ℂ)
    (MvPolynomial.psum_eq_mul_esymm_sub_sum (Fin n) ℂ k hk)
  simp only [map_sub, map_mul, map_neg, map_one, map_pow, map_natCast, map_sum,
    MvPolynomial.psum, MvPolynomial.aeval_X, MvPolynomial.aeval_esymm_eq_multiset_esymm] at h
  rw [h, Finset.sum_filter]
  exact congrArg (_ - ·) (Finset.Nat.sum_antidiagonal_eq_sum_range_succ
    (fun a b => if 0 < a ∧ a < k then (-1 : ℂ) ^ a * elementary z a *
      (∑ i, z i ^ b) else 0) k)

private theorem elementary_zero {n : ℕ} (z : Configuration n) : elementary z 0 = 1 := by
  simpa only [MvPolynomial.aeval_esymm_eq_multiset_esymm, elementary, layerElements, map_one]
    using congrArg (MvPolynomial.aeval z : MvPolynomial (Fin n) ℂ →ₐ[ℂ] ℂ)
      (MvPolynomial.esymm_zero (Fin n) ℂ)

private theorem six_roots_cube_quadratic (z : Configuration 6)
    (h1 : ∑ i, z i ^ 1 = 0) (h2 : ∑ i, z i ^ 2 = 0)
    (h4 : ∑ i, z i ^ 4 = 0) (h5 : ∑ i, z i ^ 5 = 0) (i : Fin 6) :
    (z i ^ 3) ^ 2 - elementary z 3 * z i ^ 3 + elementary z 6 = 0 := by
  have hN1 := scalar_newton z (k := 1) (by decide)
  have hN2 := scalar_newton z (k := 2) (by decide)
  have hN4 := scalar_newton z (k := 4) (by decide)
  have hN5 := scalar_newton z (k := 5) (by decide)
  have h1' : ∑ i, z i = 0 := by simpa only [pow_one] using h1
  have hE1 : elementary z 1 = 0 := by
    norm_num [Finset.sum_range_succ, h1'] at hN1
    exact hN1.symm
  have hE2 : elementary z 2 = 0 := by
    norm_num [Finset.sum_range_succ, hE1] at hN2
    linear_combination (hN2 - h2) / 2
  have hE4 : elementary z 4 = 0 := by
    norm_num [Finset.sum_range_succ, hE1, hE2, h1'] at hN4
    linear_combination (hN4 - h4) / 4
  have hE5 : elementary z 5 = 0 := by
    norm_num [Finset.sum_range_succ, hE1, hE2, hE4, h2] at hN5
    linear_combination (h5 - hN5) / 5
  let s := layerElements z
  have hcard : s.card = 6 := by simp [s, layerElements]
  have hroot : Polynomial.eval (z i) (s.map (fun t => Polynomial.X - Polynomial.C t)).prod = 0 := by
    rw [Polynomial.eval_multiset_prod, Multiset.map_map, Multiset.prod_eq_zero]
    exact Multiset.mem_map.mpr ⟨z i, Multiset.mem_map.mpr ⟨i, Finset.mem_univ i, rfl⟩, by simp⟩
  rw [Multiset.prod_X_sub_X_eq_sum_esymm, Polynomial.eval_finset_sum, hcard] at hroot
  have hs (k : ℕ) : s.esymm k = elementary z k := rfl
  norm_num [Finset.sum_range_succ, hs, elementary_zero, hE1, hE2, hE4, hE5] at hroot
  linear_combination hroot

private def cubeFiber {n : ℕ} (z : Configuration n) (w : ℂ) : Finset (Fin n) :=
  open Classical in Finset.univ.filter (fun i => z i ^ 3 = w)

private theorem mem_cubeFiber {n : ℕ} {z : Configuration n} {i : Fin n} {w : ℂ} :
    i ∈ cubeFiber z w ↔ z i ^ 3 = w := by
  classical
  simp [cubeFiber]

private theorem cubeFiber_card_le_three {n : ℕ} {z : Configuration n}
    (hz : Function.Injective z) (w : ℂ) : (cubeFiber z w).card ≤ 3 := by
  classical
  by_cases hs : (cubeFiber z w).Nonempty
  · obtain ⟨i, hi⟩ := hs
    let g : Configuration 3 := ![z i, omega * z i, omega ^ 2 * z i]
    have hr (j : ↥(cubeFiber z w)) : ∃ k, z j.val = g k := by
      rcases cube_roots_same_orbit ((mem_cubeFiber.mp j.property).trans
        (mem_cubeFiber.mp hi).symm) with h | h | h
      · exact ⟨0, h⟩
      · exact ⟨1, h⟩
      · exact ⟨2, h⟩
    choose f hf using hr
    have hf' : Function.Injective f := by
      intro a b he
      apply Subtype.ext
      apply hz
      rw [hf a, hf b, he]
    simpa using Fintype.card_le_of_injective f hf'
  · simp [Finset.not_nonempty_iff_eq_empty.mp hs]

private theorem quadratic_root_cases {A B x y u : ℂ} (hxy : x ≠ y)
    (hx : x ^ 2 - A * x + B = 0) (hy : y ^ 2 - A * y + B = 0)
    (hu : u ^ 2 - A * u + B = 0) :
    u = x ∨ u = y := by
  have hf : (x - y) * (x + y - A) = 0 := by linear_combination hx - hy
  have hs := sub_eq_zero.mp ((mul_eq_zero.mp hf).resolve_left (sub_ne_zero.mpr hxy))
  have hp : B = x * y := by linear_combination hx - x * hs
  have he : (u - x) * (u - y) = 0 := by linear_combination hu - u * hs - hp
  exact (mul_eq_zero.mp he).imp sub_eq_zero.mp sub_eq_zero.mp

/-- The vanishing moments force exactly two cubic images, each attained
    three times. The cubic moment is three times their sum. -/
theorem six_layer_cubic_images {z : Configuration 6} (hz : Function.Injective z)
    (h1 : ∑ i, z i ^ 1 = 0) (h2 : ∑ i, z i ^ 2 = 0)
    (h4 : ∑ i, z i ^ 4 = 0) (h5 : ∑ i, z i ^ 5 = 0) :
    ∃ j : Fin 6, z 0 ^ 3 ≠ z j ^ 3 ∧
      (∀ i, z i ^ 3 = z 0 ^ 3 ∨ z i ^ 3 = z j ^ 3) ∧
      (∑ i, z i ^ 3) = 3 * (z 0 ^ 3 + z j ^ 3) := by
  classical
  have he : ∃ j, z 0 ^ 3 ≠ z j ^ 3 := by
    by_contra! h
    have he : cubeFiber z (z 0 ^ 3) = Finset.univ := by
      ext j
      simp [mem_cubeFiber, (h j).symm]
    have hk := cubeFiber_card_le_three hz (z 0 ^ 3)
    rw [he] at hk
    norm_num at hk
  obtain ⟨j, hj⟩ := he
  have ho (i : Fin 6) : z i ^ 3 = z 0 ^ 3 ∨ z i ^ 3 = z j ^ 3 :=
    quadratic_root_cases hj (six_roots_cube_quadratic z h1 h2 h4 h5 0)
      (six_roots_cube_quadratic z h1 h2 h4 h5 j) (six_roots_cube_quadratic z h1 h2 h4 h5 i)
  refine ⟨j, hj, ho, ?_⟩
  let a := cubeFiber z (z 0 ^ 3)
  let b := cubeFiber z (z j ^ 3)
  have hdis : Disjoint a b := by
    rw [Finset.disjoint_left]
    intro i ha hb
    exact hj ((mem_cubeFiber.mp ha).symm.trans (mem_cubeFiber.mp hb))
  have hunion : a ∪ b = Finset.univ := by
    ext i
    simpa only [a, b, Finset.mem_union, mem_cubeFiber, Finset.mem_univ, iff_true] using ho i
  have hc := Finset.card_union_of_disjoint hdis
  rw [hunion, Finset.card_univ, Fintype.card_fin] at hc
  have ha := cubeFiber_card_le_three hz (z 0 ^ 3)
  have hb := cubeFiber_card_le_three hz (z j ^ 3)
  change a.card ≤ 3 at ha
  change b.card ≤ 3 at hb
  have ha' : a.card = 3 := by omega
  have hb' : b.card = 3 := by omega
  have hsa : (∑ i ∈ a, z i ^ 3) = 3 * z 0 ^ 3 := by
    calc
      _ = ∑ _i ∈ a, z 0 ^ 3 := Finset.sum_congr rfl (fun _ hi => mem_cubeFiber.mp hi)
      _ = _ := by simp [ha', nsmul_eq_mul]
  have hsb : (∑ i ∈ b, z i ^ 3) = 3 * z j ^ 3 := by
    calc
      _ = ∑ _i ∈ b, z j ^ 3 := Finset.sum_congr rfl (fun _ hi => mem_cubeFiber.mp hi)
      _ = _ := by simp [hb', nsmul_eq_mul]
  rw [← hunion, Finset.sum_union hdis, hsa, hsb]
  ring

end Mordell
