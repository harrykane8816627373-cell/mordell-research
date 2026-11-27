import Mordell.Normalization
import Mordell.Rigidity
import Mathlib.Data.Fintype.EquivFin

/-! Extraction of the radius layers of a finite configuration. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def radiusLayer {n : ℕ} (z : Configuration n) (r : ℝ) : Finset (Fin n) :=
  open Classical in Finset.univ.filter (fun i => Complex.normSq (z i) = r)

theorem mem_radiusLayer {n : ℕ} {z : Configuration n} {r : ℝ} {i : Fin n} :
    i ∈ radiusLayer z r ↔ Complex.normSq (z i) = r := by
  classical
  simp [radiusLayer]

def layerMoment {n : ℕ} (z : Configuration n) (r : ℝ) (k : ℕ) : ℂ :=
  ∑ i ∈ radiusLayer z r, z i ^ k

theorem sum_three_layers {n : ℕ} {z : Configuration n} {a b c : ℝ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (hs : ∀ i, Complex.normSq (z i) = a ∨ Complex.normSq (z i) = b ∨ Complex.normSq (z i) = c)
    {M : Type*} [AddCommMonoid M] (f : Fin n → M) :
    (∑ i ∈ radiusLayer z a, f i) + (∑ i ∈ radiusLayer z b, f i) +
      (∑ i ∈ radiusLayer z c, f i) = ∑ i, f i := by
  classical
  simp only [radiusLayer, Finset.sum_filter, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rcases hs i with h | h | h
  · simp [h, hab, hac]
  · simp [h, hab.symm, hbc]
  · simp [h, hac.symm, hbc.symm]

theorem mass_of_three_layers {n : ℕ} {z : Configuration n} {a b c : ℝ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (hs : ∀ i, Complex.normSq (z i) = a ∨ Complex.normSq (z i) = b ∨ Complex.normSq (z i) = c) :
    ((radiusLayer z a).card : ℝ) * a + (radiusLayer z b).card * b + (radiusLayer z c).card * c =
      mass z := by
  have h (r : ℝ) : (∑ i ∈ radiusLayer z r, Complex.normSq (z i)) =
      ((radiusLayer z r).card : ℝ) * r := by
    calc
      _ = ∑ _i ∈ radiusLayer z r, r :=
        Finset.sum_congr rfl fun _ hi => mem_radiusLayer.mp hi
      _ = _ := by simp
  simpa only [h, mass] using sum_three_layers hab hac hbc hs (fun i => Complex.normSq (z i))

theorem count_three_layers {n : ℕ} {z : Configuration n} {a b c : ℝ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (hs : ∀ i, Complex.normSq (z i) = a ∨ Complex.normSq (z i) = b ∨ Complex.normSq (z i) = c) :
    (radiusLayer z a).card + (radiusLayer z b).card + (radiusLayer z c).card = n := by
  simpa using sum_three_layers hab hac hbc hs (fun _ => (1 : ℕ))

def enumerateLayer {n k : ℕ} {z : Configuration n} {r : ℝ}
    (h : (radiusLayer z r).card = k) : Fin k ≃ ↥(radiusLayer z r) :=
  (Fintype.equivFinOfCardEq (by simpa using h)).symm

def listedLayer {n k : ℕ} {z : Configuration n} {r : ℝ}
    (h : (radiusLayer z r).card = k) : Configuration k :=
  fun i => z (enumerateLayer h i).val

theorem normSq_listedLayer {n k : ℕ} {z : Configuration n} {r : ℝ}
    (h : (radiusLayer z r).card = k) (i : Fin k) :
    Complex.normSq (listedLayer h i) = r := by
  change Complex.normSq (z (enumerateLayer h i).val) = r
  exact mem_radiusLayer.mp (enumerateLayer h i).property

theorem moment_listedLayer {n k : ℕ} {z : Configuration n} {r : ℝ}
    (h : (radiusLayer z r).card = k) (d : ℕ) :
    (∑ i, listedLayer h i ^ d) = layerMoment z r d := by
  calc
    _ = ∑ j : ↥(radiusLayer z r), z j.val ^ d :=
      Equiv.sum_comp (enumerateLayer h) (fun j => z j.val ^ d)
    _ = _ := Finset.sum_coe_sort (radiusLayer z r) (fun j => z j ^ d)

theorem injective_listedLayer {n k : ℕ} {z : Configuration n} {r : ℝ}
    (h : (radiusLayer z r).card = k) (hz : Function.Injective z) :
    Function.Injective (listedLayer h) :=
  hz.comp (Subtype.val_injective.comp (enumerateLayer h).injective)

theorem mem_range_listedLayer {n k : ℕ} {z : Configuration n} {r : ℝ}
    (h : (radiusLayer z r).card = k) {i : Fin n} (hi : Complex.normSq (z i) = r) :
    z i ∈ Set.range (listedLayer h) := by
  let j : ↥(radiusLayer z r) := ⟨i, mem_radiusLayer.mpr hi⟩
  refine ⟨(enumerateLayer h).symm j, ?_⟩
  simp [listedLayer, j]

theorem sameOrbit_of_range_subset {n : ℕ} {z w : Configuration n} {s : ℂ}
    (hz : Function.Injective z) (hs : Complex.normSq s = 1)
    (h : Set.range z ⊆ Set.range (scale s w)) : SameOrbit z w := by
  classical
  have hchoose (i : Fin n) : ∃ j, z i = s * w j := by
    obtain ⟨j, hj⟩ := h (Set.mem_range_self i)
    exact ⟨j, hj.symm⟩
  choose f hf using hchoose
  have hi : Function.Injective f := by
    intro i j hij
    apply hz
    rw [hf i, hf j, hij]
  exact ⟨s, Equiv.ofBijective f ⟨hi, Finite.surjective_of_injective hi⟩, hs, hf⟩

/-- A projected feature vector with one removed direction determines a
    common complex multiple of the reference moment vector. -/
theorem proportional_of_projected_zero {ι : Type*} (pivot : ι) (m ref : ι → ℂ)
    (href : ref pivot ≠ 0)
    (hm : ∀ j, j ≠ pivot → m j - (ref j / ref pivot) * m pivot = 0) :
    ∃ η : ℂ, ∀ j, m j = η * ref j := by
  refine ⟨m pivot / ref pivot, fun j => ?_⟩
  by_cases hj : j = pivot
  · subst j
    rw [div_mul_cancel₀ _ href]
  · rw [sub_eq_zero.mp (hm j hj)]
    ring

def weightedMoment {n : ℕ} (z : Configuration n) (shift : ℝ) (d k : ℕ) : ℂ :=
  ∑ i, ((Complex.normSq (z i) - shift : ℝ) : ℂ) ^ d * z i ^ k

theorem weightedMoment_eq_of_three_layers {n : ℕ} {z : Configuration n} {a b c : ℝ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (hs : ∀ i, Complex.normSq (z i) = a ∨ Complex.normSq (z i) = b ∨ Complex.normSq (z i) = c)
    (shift : ℝ) (d k : ℕ) :
    weightedMoment z shift d k = ((a - shift : ℝ) : ℂ) ^ d * layerMoment z a k +
      ((b - shift : ℝ) : ℂ) ^ d * layerMoment z b k +
      ((c - shift : ℝ) : ℂ) ^ d * layerMoment z c k := by
  have hr (r : ℝ) :
      (∑ i ∈ radiusLayer z r, ((Complex.normSq (z i) - shift : ℝ) : ℂ) ^ d * z i ^ k) =
        ((r - shift : ℝ) : ℂ) ^ d * layerMoment z r k := by
    rw [layerMoment, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i hi => by rw [mem_radiusLayer.mp hi]
  simpa only [weightedMoment, hr] using
    (sum_three_layers hab hac hbc hs (fun i =>
      ((Complex.normSq (z i) - shift : ℝ) : ℂ) ^ d * z i ^ k)).symm

/-- The first three radial coefficients of a harmonic determine its layer moments. -/
theorem layerMoments_eq_of_weightedMoments {n : ℕ} {z : Configuration n} {a b c : ℝ}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (hs : ∀ i, Complex.normSq (z i) = a ∨ Complex.normSq (z i) = b ∨ Complex.normSq (z i) = c)
    (shift : ℝ) (k : ℕ) (A B C : ℂ)
    (h : ∀ d : Fin 3, weightedMoment z shift d.val k =
      ((a - shift : ℝ) : ℂ) ^ d.val * A + ((b - shift : ℝ) : ℂ) ^ d.val * B +
      ((c - shift : ℝ) : ℂ) ^ d.val * C) :
    layerMoment z a k = A ∧ layerMoment z b k = B ∧ layerMoment z c k = C := by
  have hd {r s : ℝ} (hne : r ≠ s) : ((r - shift : ℝ) : ℂ) ≠ ((s - shift : ℝ) : ℂ) := by
    intro he
    have hr : r - shift = s - shift := by exact_mod_cast he
    exact hne (sub_left_inj.mp hr)
  have he (d : Fin 3) := (weightedMoment_eq_of_three_layers hab hac hbc hs shift d.val k).symm.trans (h d)
  have h0 := he 0
  have h1 := he 1
  have h2 := he 2
  obtain ⟨ha, hb, hc⟩ := three_node_vandermonde (hd hab) (hd hac) (hd hbc)
    (A := layerMoment z a k - A) (B := layerMoment z b k - B) (C := layerMoment z c k - C)
    (by norm_num at h0; linear_combination h0)
    (by simp only [show (1 : Fin 3).val = 1 by decide, pow_one] at h1; linear_combination h1)
    (by simp only [show (2 : Fin 3).val = 2 by decide] at h2; linear_combination h2)
  exact ⟨sub_eq_zero.mp ha, sub_eq_zero.mp hb, sub_eq_zero.mp hc⟩

end Mordell
