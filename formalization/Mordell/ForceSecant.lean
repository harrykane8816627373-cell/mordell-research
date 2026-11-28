import Mordell.GaugeRoot
import Mathlib.Data.Matrix.Mul

/-! Exact matrix secants for the Cartesian force on the rotation slice. -/

open scoped BigOperators Matrix

noncomputable section

namespace Mordell

def gaugeAxis {n : ℕ} : Fin (n * 2 + 1) → Fin (n + 1) × Fin 2 :=
  Fin.cases (0, 0) (fun k =>
    let ij : Fin n × Fin 2 := finProdFinEquiv.symm k
    (ij.1.succ, ij.2))

theorem gaugeForce_apply {n : ℕ} (s : ℝ) (x : GaugeCoordinates n) (p : Fin (n * 2 + 1)) :
    gaugeForce s x p = complexComponent (cartesianResidual s (gaugeConfiguration x) (gaugeAxis p).1)
      (gaugeAxis p).2 := by
  refine Fin.cases ?_ (fun k => ?_) p
  · simp [gaugeForce, gaugeProjection, gaugeAxis, complexComponent]
  · rfl

def pairDifferences {n : ℕ} (x : GaugeCoordinates n) (i j : Fin (n + 1)) (a : Fin 2) : ℝ :=
  gaugeRead x i a - gaugeRead x j a

def pairSecant {n : ℕ} (x y : GaugeCoordinates n) (i j : Fin (n + 1)) (a b : Fin 2) : ℝ :=
  if i = j then 0 else
    let d := distanceSq (gaugeConfiguration x) i j
    let e := distanceSq (gaugeConfiguration y) i j
    (if a = b then d⁻¹ else 0) -
      pairDifferences y i j a * (pairDifferences x i j b + pairDifferences y i j b) / (d * e)

private theorem pairSecant_apply {n : ℕ} (x y : GaugeCoordinates n) (i j : Fin (n + 1)) (a : Fin 2)
    (hx : Function.Injective (gaugeConfiguration x)) (hy : Function.Injective (gaugeConfiguration y)) :
    (∑ b, pairSecant x y i j a b * (pairDifferences x i j b - pairDifferences y i j b)) =
      if i = j then 0 else pairDifferences x i j a / distanceSq (gaugeConfiguration x) i j -
        pairDifferences y i j a / distanceSq (gaugeConfiguration y) i j := by
  by_cases hij : i = j
  · simp [pairSecant, hij]
  let d := distanceSq (gaugeConfiguration x) i j
  let e := distanceSq (gaugeConfiguration y) i j
  have hd : d ≠ 0 := (distanceSq_pos hx hij).ne'
  have he : e ≠ 0 := (distanceSq_pos hy hij).ne'
  have hdiff : (∑ b, (pairDifferences x i j b + pairDifferences y i j b) *
      (pairDifferences x i j b - pairDifferences y i j b)) = d - e := by
    simp only [Fin.sum_univ_two, pairDifferences, d, e, distanceSq, Complex.normSq_apply,
      Complex.sub_re, Complex.sub_im, gaugeConfiguration]
    ring
  have hdiag : (∑ b : Fin 2, (if a = b then d⁻¹ else 0) *
      (pairDifferences x i j b - pairDifferences y i j b)) =
      (pairDifferences x i j a - pairDifferences y i j a) / d := by
    simp_rw [ite_mul, zero_mul]
    simp [div_eq_mul_inv, mul_comm]
  have hsecond : (∑ b, (pairDifferences y i j a * (pairDifferences x i j b + pairDifferences y i j b) /
      (d * e)) * (pairDifferences x i j b - pairDifferences y i j b)) =
      pairDifferences y i j a / (d * e) * (d - e) := by
    rw [← hdiff, Finset.mul_sum]
    exact Finset.sum_congr rfl fun _ _ => by ring
  simp only [pairSecant, if_neg hij, sub_mul, Finset.sum_sub_distrib]
  change (∑ b : Fin 2, (if a = b then d⁻¹ else 0) *
      (pairDifferences x i j b - pairDifferences y i j b)) -
    (∑ b : Fin 2, (pairDifferences y i j a * (pairDifferences x i j b + pairDifferences y i j b) /
      (d * e)) * (pairDifferences x i j b - pairDifferences y i j b)) =
    pairDifferences x i j a / d - pairDifferences y i j a / e
  rw [hdiag, hsecond]
  field_simp [hd, he]
  ring

def gaugeIncidence {S : Type*} [Zero S] [One S] {n : ℕ}
    (i : Fin (n + 1)) (a : Fin 2) (p : Fin (n * 2 + 1)) : S :=
  gaugeRead (fun q => if q = p then 1 else 0) i a

theorem gaugeIncidence_cast {n : ℕ} (i : Fin (n + 1)) (a : Fin 2) (p : Fin (n * 2 + 1)) :
    (↑(gaugeIncidence (S := ℤ) i a p) : ℝ) = gaugeIncidence i a p := by
  refine Fin.cases ?_ (fun j => ?_) i
  · by_cases ha : a = 0 <;> by_cases hp : (0 : Fin (n * 2 + 1)) = p <;>
      simp [gaugeIncidence, gaugeRead, ha, hp]
  · by_cases hp : (finProdFinEquiv (j, a)).succ = p <;> simp [gaugeIncidence, gaugeRead, hp]

theorem abs_gaugeIncidence_sub_le_one {n : ℕ} (i j : Fin (n + 1)) (a : Fin 2) (p : Fin (n * 2 + 1)) :
    |(gaugeIncidence i a p : ℝ) - gaugeIncidence j a p| ≤ 1 := by
  have h (k : Fin (n + 1)) : (gaugeIncidence k a p : ℝ) = 0 ∨ (gaugeIncidence k a p : ℝ) = 1 := by
    refine Fin.cases ?_ (fun l => ?_) k
    · by_cases ha : a = 0 <;> by_cases hp : (0 : Fin (n * 2 + 1)) = p <;>
        simp [gaugeIncidence, gaugeRead, ha, hp]
    · by_cases hp : (finProdFinEquiv (l, a)).succ = p <;> simp [gaugeIncidence, gaugeRead, hp]
  rcases h i with hi | hi <;> rcases h j with hj | hj <;> norm_num [hi, hj]

private theorem incidence_apply {n : ℕ} (x : GaugeCoordinates n) (i : Fin (n + 1)) (a : Fin 2) :
    ∑ p, gaugeIncidence i a p * x p = gaugeRead x i a := by
  refine Fin.cases ?_ (fun j => ?_) i
  · by_cases h : a = 0
    · simp [gaugeIncidence, gaugeRead, h, ite_mul]
    · simp [gaugeIncidence, gaugeRead, h]
  · simp [gaugeIncidence, gaugeRead, ite_mul]

private theorem gaugeRead_sub {n : ℕ} (x y : GaugeCoordinates n) (i : Fin (n + 1)) (a : Fin 2) :
    gaugeRead (x - y) i a = gaugeRead x i a - gaugeRead y i a := by
  refine Fin.cases ?_ (fun j => ?_) i
  · by_cases h : a = 0 <;> simp [gaugeRead, h]
  · rfl

def gaugeSecant {n : ℕ} (s : ℝ) (x y : GaugeCoordinates n) : Matrix (Fin (n * 2 + 1)) (Fin (n * 2 + 1)) ℝ :=
  fun p q =>
    let i := (gaugeAxis p).1
    let a := (gaugeAxis p).2
    (∑ j, ∑ b, pairSecant x y i j a b * (gaugeIncidence i b q - gaugeIncidence j b q)) -
      s * gaugeIncidence i a q

theorem gaugeForce_sub_eq_secant {n : ℕ} (s : ℝ) (x y : GaugeCoordinates n)
    (hx : Function.Injective (gaugeConfiguration x)) (hy : Function.Injective (gaugeConfiguration y)) :
    gaugeForce s x - gaugeForce s y = gaugeSecant s x y *ᵥ (x - y) := by
  funext p
  let i := (gaugeAxis p).1
  let a := (gaugeAxis p).2
  have hi : (gaugeAxis p).1 = i := rfl
  have ha : (gaugeAxis p).2 = a := rfl
  simp only [Pi.sub_apply, gaugeForce_apply, hi, ha, cartesianResidual_component,
    component_gaugeConfiguration]
  have hmatrix : (gaugeSecant s x y *ᵥ (x - y)) p =
      (∑ j, ∑ b, pairSecant x y i j a b * (gaugeRead (x - y) i b - gaugeRead (x - y) j b)) -
        s * gaugeRead (x - y) i a := by
    simp only [Matrix.mulVec, dotProduct, gaugeSecant, hi, ha, sub_mul, Finset.sum_sub_distrib, Finset.sum_mul]
    rw [show (∑ q, s * gaugeIncidence i a q * (x - y) q) = s * gaugeRead (x - y) i a by
      simp only [mul_assoc, ← Finset.mul_sum, incidence_apply]]
    congr 1
    calc
      _ = ∑ j, ∑ b, ∑ q, pairSecant x y i j a b * ((gaugeIncidence i b q - gaugeIncidence j b q) * (x - y) q) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun _ _ => by simp only [mul_assoc]; exact Finset.sum_comm
      _ = _ := by
        simp only [← Finset.mul_sum, sub_mul, Finset.sum_sub_distrib, incidence_apply]
  rw [hmatrix]
  have hq (j : Fin (n + 1)) (b : Fin 2) : gaugeRead (x - y) i b - gaugeRead (x - y) j b =
      pairDifferences x i j b - pairDifferences y i j b := by
    simp only [gaugeRead_sub, pairDifferences]
    ring
  simp only [hq]
  have hpair (j : Fin (n + 1)) := pairSecant_apply x y i j a hx hy
  simp only [hpair, gaugeRead_sub]
  have he (j : Fin (n + 1)) : (if i = j then 0 else pairDifferences x i j a /
      distanceSq (gaugeConfiguration x) i j - pairDifferences y i j a / distanceSq (gaugeConfiguration y) i j) =
      (if i = j then 0 else pairDifferences x i j a / distanceSq (gaugeConfiguration x) i j) -
      (if i = j then 0 else pairDifferences y i j a / distanceSq (gaugeConfiguration y) i j) := by
    split_ifs <;> simp
  simp only [he, Finset.sum_sub_distrib]
  simp only [pairDifferences]
  ring

end Mordell
