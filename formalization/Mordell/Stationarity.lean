import Mordell.Normalization
import Mathlib.Algebra.Ring.GeomSum

/-! The penalized objective and algebraic force moments shared by both proofs. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

def penalizedLogProduct {n : ℕ} (q : ℝ) (z : Configuration n) : ℝ :=
  Real.log (distanceProduct z) - (pairCount n : ℝ) / q * mass z

/-- A constrained maximizer already dominates on every positive mass after
    subtracting the quadratic penalty. -/
theorem penalizedLogProduct_le_of_maximizer {n : ℕ} {q : ℝ} {z w : Configuration n}
    (hq : 0 < q) (hz : IsMaximizer q z) (hw : Function.Injective w) (hmw : 0 < mass w) :
    penalizedLogProduct q w ≤ penalizedLogProduct q z := by
  let s : ℂ := (Real.sqrt (q / mass w) : ℝ)
  have hs : Complex.normSq s = q / mass w := by
    dsimp [s]
    rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt (div_pos hq hmw).le]
  have hscaled : mass (scale s w) = q := by
    rw [mass_scale, hs, div_mul_cancel₀ _ hmw.ne']
  have hN : 0 ≤ (pairCount n : ℝ) := Nat.cast_nonneg _
  have hp := distanceProduct_pos hw
  have hsp : 0 < (q / mass w) ^ pairCount n := pow_pos (div_pos hq hmw) _
  have he := hz.2 _ hscaled
  rw [distanceProduct_scale, hs] at he
  have h := (Real.log_le_log_iff (mul_pos hsp hp) ((mul_pos hsp hp).trans_le he)).2 he
  rw [Real.log_mul hsp.ne' hp.ne', Real.log_pow, Real.log_div hq.ne' hmw.ne'] at h
  have hl := Real.log_le_sub_one_of_pos (div_pos hmw hq)
  rw [Real.log_div hmw.ne' hq.ne'] at hl
  have hmajor := mul_le_mul_of_nonneg_left hl hN
  unfold penalizedLogProduct
  rw [hz.1]
  have hcancel : (pairCount n : ℝ) / q * q = pairCount n := div_mul_cancel₀ _ hq.ne'
  rw [hcancel]
  have hnormalize : (pairCount n : ℝ) / q * mass w = (pairCount n : ℝ) * (mass w / q) := by ring
  rw [hnormalize]
  linarith only [h, hmajor]

def force {n : ℕ} (z : Configuration n) (i : Fin n) : ℂ :=
  ∑ j, if i = j then 0 else (z i - z j)⁻¹

def Stationary {n : ℕ} (multiplier : ℂ) (z : Configuration n) : Prop :=
  ∀ i, force z i = multiplier * conj (z i)

theorem Stationary.scale {n : ℕ} {multiplier : ℂ} {z : Configuration n}
    (h : Stationary multiplier z) (s : ℂ) :
    Stationary (multiplier / (Complex.normSq s : ℂ)) (Mordell.scale s z) := by
  have hf (i : Fin n) : force (Mordell.scale s z) i = s⁻¹ * force z i := by
    simp only [force, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    by_cases hij : i = j
    · simp [hij]
    · simp only [if_neg hij, Mordell.scale, ← mul_sub, mul_inv_rev]
      ring
  intro i
  rw [hf, h i, Complex.inv_def]
  simp only [Mordell.scale, map_mul, ← Complex.ofReal_inv, div_eq_mul_inv]
  ring

def powerMoment {n : ℕ} (z : Configuration n) (k : ℕ) : ℂ := ∑ i, z i ^ k

theorem force_pairing {n : ℕ} (z : Configuration n) (v : Fin n → ℂ) :
    2 * (∑ i, v i * force z i) =
      ∑ i, ∑ j, if i = j then 0 else (v i - v j) * (z i - z j)⁻¹ := by
  classical
  let f := fun i j => if i = j then (0 : ℂ) else (z i - z j)⁻¹
  have ha (i j : Fin n) : f j i = -f i j := by
    by_cases h : i = j
    · subst j
      simp [f]
    · simp only [f, if_neg h, if_neg (Ne.symm h)]
      rw [show z j - z i = -(z i - z j) by ring, inv_neg]
  have hswap : (∑ i, ∑ j, v j * f i j) = -(∑ i, v i * force z i) := by
    rw [Finset.sum_comm, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    change (∑ j, v i * f j i) = -(v i * ∑ j, f i j)
    rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun j _ => by rw [ha i j, mul_neg]
  have hpoint (i j : Fin n) :
      (if i = j then (0 : ℂ) else (v i - v j) * (z i - z j)⁻¹) =
        v i * f i j - v j * f i j := by
    split_ifs with h <;> simp [f, h, sub_mul]
  simp only [hpoint, Finset.sum_sub_distrib, ← Finset.mul_sum, hswap]
  change 2 * (∑ i, v i * force z i) = (∑ i, v i * force z i) - -(∑ i, v i * force z i)
  ring

theorem centered_of_stationary {n : ℕ} {z : Configuration n} {multiplier : ℂ}
    (hmultiplier : multiplier ≠ 0) (hz : Stationary multiplier z) :
    ∑ i, z i = 0 := by
  have h := force_pairing z (fun _ => 1)
  unfold Stationary at hz
  simp only [one_mul, sub_self, zero_mul, ite_self, Finset.sum_const_zero, hz, ← Finset.mul_sum,
    ← map_sum] at h
  rcases mul_eq_zero.mp h with h2 | hh
  · norm_num at h2
  · have hh' := (mul_eq_zero.mp hh).resolve_left hmultiplier
    have hsum := congrArg (starRingEnd ℂ) hh'
    simpa using hsum

/-- The divided-difference Ward identity with the diagonal subtracted
    from the ordered sum. -/
theorem ward_identity {n : ℕ} {z : Configuration n} {multiplier : ℂ}
    (hz : Function.Injective z) (hf : Stationary multiplier z) (k : ℕ) :
    2 * multiplier * (∑ i, (Complex.normSq (z i) : ℂ) * z i ^ k) =
      ∑ r ∈ Finset.range (k + 1), (powerMoment z r * powerMoment z (k - r) - powerMoment z k) := by
  classical
  let summand := fun i j r => z i ^ r * z j ^ (k - r)
  have hdiv (i j : Fin n) (hij : i ≠ j) :
      (z i ^ (k + 1) - z j ^ (k + 1)) * (z i - z j)⁻¹ =
        ∑ r ∈ Finset.range (k + 1), summand i j r := by
    rw [← geom_sum₂_mul (z i) (z j) (k + 1), mul_assoc,
      mul_inv_cancel₀ (sub_ne_zero.mpr (hz.ne hij)), mul_one]
    rfl
  have hdiag (i : Fin n) (r : ℕ) (hr : r ∈ Finset.range (k + 1)) :
      summand i i r = z i ^ k := by
    have hr' : r ≤ k := by have := Finset.mem_range.mp hr; omega
    simp only [summand, ← pow_add, Nat.add_sub_of_le hr']
  have hsplit (i j : Fin n) :
      (if i = j then (0 : ℂ) else (z i ^ (k + 1) - z j ^ (k + 1)) * (z i - z j)⁻¹) =
        ∑ r ∈ Finset.range (k + 1), (summand i j r - if j = i then z i ^ k else 0) := by
    by_cases h : i = j
    · subst j
      simp only [if_true]
      rw [Finset.sum_eq_zero fun r hr => sub_eq_zero.mpr (hdiag i r hr)]
    · simp [h, (Ne.symm h), hdiv i j h]
  have hleft :
      (∑ i, z i ^ (k + 1) * force z i) =
        multiplier * ∑ i, (Complex.normSq (z i) : ℂ) * z i ^ k := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [hf i, Complex.normSq_eq_conj_mul_self, pow_succ]
    ring
  calc
    _ = 2 * (∑ i, z i ^ (k + 1) * force z i) := by rw [hleft, mul_assoc]
    _ = ∑ i, ∑ j, ∑ r ∈ Finset.range (k + 1),
        (summand i j r - if j = i then z i ^ k else 0) := by
      simp only [force_pairing, hsplit]
    _ = ∑ r ∈ Finset.range (k + 1), ∑ i, ∑ j,
        (summand i j r - if j = i then z i ^ k else 0) := by
      calc
        _ = ∑ i, ∑ r ∈ Finset.range (k + 1), ∑ j,
            (summand i j r - if j = i then z i ^ k else 0) :=
          Finset.sum_congr rfl fun _ _ => Finset.sum_comm
        _ = _ := Finset.sum_comm
    _ = _ := by
      refine Finset.sum_congr rfl fun r _ => ?_
      simp only [summand, Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.sum_mul,
        Finset.sum_ite_eq', Finset.mem_univ, if_true, powerMoment]

def forceResidual {n : ℕ} (multiplier : ℝ) (z : Configuration n) (i : Fin n) : ℂ :=
  force z i - (multiplier : ℂ) * conj (z i)

theorem weighted_residual_im_zero {n : ℕ} {z : Configuration n}
    (hz : Function.Injective z) (multiplier : ℝ) :
    (∑ i, z i * forceResidual multiplier z i).im = 0 := by
  have hpoint (i j : Fin n) :
      (if i = j then (0 : ℂ) else (z i - z j) * (z i - z j)⁻¹) =
        if i = j then 0 else 1 := by
    by_cases h : i = j
    · simp [h]
    · simp [h, sub_ne_zero.mpr (hz.ne h)]
  have him (f : Fin n → ℂ) : (∑ i, f i).im = ∑ i, (f i).im :=
    map_sum Complex.imAddGroupHom f Finset.univ
  have hsum := force_pairing z z
  have hsum' : 2 * (∑ i, z i * force z i) =
      ∑ i : Fin n, ∑ j : Fin n, if i = j then (0 : ℂ) else 1 := by
    simpa only [hpoint] using hsum
  have h01 (i j : Fin n) : (if i = j then (0 : ℂ) else 1).im = 0 := by split <;> simp
  have hforce : ((2 : ℂ) * (∑ i, z i * force z i)).im =
      (∑ i : Fin n, ∑ j : Fin n, if i = j then (0 : ℂ) else 1).im :=
    congrArg Complex.im hsum'
  have hzero :
      (∑ i : Fin n, ∑ j : Fin n, if i = j then (0 : ℂ) else 1).im = 0 := by
    simp only [him, h01, Finset.sum_const_zero]
  rw [hzero, Complex.mul_im] at hforce
  norm_num at hforce
  have hr (i : Fin n) : z i * ((multiplier : ℂ) * conj (z i)) =
      (multiplier : ℂ) * Complex.normSq (z i) := by rw [← Complex.mul_conj]; ring
  simp only [forceResidual, mul_sub, hr, Finset.sum_sub_distrib, Complex.sub_im, him,
    ← Complex.ofReal_mul, Complex.ofReal_im, Finset.sum_const_zero, sub_zero]
  simpa only [Complex.mul_im] using hforce

/-- In the real-axis gauge, the weighted residual identity recovers the
    omitted imaginary force equation. -/
theorem stationary_of_gauge_residual {n : ℕ} {z : Configuration n} {multiplier : ℝ}
    (hz : Function.Injective z) (pivot : Fin n) (hreal : (z pivot).im = 0)
    (hnonzero : (z pivot).re ≠ 0) (hpr : (forceResidual multiplier z pivot).re = 0)
    (hother : ∀ i, i ≠ pivot → forceResidual multiplier z i = 0) :
    Stationary (multiplier : ℂ) z := by
  have h := weighted_residual_im_zero hz multiplier
  have hsingle : (∑ i, z i * forceResidual multiplier z i) =
      z pivot * forceResidual multiplier z pivot := by
    apply Finset.sum_eq_single pivot
    · intro i _ hi
      rw [hother i hi, mul_zero]
    · simp
  rw [hsingle, Complex.mul_im, hreal, zero_mul, add_zero] at h
  have hpi := (mul_eq_zero.mp h).resolve_left hnonzero
  have hp : forceResidual multiplier z pivot = 0 :=
    Complex.ext (by simpa using hpr) (by simpa using hpi)
  intro i
  have hres : forceResidual multiplier z i = 0 := by
    by_cases hi : i = pivot
    · simpa only [hi] using hp
    · exact hother i hi
  exact sub_eq_zero.mp hres

end Mordell
