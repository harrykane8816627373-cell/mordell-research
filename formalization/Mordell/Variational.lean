import Mordell.PairGeometry
import Mathlib.Analysis.Calculus.LocalExtr.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Deriv

/-! First-order calculus for the finite distance product. -/

open scoped BigOperators ComplexConjugate Topology

noncomputable section

namespace Mordell

def affinePath {n : ℕ} (z v : Configuration n) (t : ℝ) : Configuration n :=
  fun i => z i + (t : ℂ) * v i

@[simp] theorem affinePath_zero {n : ℕ} (z v : Configuration n) : affinePath z v 0 = z := by
  ext i
  simp [affinePath]

theorem continuous_affinePath {n : ℕ} (z v : Configuration n) : Continuous (affinePath z v) := by
  apply continuous_pi
  intro i
  change Continuous (fun t : ℝ => z i + (t : ℂ) * v i)
  fun_prop

private theorem affinePath_valid_near_zero {n : ℕ} {z : Configuration n}
    (hz : Function.Injective z) (hq : 0 < mass z) (v : Configuration n) :
    ∀ᶠ t in 𝓝 (0 : ℝ), Function.Injective (affinePath z v t) ∧ 0 < mass (affinePath z v t) := by
  have hp : IsOpen {t : ℝ | 0 < distanceProduct (affinePath z v t)} :=
    isOpen_lt continuous_const (continuous_distanceProduct.comp (continuous_affinePath z v))
  have hm : IsOpen {t : ℝ | 0 < mass (affinePath z v t)} :=
    isOpen_lt continuous_const (continuous_mass.comp (continuous_affinePath z v))
  have hp0 : (0 : ℝ) ∈ {t : ℝ | 0 < distanceProduct (affinePath z v t)} := by
    simpa using distanceProduct_pos hz
  have hm0 : (0 : ℝ) ∈ {t : ℝ | 0 < mass (affinePath z v t)} := by simpa using hq
  filter_upwards [hp.mem_nhds hp0, hm.mem_nhds hm0] with t ht ht'
  exact ⟨injective_of_distanceProduct_pos ht, ht'⟩

private theorem normSq_affine_hasDerivAt (a v : ℂ) :
    HasDerivAt (fun t : ℝ => Complex.normSq (a + (t : ℂ) * v))
      (2 * realDot v a) 0 := by
  have hr := (((hasDerivAt_id (0 : ℝ)).mul_const v.re).const_add a.re).pow 2
  have hi := (((hasDerivAt_id (0 : ℝ)).mul_const v.im).const_add a.im).pow 2
  convert hr.add hi using 1
  · funext t
    simp only [Complex.normSq_apply, Complex.add_re, Complex.add_im, Complex.mul_re,
      Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero, sub_zero, id_eq,
      Pi.add_apply, Pi.pow_apply]
    ring
  · simp only [id_eq, zero_mul, add_zero, realDot,
      Complex.mul_re, Complex.conj_re, Complex.conj_im]
    norm_num
    ring

private theorem mass_affine_hasDerivAt {n : ℕ} (z v : Configuration n) :
    HasDerivAt (fun t => mass (affinePath z v t))
      (2 * ∑ i, realDot (v i) (z i)) 0 := by
  rw [Finset.mul_sum]
  exact HasDerivAt.fun_sum fun i _ => normSq_affine_hasDerivAt (z i) (v i)

private theorem log_normSq_affine_hasDerivAt {a : ℂ} (ha : a ≠ 0) (v : ℂ) :
    HasDerivAt (fun t : ℝ => Real.log (Complex.normSq (a + (t : ℂ) * v)))
      (2 * (v / a).re) 0 := by
  have h := normSq_affine_hasDerivAt a v
  have hn : Complex.normSq (a + ((0 : ℝ) : ℂ) * v) ≠ 0 := by
    simpa using (Complex.normSq_pos.mpr ha).ne'
  convert h.log hn using 1
  simp only [Complex.ofReal_zero, zero_mul, add_zero, Complex.div_re, realDot, Complex.mul_re,
    Complex.conj_im, Complex.conj_re]
  ring

private def logDirection {n : ℕ} (z v : Configuration n) (i j : Fin n) : ℝ :=
  ((v i - v j) / (z i - z j)).re

private theorem logDistance_affine_hasDerivAt {n : ℕ} {z v : Configuration n} {i j : Fin n}
    (h : z i ≠ z j) :
    HasDerivAt (fun t => Real.log (distanceSq (affinePath z v t) i j))
      (2 * logDirection z v i j) 0 := by
  convert log_normSq_affine_hasDerivAt (sub_ne_zero.mpr h) (v i - v j) using 1
  funext t
  simp only [affinePath, distanceSq]
  congr 2
  ring

private theorem pairSum_logDirection_eq {n : ℕ} (z v : Configuration n) :
    pairSum (logDirection z v) = (∑ i, v i * force z i).re := by
  have hsym (i j : Fin n) : logDirection z v i j = logDirection z v j i := by
    unfold logDirection
    rw [show v j - v i = -(v i - v j) by ring, show z j - z i = -(z i - z j) by ring,
      neg_div_neg_eq]
  have hdiag (i : Fin n) : logDirection z v i i = 0 := by simp [logDirection]
  have h := twice_pairSum (logDirection z v) hsym
  rw [show (∑ i, logDirection z v i i) = 0 by simp [hdiag], sub_zero] at h
  have hr := congrArg Complex.re (force_pairing z v)
  have hreal (f : Fin n → ℂ) : (∑ i, f i).re = ∑ i, (f i).re :=
    map_sum Complex.reAddGroupHom f Finset.univ
  have hentry (i j : Fin n) :
      (if i = j then (0 : ℂ) else (v i - v j) * (z i - z j)⁻¹).re = logDirection z v i j := by
    by_cases hi : i = j
    · subst j
      simp [hdiag]
    · simp only [if_neg hi, logDirection, div_eq_mul_inv]
  change ((2 : ℂ) * (∑ i, v i * force z i)).re =
      (∑ i : Fin n, ∑ j : Fin n, if i = j then (0 : ℂ) else (v i - v j) * (z i - z j)⁻¹).re at hr
  rw [hreal] at hr
  simp only [hreal, hentry] at hr
  rw [Complex.mul_re] at hr
  norm_num at hr
  simp only [hreal, Complex.mul_re, Finset.sum_sub_distrib]
  linarith

private theorem log_product_affine_hasDerivAt {n : ℕ} {z : Configuration n}
    (hz : Function.Injective z) (hq : 0 < mass z) (v : Configuration n) :
    HasDerivAt (fun t => Real.log (distanceProduct (affinePath z v t)))
      (2 * (∑ i, v i * force z i).re) 0 := by
  have hpair : HasDerivAt (fun t => pairSum (fun i j => Real.log (distanceSq (affinePath z v t) i j)))
      (2 * pairSum (logDirection z v)) 0 := by
    have hscale (i j : Fin n) :
        (if i < j then 2 * logDirection z v i j else 0) =
          2 * (if i < j then logDirection z v i j else 0) := by split <;> simp
    have hs : HasDerivAt (fun t => pairSum (fun i j => Real.log (distanceSq (affinePath z v t) i j)))
        (∑ i : Fin n, ∑ j : Fin n, if i < j then 2 * logDirection z v i j else 0) 0 := by
      apply HasDerivAt.fun_sum
      intro i _
      apply HasDerivAt.fun_sum
      intro j _
      by_cases hij : i < j
      · simp only [hij, if_true]
        exact logDistance_affine_hasDerivAt (hz.ne hij.ne)
      · simp only [hij, if_false]
        exact hasDerivAt_const _ _
    simpa only [hscale, ← Finset.mul_sum, pairSum] using hs
  rw [pairSum_logDirection_eq] at hpair
  apply hpair.congr_of_eventuallyEq
  filter_upwards [affinePath_valid_near_zero hz hq v] with t ht
  exact log_distanceProduct ht.1

private theorem penalized_affine_hasDerivAt {n : ℕ} {q : ℝ} {z : Configuration n}
    (hz : Function.Injective z) (hq : 0 < mass z) (v : Configuration n) :
    HasDerivAt (fun t => penalizedLogProduct q (affinePath z v t))
      (2 * (∑ i, v i * forceResidual ((pairCount n : ℝ) / q) z i).re) 0 := by
  let m := (pairCount n : ℝ) / q
  have h := (log_product_affine_hasDerivAt hz hq v).sub
    ((mass_affine_hasDerivAt z v).const_mul m)
  have him (f : Fin n → ℂ) : (∑ i, f i).re = ∑ i, (f i).re :=
    map_sum Complex.reAddGroupHom f Finset.univ
  convert h using 1
  · simp only [forceResidual, mul_sub, Finset.sum_sub_distrib, Complex.sub_re, him]
    have hterm (i : Fin n) : (v i * ((m : ℂ) * conj (z i))).re = m * realDot (v i) (z i) := by
      rw [show v i * ((m : ℂ) * conj (z i)) = (m : ℂ) * (v i * conj (z i)) by ring]
      simp [realDot]
    rw [Finset.sum_congr rfl (fun i _ => hterm i), ← Finset.mul_sum]
    ring

/-- Every positive-product maximizer satisfies the full complex force
    equations. The test direction is the conjugate force residual. -/
private theorem stationary_of_injective_maximizer {n : ℕ} {q : ℝ} {z : Configuration n}
    (hq : 0 < q) (hz : IsMaximizer q z) (hi : Function.Injective z) :
    Stationary (((pairCount n : ℝ) / q : ℝ) : ℂ) z := by
  let m := (pairCount n : ℝ) / q
  let v (i : Fin n) := conj (forceResidual m z i)
  have hmass : 0 < mass z := hz.1.symm ▸ hq
  have hl : IsLocalMax (fun t => penalizedLogProduct q (affinePath z v t)) 0 := by
    filter_upwards [affinePath_valid_near_zero hi hmass v] with t ht
    simpa only [affinePath_zero] using penalizedLogProduct_le_of_maximizer hq hz ht.1 ht.2
  have he := hl.hasDerivAt_eq_zero (penalized_affine_hasDerivAt (q := q) hi hmass v)
  have hreal (f : Fin n → ℂ) : (∑ i, f i).re = ∑ i, (f i).re :=
    map_sum Complex.reAddGroupHom f Finset.univ
  have hp (i : Fin n) :
      (v i * forceResidual m z i).re = Complex.normSq (forceResidual m z i) := by
    change (conj (forceResidual m z i) * forceResidual m z i).re = _
    rw [← Complex.normSq_eq_conj_mul_self, Complex.ofReal_re]
  change 2 * (∑ i, v i * forceResidual m z i).re = 0 at he
  rw [hreal, Finset.sum_congr rfl (fun i _ => hp i)] at he
  have ht : ∑ i, Complex.normSq (forceResidual m z i) = 0 := by linarith
  intro i
  have hn : Complex.normSq (forceResidual m z i) ≤ ∑ j, Complex.normSq (forceResidual m z j) :=
    Finset.single_le_sum (fun _ _ => Complex.normSq_nonneg _) (Finset.mem_univ i)
  rw [ht] at hn
  exact sub_eq_zero.mp (Complex.normSq_eq_zero.mp
    (le_antisymm hn (Complex.normSq_nonneg _)))

theorem injective_of_maximizer {n : ℕ} {q : ℝ} {z : Configuration n}
    (hn : 1 < n) (hq : 0 < q) (hz : IsMaximizer q z) : Function.Injective z := by
  let u : Configuration n := fun i => (i.val : ℂ)
  have hu : Function.Injective u := by
    intro i j hij
    apply Fin.ext
    change (i.val : ℂ) = (j.val : ℂ) at hij
    exact_mod_cast hij
  let i : Fin n := ⟨1, hn⟩
  have hu1 : u i = 1 := by simp [u, i]
  have hmu : 0 < mass u := by
    have h := Finset.single_le_sum (s := Finset.univ)
      (f := fun j => Complex.normSq (u j)) (fun j _ => Complex.normSq_nonneg (u j))
      (Finset.mem_univ i)
    change Complex.normSq (u i) ≤ mass u at h
    rw [hu1, map_one] at h
    exact zero_lt_one.trans_le h
  let s : ℂ := (Real.sqrt (q / mass u) : ℝ)
  have hs : Complex.normSq s = q / mass u := by
    change Complex.normSq ((Real.sqrt (q / mass u) : ℝ) : ℂ) = _
    rw [Complex.normSq_ofReal, ← sq, Real.sq_sqrt (div_pos hq hmu).le]
  have hsn : s ≠ 0 := Complex.normSq_pos.mp (hs.symm ▸ div_pos hq hmu)
  have hi : Function.Injective (scale s u) := fun _ _ h => hu (mul_left_cancel₀ hsn h)
  have hm : mass (scale s u) = q := by
    rw [mass_scale, hs, div_mul_cancel₀ _ hmu.ne']
  exact injective_of_distanceProduct_pos ((distanceProduct_pos hi).trans_le (hz.2 _ hm))

/-- Every maximizer at positive mass satisfies the complex first-order equations. -/
theorem stationary_of_maximizer {n : ℕ} {q : ℝ} {z : Configuration n}
    (hn : 1 < n) (hq : 0 < q) (hz : IsMaximizer q z) :
    Stationary (((pairCount n : ℝ) / q : ℝ) : ℂ) z :=
  stationary_of_injective_maximizer hq hz (injective_of_maximizer hn hq hz)

theorem stationary_at_mass_n_sub_one {n : ℕ} {z : Configuration n} (hn : 1 < n)
    (hz : IsMaximizer ((n : ℝ) - 1) z) : Stationary ((n : ℂ) / 2) z := by
  have hn' : 0 < (n : ℝ) - 1 := sub_pos.mpr (by exact_mod_cast hn)
  have hcount : ((pairCount n : ℝ) / ((n : ℝ) - 1)) = (n : ℝ) / 2 := by
    have h := twice_pairSum (fun _ _ : Fin n => (1 : ℝ)) (fun _ _ => rfl)
    have hc : (pairCount n : ℝ) = pairSum (fun _ _ : Fin n => (1 : ℝ)) := by
      simp only [pairCount, pairSum, Nat.cast_sum]
      refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
      split_ifs <;> norm_num
    rw [← hc] at h
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one] at h
    apply (div_eq_iff hn'.ne').mpr
    linear_combination h / 2
  have h := stationary_of_maximizer hn hn' hz
  rw [hcount, Complex.ofReal_div, Complex.ofReal_natCast] at h
  norm_num at h ⊢
  exact h

end Mordell
