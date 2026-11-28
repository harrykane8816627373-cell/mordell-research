import Mordell.PairGeometry
import Mathlib.Analysis.Calculus.DerivativeTest
import Mathlib.Analysis.Calculus.LocalExtr.Basic
import Mathlib.Analysis.Complex.RealDeriv
import Mathlib.Analysis.SpecialFunctions.Log.Deriv

/-! Variational calculus for the finite distance product. -/

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

private theorem normSq_affine_hasDerivAt (a v : ℂ) (t : ℝ) :
    HasDerivAt (fun t : ℝ => Complex.normSq (a + (t : ℂ) * v))
      (2 * realDot v (a + (t : ℂ) * v)) t := by
  have hr := (((hasDerivAt_id t).mul_const v.re).const_add a.re).pow 2
  have hi := (((hasDerivAt_id t).mul_const v.im).const_add a.im).pow 2
  convert hr.add hi using 1
  · funext t
    simp only [Complex.normSq_apply, Complex.add_re, Complex.add_im, Complex.mul_re,
      Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero, sub_zero, id_eq,
      Pi.add_apply, Pi.pow_apply]
    ring
  · simp only [id_eq, realDot, Complex.mul_re, Complex.conj_re, Complex.conj_im, Complex.add_re,
      Complex.add_im, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero, sub_zero]
    ring

private theorem mass_affine_hasDerivAt {n : ℕ} (z v : Configuration n) (t : ℝ) :
    HasDerivAt (fun t => mass (affinePath z v t))
      (2 * ∑ i, realDot (v i) (affinePath z v t i)) t := by
  rw [Finset.mul_sum]
  exact HasDerivAt.fun_sum fun i _ => normSq_affine_hasDerivAt (z i) (v i) t

private theorem log_normSq_affine_hasDerivAt (a v : ℂ) {t : ℝ} (ha : a + (t : ℂ) * v ≠ 0) :
    HasDerivAt (fun t : ℝ => Real.log (Complex.normSq (a + (t : ℂ) * v)))
      (2 * (v / (a + (t : ℂ) * v)).re) t := by
  have h := normSq_affine_hasDerivAt a v t
  have hn := (Complex.normSq_pos.mpr ha).ne'
  convert h.log hn using 1
  generalize a + (t : ℂ) * v = w
  simp only [Complex.div_re, realDot, Complex.mul_re,
    Complex.conj_im, Complex.conj_re]
  ring

private def logDirection {n : ℕ} (z v : Configuration n) (i j : Fin n) : ℝ :=
  ((v i - v j) / (z i - z j)).re

private theorem affinePath_sub {n : ℕ} (z v : Configuration n) (t : ℝ) (i j : Fin n) :
    affinePath z v t i - affinePath z v t j = z i - z j + (t : ℂ) * (v i - v j) := by
  unfold affinePath
  ring

private theorem logDistance_affine_hasDerivAt {n : ℕ} {z v : Configuration n} {i j : Fin n} {t : ℝ}
    (h : affinePath z v t i ≠ affinePath z v t j) :
    HasDerivAt (fun s => Real.log (distanceSq (affinePath z v s) i j))
      (2 * logDirection (affinePath z v t) v i j) t := by
  have hn : z i - z j + (t : ℂ) * (v i - v j) ≠ 0 := by
    rw [← affinePath_sub]
    exact sub_ne_zero.mpr h
  simpa only [distanceSq, logDirection, affinePath_sub] using
    log_normSq_affine_hasDerivAt (z i - z j) (v i - v j) hn

private theorem pairSum_hasDerivAt {n : ℕ} {f : ℝ → Fin n → Fin n → ℝ} {d : Fin n → Fin n → ℝ} {t : ℝ}
    (h : ∀ i j, i < j → HasDerivAt (fun s => f s i j) (d i j) t) :
    HasDerivAt (fun s => pairSum (f s)) (pairSum d) t := by
  unfold pairSum
  apply HasDerivAt.fun_sum
  intro i _
  apply HasDerivAt.fun_sum
  intro j _
  by_cases hij : i < j
  · simpa only [hij, if_true] using h i j hij
  · simpa only [hij, if_false] using hasDerivAt_const t 0

private theorem pairSum_const_mul {n : ℕ} (c : ℝ) (f : Fin n → Fin n → ℝ) :
    pairSum (fun i j => c * f i j) = c * pairSum f := by
  unfold pairSum
  simp only [Finset.mul_sum, mul_ite, mul_zero]

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
    have hd (i j : Fin n) (hij : i < j) :=
      logDistance_affine_hasDerivAt (z := z) (v := v) (i := i) (j := j) (t := 0)
        (by simpa only [affinePath_zero] using hz.ne hij.ne)
    simpa only [pairSum_const_mul, affinePath_zero] using pairSum_hasDerivAt hd
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
    ((mass_affine_hasDerivAt z v 0).const_mul m)
  have him (f : Fin n → ℂ) : (∑ i, f i).re = ∑ i, (f i).re :=
    map_sum Complex.reAddGroupHom f Finset.univ
  convert h using 1
  · simp only [forceResidual, mul_sub, Finset.sum_sub_distrib, Complex.sub_re, him, affinePath_zero]
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

private theorem logDirection_affine_hasDerivAt {n : ℕ} {z v : Configuration n} {i j : Fin n}
    (h : z i ≠ z j) :
    HasDerivAt (fun t : ℝ => logDirection (affinePath z v t) v i j)
      (-(((v i - v j) / (z i - z j)) ^ 2).re) 0 := by
  let a := z i - z j
  let b := v i - v j
  have ha : a ≠ 0 := sub_ne_zero.mpr h
  have hd : HasDerivAt (fun w : ℂ => a + w * b) b ((0 : ℝ) : ℂ) := by
    simpa using ((hasDerivAt_id ((0 : ℝ) : ℂ)).mul_const b).const_add a
  have hq := ((hasDerivAt_const ((0 : ℝ) : ℂ) b).div hd (by simpa using ha)).real_of_complex
  simp only [logDirection, affinePath_sub, div_pow]
  simpa only [Pi.div_apply, Complex.ofReal_zero, zero_mul, add_zero, zero_sub, neg_div,
    Complex.neg_re, pow_two, a, b] using hq

private theorem realDot_affine (a v : ℂ) (t : ℝ) :
    realDot v (a + (t : ℂ) * v) = realDot v a + t * Complex.normSq v := by
  simp only [realDot, Complex.mul_re, Complex.conj_re, Complex.conj_im, Complex.add_re,
    Complex.add_im, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero,
    sub_zero, Complex.normSq_apply]
  ring

/-- At a positive-mass maximizer the real second variation gives a nonnegative
    quadratic form in every complex displacement. -/
theorem secondVariation_nonneg {n : ℕ} {q : ℝ} {z : Configuration n} (hn : 1 < n)
    (hq : 0 < q) (hz : IsMaximizer q z) (v : Configuration n) :
    0 ≤ (pairCount n : ℝ) / q * mass v +
      pairSum (fun i j => (((v i - v j) / (z i - z j)) ^ 2).re) := by
  let m := (pairCount n : ℝ) / q
  let f (t : ℝ) := pairSum (fun i j => Real.log (distanceSq (affinePath z v t) i j)) -
    m * mass (affinePath z v t)
  let g (t : ℝ) := 2 * pairSum (logDirection (affinePath z v t) v) -
    2 * m * ((∑ i, realDot (v i) (z i)) + t * mass v)
  let Q := pairSum (fun i j => (((v i - v j) / (z i - z j)) ^ 2).re)
  have hi := injective_of_maximizer hn hq hz
  have hm : 0 < mass z := hz.1.symm ▸ hq
  have hv := affinePath_valid_near_zero hi hm v
  have hl : IsLocalMax f 0 := by
    filter_upwards [hv] with t ht
    simpa only [f, ← log_distanceProduct ht.1, affinePath_zero, ← log_distanceProduct hi] using
      penalizedLogProduct_le_of_maximizer hq hz ht.1 ht.2
  have hf (t : ℝ) (ht : Function.Injective (affinePath z v t)) : HasDerivAt f (g t) t := by
    have hpair := pairSum_hasDerivAt (fun i j hij =>
      logDistance_affine_hasDerivAt (t := t) (ht.ne (ne_of_lt hij)))
    have hd := hpair.sub ((mass_affine_hasDerivAt z v t).const_mul m)
    convert hd using 1
    simp only [g, pairSum_const_mul, affinePath, realDot_affine,
      Finset.sum_add_distrib, ← Finset.mul_sum, mass]
    ring
  have he : deriv f =ᶠ[𝓝 0] g := hv.mono fun t ht => (hf t ht.1).deriv
  have hg : HasDerivAt g (-2 * (m * mass v + Q)) 0 := by
    have hpair := pairSum_hasDerivAt (fun i j hij =>
      logDirection_affine_hasDerivAt (v := v) (hi.ne (ne_of_lt hij)))
    have hsum := (((hasDerivAt_id (0 : ℝ)).mul_const (mass v)).const_add
      (∑ i, realDot (v i) (z i))).const_mul (2 * m)
    have hd := (hpair.const_mul 2).sub hsum
    have hQ : pairSum (fun i j => -(((v i - v j) / (z i - z j)) ^ 2).re) = -Q := by
      simpa only [neg_one_mul, Q] using pairSum_const_mul (-1) (fun i j =>
        (((v i - v j) / (z i - z j)) ^ 2).re)
    convert hd using 1
    simp only [hQ, one_mul]
    ring
  have hs : deriv (deriv f) 0 = -2 * (m * mass v + Q) :=
    he.deriv_eq.trans hg.deriv
  by_contra! h
  have hpos : 0 < deriv (deriv f) 0 := by
    rw [hs]
    change m * mass v + Q < 0 at h
    linarith
  have hmin := isLocalMin_of_deriv_deriv_pos hpos hl.deriv_eq_zero (hf 0 (by simpa using hi)).continuousAt
  have hc : f =ᶠ[𝓝 0] fun _ => f 0 := by
    filter_upwards [hl, hmin] with t ht ht'
    exact le_antisymm ht ht'
  have hzero := hc.deriv.deriv_eq
  simp only [deriv_const', deriv_const] at hzero
  exact hpos.ne' hzero

end Mordell
