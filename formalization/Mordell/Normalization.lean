import Mordell.Basic

/-! Compactness, scaling, and transfer of sharp variational bounds. -/

open scoped BigOperators

noncomputable section

namespace Mordell

def scale {n : ℕ} (s : ℂ) (z : Configuration n) : Configuration n :=
  fun i => s * z i

theorem mass_scale {n : ℕ} (s : ℂ) (z : Configuration n) :
    mass (scale s z) = Complex.normSq s * mass z := by
  simp only [mass, scale, Complex.normSq_mul, Finset.mul_sum]

theorem distanceSq_scale {n : ℕ} (s : ℂ) (z : Configuration n) (i j : Fin n) :
    distanceSq (scale s z) i j = Complex.normSq s * distanceSq z i j := by
  simp only [distanceSq, scale, ← mul_sub, Complex.normSq_mul]

def pairCount (n : ℕ) : ℕ :=
  ∑ i : Fin n, ∑ j : Fin n, if i < j then 1 else 0

theorem distanceProduct_scale {n : ℕ} (s : ℂ) (z : Configuration n) :
    distanceProduct (scale s z) =
      (Complex.normSq s) ^ pairCount n * distanceProduct z := by
  classical
  have hp (i j : Fin n) :
      (if i < j then distanceSq (scale s z) i j else 1) =
        (Complex.normSq s) ^ (if i < j then (1 : ℕ) else 0) *
          (if i < j then distanceSq z i j else 1) := by
    split_ifs <;> simp [distanceSq_scale]
  simp only [distanceProduct, hp, Finset.prod_mul_distrib, Finset.prod_pow_eq_pow_sum,
    pairCount]

theorem continuous_mass {n : ℕ} : Continuous (mass (n := n)) := by
  change Continuous (fun z : Configuration n => ∑ i, Complex.normSq (z i))
  simp only [Complex.normSq_eq_norm_sq]
  fun_prop

theorem continuous_distanceProduct {n : ℕ} : Continuous (distanceProduct (n := n)) := by
  classical
  unfold distanceProduct
  refine continuous_finset_prod _ fun i _ => continuous_finset_prod _ fun j _ => ?_
  by_cases h : i < j
  · simp only [h, if_true, distanceSq, Complex.normSq_eq_norm_sq]
    fun_prop
  · simp only [h, if_false]
    exact continuous_const

theorem compact_massLevel {n : ℕ} {q : ℝ} (hq : 0 ≤ q) :
    IsCompact {z : Configuration n | mass z = q} := by
  let ball : Fin n → Set ℂ := fun _ => Metric.closedBall 0 (Real.sqrt q)
  have hc : IsClosed {z : Configuration n | mass z = q} :=
    isClosed_eq continuous_mass continuous_const
  apply (isCompact_univ_pi (s := ball) fun _ => isCompact_closedBall _ _).of_isClosed_subset hc
  intro z hz
  rw [Set.mem_univ_pi]
  intro i
  change dist (z i) 0 ≤ Real.sqrt q
  rw [dist_zero_right]
  have hsum : Complex.normSq (z i) ≤ mass z :=
    Finset.single_le_sum (fun j _ => Complex.normSq_nonneg (z j)) (Finset.mem_univ i)
  have hz' : mass z = q := hz
  rw [hz', Complex.normSq_eq_norm_sq] at hsum
  have hs := Real.sq_sqrt hq
  nlinarith [norm_nonneg (z i), Real.sqrt_nonneg q]

theorem exists_maximizer {n : ℕ} {q : ℝ} (hq : 0 ≤ q)
    (w : Configuration n) (hw : mass w = q) :
    ∃ z : Configuration n, IsMaximizer q z := by
  obtain ⟨z, hz, hmax⟩ := (compact_massLevel (n := n) hq).exists_isMaxOn
    ⟨w, hw⟩ (continuous_distanceProduct (n := n)).continuousOn
  exact ⟨z, hz, fun u hu => hmax (show u ∈ {z : Configuration n | mass z = q} from hu)⟩

def SameOrbit {n : ℕ} (z w : Configuration n) : Prop :=
  ∃ (s : ℂ) (e : Fin n ≃ Fin n), Complex.normSq s = 1 ∧
    ∀ i, z i = s * w (e i)

theorem mass_relabel {n : ℕ} (z : Configuration n) (e : Fin n ≃ Fin n) :
    mass (z ∘ e) = mass z :=
  Equiv.sum_comp e (fun i => Complex.normSq (z i))

theorem distanceSq_symm {n : ℕ} (z : Configuration n) (i j : Fin n) :
    distanceSq z i j = distanceSq z j i := by
  unfold distanceSq
  rw [show z j - z i = -(z i - z j) by ring, Complex.normSq_neg]

theorem distanceProduct_relabel {n : ℕ} (z : Configuration n) (e : Fin n ≃ Fin n) :
    distanceProduct (z ∘ e) = distanceProduct z := by
  by_cases hi : Function.Injective z
  · have hi' := hi.comp e.injective
    apply Real.log_injOn_pos (distanceProduct_pos hi') (distanceProduct_pos hi)
    rw [log_distanceProduct hi', log_distanceProduct hi]
    let f (i j : Fin n) := Real.log (distanceSq z i j)
    have hf (i j : Fin n) : f i j = f j i := congrArg Real.log (distanceSq_symm z i j)
    have hsum := twice_pairSum (fun i j => f (e i) (e j)) (fun i j => hf (e i) (e j))
    have hordered : (∑ i, ∑ j, f (e i) (e j)) = ∑ i, ∑ j, f i j := by
      calc
        _ = ∑ i, ∑ j, f (e i) j :=
          Finset.sum_congr rfl (fun i _ => Equiv.sum_comp e (f (e i)))
        _ = _ := Equiv.sum_comp e (fun i => ∑ j, f i j)
    have hdiag : (∑ i, f (e i) (e i)) = ∑ i, f i i := Equiv.sum_comp e (fun i => f i i)
    rw [hordered, hdiag, ← twice_pairSum f hf] at hsum
    change pairSum (fun i j => f (e i) (e j)) = pairSum f
    linarith
  · have hj : ¬ Function.Injective (z ∘ e) := by
      intro h
      apply hi
      intro i j hij
      simpa using congrArg e (h (a₁ := e.symm i) (a₂ := e.symm j) (by simpa using hij))
    have hzero (w : Configuration n) (hw : ¬ Function.Injective w) :
        distanceProduct w = 0 :=
      le_antisymm (le_of_not_gt (fun h => hw (injective_of_distanceProduct_pos h)))
        (distanceProduct_nonneg w)
    rw [hzero _ hi, hzero _ hj]

theorem SameOrbit.mass_eq {n : ℕ} {z w : Configuration n} (h : SameOrbit z w) :
    mass z = mass w := by
  obtain ⟨s, e, hs, he⟩ := h
  have hrep : z = scale s (w ∘ e) := funext he
  rw [hrep, mass_scale, hs, one_mul, mass_relabel]

theorem SameOrbit.distanceProduct_eq {n : ℕ} {z w : Configuration n} (h : SameOrbit z w) :
    distanceProduct z = distanceProduct w := by
  obtain ⟨s, e, hs, he⟩ := h
  have hrep : z = scale s (w ∘ e) := funext he
  rw [hrep, distanceProduct_scale, hs, one_pow, one_mul, distanceProduct_relabel]

theorem SameOrbit.scale {n : ℕ} {z w : Configuration n} (h : SameOrbit z w) (s : ℂ) :
    SameOrbit (scale s z) (scale s w) := by
  obtain ⟨r, e, hr, he⟩ := h
  refine ⟨r, e, hr, fun i => ?_⟩
  change s * z i = r * (s * w (e i))
  rw [he]
  ring

def LogCertificate {n : ℕ} (q : ℝ) (w : Configuration n) : Prop :=
  ∀ z : Configuration n, IsMaximizer q z → Function.Injective z →
    Real.log (distanceProduct z) ≤ Real.log (distanceProduct w) ∧
      (Real.log (distanceProduct z) = Real.log (distanceProduct w) → SameOrbit z w)

/-- A logarithmic certificate on constrained maximizers implies the full
    sharp bound and its equality classification after rescaling. -/
theorem sharp_of_log_certificate {n : ℕ} {q : ℝ} {w : Configuration n}
    (hw : mass w = q) (hwi : Function.Injective w) {s : ℂ} (hs : s ≠ 0)
    (hcert : LogCertificate q w) :
    IsMaximizer (Complex.normSq s * q) (scale s w) ∧
      (∀ z : Configuration n, IsMaximizer (Complex.normSq s * q) z →
        SameOrbit z (scale s w)) := by
  have hq : 0 ≤ q := hw ▸ Finset.sum_nonneg (fun i _ => Complex.normSq_nonneg (w i))
  have hwpos := distanceProduct_pos hwi
  obtain ⟨m, hm⟩ := exists_maximizer hq w hw
  have hmpos := hwpos.trans_le (hm.2 w hw)
  obtain ⟨hlog, _⟩ := hcert m hm (injective_of_distanceProduct_pos hmpos)
  have hmle := (Real.log_le_log_iff hmpos hwpos).1 hlog
  have hwmax : IsMaximizer q w := ⟨hw, fun z hz => (hm.2 z hz).trans hmle⟩
  have hspos := Complex.normSq_pos.mpr hs
  have hrestore (z : Configuration n) : scale s (scale s⁻¹ z) = z := by
    funext i
    simp [scale, hs]
  have hundo (z : Configuration n) (hz : mass z = Complex.normSq s * q) :
      mass (scale s⁻¹ z) = q := by
    rw [mass_scale, Complex.normSq_inv, hz, ← mul_assoc, inv_mul_cancel₀ hspos.ne', one_mul]
  refine ⟨⟨by rw [mass_scale, hw], fun z hz => ?_⟩, fun z hz => ?_⟩
  · have he := mul_le_mul_of_nonneg_left (hwmax.2 _ (hundo z hz))
      (pow_nonneg hspos.le (pairCount n))
    rw [← distanceProduct_scale s, ← distanceProduct_scale s, hrestore] at he
    exact he
  · have hzq : IsMaximizer q (scale s⁻¹ z) := by
      refine ⟨hundo z hz.1, fun u hu => ?_⟩
      have hu' : mass (scale s u) = Complex.normSq s * q := by rw [mass_scale, hu]
      have he : distanceProduct (scale s u) ≤ distanceProduct (scale s (scale s⁻¹ z)) := by
        rw [hrestore]
        exact hz.2 _ hu'
      rw [distanceProduct_scale, distanceProduct_scale] at he
      exact le_of_mul_le_mul_left he (pow_pos hspos _)
    have hge := hzq.2 _ hw
    have hle := hwmax.2 _ hzq.1
    have hi := injective_of_distanceProduct_pos (hwpos.trans_le hge)
    obtain ⟨_, hr⟩ := hcert _ hzq hi
    have h := hr (congrArg Real.log (le_antisymm hle hge))
    simpa only [hrestore] using h.scale s

end Mordell
