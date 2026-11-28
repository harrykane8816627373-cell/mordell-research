import Mordell.Stationarity
import Mathlib.Logic.Equiv.Fin.Basic

/-! The real coordinate order used by the archived stationary-root certificates. -/

open scoped BigOperators ComplexConjugate

noncomputable section

namespace Mordell

abbrev GaugeCoordinates (n : ℕ) := Fin (n * 2 + 1) → ℝ

def gaugeRead {S : Type*} [Zero S] {n : ℕ} (x : Fin (n * 2 + 1) → S)
    (i : Fin (n + 1)) (a : Fin 2) : S :=
  Fin.cases (if a = 0 then x 0 else 0) (fun j => x (finProdFinEquiv (j, a)).succ) i

/-- `x₀`, then consecutive real/imaginary pairs, with `y₀ = 0`. -/
def gaugeConfiguration {n : ℕ} (x : GaugeCoordinates n) : Configuration (n + 1) :=
  fun i => ⟨gaugeRead x i 0, gaugeRead x i 1⟩

def complexComponent (z : ℂ) (a : Fin 2) : ℝ := if a = 0 then z.re else z.im

theorem component_gaugeConfiguration {n : ℕ} (x : GaugeCoordinates n) (i : Fin (n + 1)) (a : Fin 2) :
    complexComponent (gaugeConfiguration x i) a = gaugeRead x i a := by
  fin_cases a <;> simp [complexComponent, gaugeConfiguration]

theorem gaugeRead_ratCast {n : ℕ} (x : Fin (n * 2 + 1) → ℚ) (i : Fin (n + 1)) (a : Fin 2) :
    gaugeRead (fun j => (x j : ℝ)) i a = ((gaugeRead (S := ℚ) x i a : ℚ) : ℝ) := by
  refine Fin.cases ?_ (fun j => ?_) i
  · by_cases h : a = 0 <;> simp [gaugeRead, h]
  · rfl

def gaugeProjection {n : ℕ} (z : Configuration (n + 1)) : GaugeCoordinates n :=
  Fin.cases (z 0).re (fun k =>
    let index : Fin n × Fin 2 := finProdFinEquiv.symm k
    complexComponent (z index.1.succ) index.2)

theorem gaugeProjection_apply_pair {n : ℕ} (z : Configuration (n + 1)) (j : Fin n) (a : Fin 2) :
    gaugeProjection z (finProdFinEquiv (j, a)).succ = complexComponent (z j.succ) a := by
  simp only [gaugeProjection, Fin.cases_succ, Equiv.symm_apply_apply]

theorem gaugeProjection_gaugeConfiguration {n : ℕ} (x : GaugeCoordinates n) :
    gaugeProjection (gaugeConfiguration x) = x := by
  funext i
  refine Fin.cases ?_ (fun j => ?_) i
  · rfl
  · obtain ⟨⟨k, a⟩, rfl⟩ := finProdFinEquiv.surjective j
    rw [gaugeProjection_apply_pair]
    fin_cases a <;> simp [component_gaugeConfiguration, gaugeRead]

/-- The Cartesian force is the conjugate of the complex force residual. -/
def cartesianResidual {n : ℕ} (s : ℝ) (z : Configuration n) (i : Fin n) : ℂ :=
  (∑ j, if i = j then 0 else (z i - z j) / (distanceSq z i j : ℂ)) - (s : ℂ) * z i

theorem cartesianResidual_eq {n : ℕ} (s : ℝ) (z : Configuration n) (i : Fin n) :
    cartesianResidual s z i = conj (forceResidual s z i) := by
  have hpair (j : Fin n) : ((z i - z j) / (distanceSq z i j : ℂ)) = conj ((z i - z j)⁻¹) := by
    conv_rhs => rw [Complex.inv_def]
    rw [map_mul, Complex.conj_conj, Complex.conj_ofReal, Complex.ofReal_inv, div_eq_mul_inv, distanceSq]
  simp only [cartesianResidual, forceResidual, force, map_sub, map_sum, map_mul,
    Complex.conj_ofReal, Complex.conj_conj, hpair]
  congr 1
  exact Finset.sum_congr rfl (fun j _ => by split <;> simp_all)

theorem cartesianResidual_component {n : ℕ} (s : ℝ) (z : Configuration n) (i : Fin n) (a : Fin 2) :
    complexComponent (cartesianResidual s z i) a =
      (∑ j, if i = j then 0 else (complexComponent (z i) a - complexComponent (z j) a) /
        distanceSq z i j) - s * complexComponent (z i) a := by
  have hre (f : Fin n → ℂ) : (∑ j, f j).re = ∑ j, (f j).re :=
    map_sum Complex.reAddGroupHom f Finset.univ
  have him (f : Fin n → ℂ) : (∑ j, f j).im = ∑ j, (f j).im :=
    map_sum Complex.imAddGroupHom f Finset.univ
  fin_cases a <;>
    simp only [cartesianResidual, complexComponent, Fin.isValue, Fin.zero_eta, Fin.reduceFinMk,
      Fin.one_eq_zero_iff, Nat.reduceEqDiff, if_true, if_false, Complex.sub_re, Complex.sub_im,
      Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul,
      sub_zero, add_zero, hre, him] <;>
    congr 1 <;> exact Finset.sum_congr rfl (fun j _ => by split <;> simp)

def gaugeForce {n : ℕ} (s : ℝ) (x : GaugeCoordinates n) : GaugeCoordinates n :=
  gaugeProjection (cartesianResidual s (gaugeConfiguration x))

/-- The gauge projection of the rational real force has precisely the full
    stationary configurations as zeros when the first coordinate is nonzero. -/
theorem gaugeForce_eq_zero_iff {n : ℕ} (s : ℝ) (x : GaugeCoordinates n)
    (hx : x 0 ≠ 0) (hno : Function.Injective (gaugeConfiguration x)) :
    gaugeForce s x = 0 ↔ Stationary (s : ℂ) (gaugeConfiguration x) := by
  constructor
  · intro hf
    have hzero (i : Fin (n * 2 + 1)) : gaugeForce s x i = 0 := congrFun hf i
    apply stationary_of_gauge_residual hno 0 (by simp [gaugeConfiguration, gaugeRead]) hx
    · simpa [gaugeForce, gaugeProjection, cartesianResidual_eq] using hzero 0
    · intro i hi
      obtain ⟨j, rfl⟩ := Fin.exists_succ_eq.mpr hi
      apply Complex.ext
      · simpa [gaugeForce, gaugeProjection_apply_pair, complexComponent, cartesianResidual_eq] using
          hzero (finProdFinEquiv (j, (0 : Fin 2))).succ
      · simpa [gaugeForce, gaugeProjection_apply_pair, complexComponent, cartesianResidual_eq] using
          hzero (finProdFinEquiv (j, (1 : Fin 2))).succ
  · intro h
    have hres : cartesianResidual s (gaugeConfiguration x) = 0 := by
      funext i
      simp [cartesianResidual_eq, forceResidual, h i]
    funext i
    simp only [gaugeForce, hres]
    refine Fin.cases ?_ (fun k => ?_) i
    · simp [gaugeProjection]
    · simp [gaugeProjection, complexComponent]

theorem mass_twelve_of_stationary {z : Configuration 12} (hz : Function.Injective z)
    (hs : Stationary (11 / 2 : ℂ) z) : mass z = 12 := by
  have h := ward_identity hz hs 0
  have hm : (∑ i, (Complex.normSq (z i) : ℂ)) = (mass z : ℂ) := by simp [mass]
  norm_num [powerMoment, Finset.sum_range_succ, hm] at h
  have he := congrArg Complex.re h
  simp [Complex.mul_re] at he
  linarith

end Mordell
