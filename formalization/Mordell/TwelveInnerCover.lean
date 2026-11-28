import Mordell.InnerCoverChecks

/-! Kernel replay of the 39-node positive-angle inner-inner subtree. -/

namespace Mordell.InnerCoverData

noncomputable section

theorem trace_matches : String.ofList tree.traceCode = traceCode := by decide +kernel

/-- The signed inner-inner transcript proves the log comparison given a lower
    Bernstein tensor for the pair kernel at the root. -/
theorem signed_inner_log_bound {H : (Fin 3 → ℝ) → ℝ} (hH : GridEnvelope degree root grid0 H)
    {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    Real.log (gapDistance x) ≤ H x := by
  have hd' : gramPairDomain ((1 / 3 : ℚ) : ℝ) true x := by simpa using hd
  exact tree.sound (by intro j; fin_cases j <;> decide) (by norm_num) archived_tree_accepted hH x hx hd'

def rationalLowerKernel (x : Fin 3 → ℝ) : ℝ :=
  tensorBernstein degree (fun i => (gridControls degree grid0 i : ℝ))
    (fun j => (x j - (root.lower j : ℝ)) / ((root.upper j : ℝ) - root.lower j))

private theorem rationalLowerKernel_envelope : GridEnvelope degree root grid0 rationalLowerKernel := by
  intro u _
  have hw (j : Fin 3) : (root.upper j : ℝ) - root.lower j ≠ 0 := by
    fin_cases j <;> norm_num [root]
  have he : (fun j => (((root.lower j : ℝ) + ((root.upper j : ℝ) - root.lower j) * u j) - root.lower j) /
      ((root.upper j : ℝ) - root.lower j)) = u := by
    funext j
    field_simp [hw j]
    ring
  unfold rationalLowerKernel
  rw [he]

/-- The concrete rational Bernstein polynomial bounds the logarithm on this signed box. -/
theorem signed_inner_rational_bound {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    Real.log (gapDistance x) ≤ rationalLowerKernel x :=
  signed_inner_log_bound rationalLowerKernel_envelope hx hd

theorem signed_inner_kernel_bound {H : (Fin 3 → ℝ) → ℝ}
    (hH : ∀ x, root.Contains x → gramPairDomain (1 / 3) true x → rationalLowerKernel x ≤ H x)
    {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    Real.log (gapDistance x) ≤ H x :=
  (signed_inner_rational_bound hx hd).trans (hH x hx hd)

/-- The physical distance floor bounds the angular coordinate more tightly than the signed root box. -/
theorem signed_inner_angular_cap {x : Fin 3 → ℝ} (hx : root.Contains x) (hd : gramPairDomain (1 / 3) true x) :
    0 ≤ x 2 ∧ x 2 ≤ (3 / 40 : ℝ) := by
  have hu : x 0 ≤ ((6 / 25 : ℚ) : ℝ) := (hx 0).2
  have hv : x 1 ≤ ((6 / 25 : ℚ) : ℝ) := (hx 1).2
  have ht := (hx 2).1
  change (root.lower 2 : ℝ) ≤ x 2 at ht
  rw [show root.lower 2 = 0 by decide +kernel, Rat.cast_zero] at ht
  have hdist := hd.2.2.2.1
  norm_num [gapDistance] at hu hv hdist
  exact ⟨ht, by linarith⟩

end
end Mordell.InnerCoverData
