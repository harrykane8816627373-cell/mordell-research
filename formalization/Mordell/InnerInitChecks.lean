import Mordell.InnerInitData

/-! Exact arithmetic for the source tail, root power coefficients and three transformations. -/

namespace Mordell.InnerCoverData

def sourceCap : Fin 3 → ℚ := ![6 / 25, 6 / 25, 3 / 40]

def zeroIndex : GridIndex degree := fun _ => 0

set_option maxRecDepth 25000
set_option maxHeartbeats 0

theorem tail_checked : (tailNumerator : ℚ) / (stageDen0 : ℚ) = discardedLower sourceCap sourceOutside := by
  decide +kernel

theorem input_checked : ∀ i : GridIndex degree, stage0 i = inputMultiplier * sourceInside i +
    (if i 0 = 0 ∧ i 1 = 0 ∧ i 2 = 0 then tailNumerator else 0) := by
  decide +kernel

theorem matrix0_checked : ∀ i j : Fin (degree 0 + 1), (axisMatrix0 i j : ℚ) / (axisDen0 : ℚ) =
    axisPowerControl ℚ (degree 0) (root.lower 0) (root.upper 0) i.val j.val := by decide +kernel

theorem matrix1_checked : ∀ i j : Fin (degree 1 + 1), (axisMatrix1 i j : ℚ) / (axisDen1 : ℚ) =
    axisPowerControl ℚ (degree 1) (root.lower 1) (root.upper 1) i.val j.val := by decide +kernel

theorem matrix2_checked : ∀ i j : Fin (degree 2 + 1), (axisMatrix2 i j : ℚ) / (axisDen2 : ℚ) =
    axisPowerControl ℚ (degree 2) (root.lower 2) (root.upper 2) i.val j.val := by decide +kernel

theorem stage1_checked : ∀ i : GridIndex degree, stage1 i = tensorAxisMap 0 axisMatrix0 stage0 i := by
  decide +kernel

theorem stage2_checked : ∀ i : GridIndex degree, stage2 i = tensorAxisMap 1 axisMatrix1 stage1 i := by
  decide +kernel

theorem stage3_checked : ∀ i : GridIndex degree, stage3 i = tensorAxisMap 2 axisMatrix2 stage2 i := by
  decide +kernel

theorem rounding_checked : ∀ i : GridIndex degree, gridLane degree grid0 i * (stageDen3 : ℤ) ≤
    stage3 i * (2 : ℤ) ^ (32 : ℕ) := by decide +kernel

end Mordell.InnerCoverData
