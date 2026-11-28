import Mordell.GramData

/-! Kernel replay of the corrected Gram boxes in both sharp proofs. -/

namespace Mordell

set_option maxRecDepth 20000 in
set_option maxHeartbeats 0 in
private theorem tenGramData_checked (i : Fin 14) : (tenGramData i).2.widthsOK ∧
    ∀ j, (tenGramData i).2.rowOK j := by
  fin_cases i <;> decide +kernel

set_option maxRecDepth 20000 in
set_option maxHeartbeats 0 in
private theorem twelveGramData_checked (i : Fin 12) : (twelveGramData i).2.widthsOK ∧
    ∀ j, (twelveGramData i).2.rowOK j := by
  fin_cases i <;> decide +kernel

/-- Twelve harmonic blocks and the two ten-point second-variation blocks. -/
theorem ten_gram_box_posDef (i : Fin 14) {A : Matrix (Fin (tenGramData i).1) (Fin (tenGramData i).1) ℝ}
    (hs : A.IsHermitian) (h : ∀ j k, ((tenGramData i).2.lower j k : ℝ) ≤ A j k ∧
      A j k ≤ ((tenGramData i).2.upper j k : ℝ)) : A.PosDef :=
  (tenGramData i).2.sound (tenGramData_checked i).1 (tenGramData_checked i).2 hs h

/-- The twelve corrected harmonic blocks of the twelve-point certificate. -/
theorem twelve_gram_box_posDef (i : Fin 12)
    {A : Matrix (Fin (twelveGramData i).1) (Fin (twelveGramData i).1) ℝ}
    (hs : A.IsHermitian) (h : ∀ j k, ((twelveGramData i).2.lower j k : ℝ) ≤ A j k ∧
      A j k ≤ ((twelveGramData i).2.upper j k : ℝ)) : A.PosDef :=
  (twelveGramData i).2.sound (twelveGramData_checked i).1 (twelveGramData_checked i).2 hs h

end Mordell
