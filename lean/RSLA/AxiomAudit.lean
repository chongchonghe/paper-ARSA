import RSLA.CFL
import RSLA.TwoStage
import RSLA.Strict
import RSLA.WeightedNorm
import RSLA.M1Full
import RSLA.M1UnequalFull
import Lean.Util.CollectAxioms

/-! Fail compilation if any RSLA theorem depends on an axiom other than Lean's
standard classical foundations. This checks dependencies transitively, including
imported theorems used in a proof. -/
open Lean in
run_cmd do
  let allowed : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  let mut count := 0
  for (name, info) in (← getEnv).constants.toList do
    if (`RSLA).isPrefixOf name && info.isTheorem then
      count := count + 1
      for ax in ← collectAxioms name do
        unless allowed.contains ax do
          throwError "{name} depends on disallowed axiom {ax}"
  logInfo m!"Axiom audit passed for {count} RSLA theorems: only propext, Classical.choice, and Quot.sound are permitted."
