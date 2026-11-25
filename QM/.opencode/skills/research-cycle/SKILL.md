---
name: research-cycle
description: Screen new quantum-research proposals, plan decisive tests, and resume or checkpoint substantial investigations without repeating rejected ideas.
---

# Evidence-first research cycle

Use only the relevant steps. This is not a mandatory report template or a reason
to reread unchanged notes. Simple questions and Git synchronization bypass it.
All workspace paths below are relative to the workspace root.

## 1. Recover the actual question

- Use `notes/INDEX.md` to find the existing topic note and its latest checkpoint.
  Read `notes/RESEARCH_CONTEXT.md` only for missing current scope.
- Preserve the user's requested object and acceptance criterion. Distinguish a
  theorem, a candidate mechanism, a diagnostic, a numerical test and a literature lead.
- Before proposing a new direction, consult the relevant section of
  `notes/matejs-notes-research-map.md` and the topic's decision ledger.
  Inspect original notebook equations only where the comparison requires them.
- Separate "already in the notebook", "new combination here" and "literature novelty
  not checked". A new name is not a new mechanism. Reopen a failed approach only
  with a changed assumption, new evidence or a genuinely different argument.
- Ask one focused question only if a missing choice changes the scientific problem.
  Otherwise proceed with a stated, reversible assumption; do not ask the user to
  repeat choices already recorded.

## 2. Screen before investing

For the candidate under consideration, identify briefly:
- evolving object and base dynamics; perturbation and its consistent derivative;
- metric/readout, what stays fixed, and what is optimized;
- exact versus effective dynamics, finite versus asymptotic time, and the proposed claim;
- nearest known obstruction and the cheapest test that could reject the claim.

For the current Lyapunov programme, use
`notes/open-problems/quantum-lyapunov-programme.md` for authoritative scope.
In particular, do not substitute conditioned or mean-field dynamics for exact
closed-system dynamics, or a finite-time indicator for the requested construction.
Readout nonclosure alone is not a rejection if exact rho is retained underneath.
Visibility gain, singular normalization and coordinate/metric changes require
interpretation; test a regular control before interpreting growth as chaos.
Probability-value correspondence does not establish tangent correspondence.

Prefer one discriminating next step over a long list of possibilities, unless
the user requests a survey. State what either outcome would and would not settle.

## 3. Verify at the appropriate cost

- Reuse existing proofs, source annotations and scripts. Check cited hypotheses;
  retrieve primary equations when a conclusion depends on them.
- For experiments, read `calculations/README.md`; use a toy case/pilot before a
  sweep. Fix controls and criteria before interpreting results. Record commands,
  parameters, seeds, tolerances, code version/hash and actual artifact paths.
- Separate feasibility, convergence and optimality. Label repairs and provisional
  results; never promote an earlier failed or unconverged run to a final result.
- Delegate only an independent check or substantial bounded task; include the
  question, scope, evidence paths, allowed edits and required validation.
- Stop on a decisive obstruction or missing essential input. Record what would
  justify reopening; do not retry the same idea with renamed variables.

## 4. Leave a restart point, not another transcript

Only quantum-research updates canonical notes. At a substantial stopping point,
update the existing topic note, retaining evidence and correcting superseded claims.
Reuse its decision/next-step section; add a short checkpoint only if needed:

> Objective and fixed scope:
> Current result and status:
> Decisive evidence / code and final artifact paths:
> Rejected route and reason, if new:
> Blocker or unresolved question:
> Next discriminating step and what it would settle:

Omit empty fields. Do not create a second note for the same result or copy whole
agent reports. Update INDEX/CATALOG only for new paths; RESEARCH_CONTEXT only for
changed programme scope or headline status. Do not rewrite them after every turn.

Before context compaction or a handoff, preserve substantial unsaved findings in
the canonical note. Compact only closed discussion; keep exact active equations,
errors and pending-agent ownership available. Compression is not durable research
memory, and a checkpoint is not evidence that an unperformed check passed.

On resume, read the topic checkpoint and referenced evidence as needed; continue
from the blocker/next step rather than replaying the investigation.
