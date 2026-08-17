# ONCVIZ-001 Phase 1 Benchmark -- Release Artifacts

Supporting artifacts for "Human-Calibrated Validation of LLM Oncology
Clinical Trial Visualizations." This is the Phase 1 calibration benchmark:
5 of 12 planned plot types (Kaplan-Meier, waterfall, forest, swimmer,
CONSORT), 20 of 46 planned taxonomy errors, 25 samples (5 reference + 20
flawed), evaluated on Claude Opus 4.8 and GPT-5.6.

## IMPORTANT: provenance of these files (read before using)

This directory was reconstructed in a later working session than the one
in which the benchmark was originally built. Two categories of file here
have different reliability:

- **Exact, verified-real data**: `raw_outputs/combined_final_strict_prompt.csv`
  is the actual 228-response output file, byte-identical to what the
  original test run produced. `semantic_scoring/semantic_score.py` and
  `scoring/score.py` / `scoring/score_unique_sample.py` were re-run
  against this real file in this session and confirmed to reproduce the
  documented headline numbers (Claude 97.0% semantic / 100% exact; GPT
  76.8% semantic / 39.4% exact; unique-sample specificity Claude 100%,
  GPT 40.0%).
- **Reconstructed from documented decisions, not re-verified against
  original code**: the injector scripts (`injectors/`), checker functions
  (`checkers/checkers.py` -- except the 3 functions marked "verbatim"),
  and the parsing rubric (`parsing/parse_response.py`) are rebuilt to the
  exact logical specification recorded during the original build (which
  was documented in detail, including specific bugs found and fixed along
  the way), but the original R/Python files themselves no longer exist on
  any accessible system. Treat these as a faithful specification, and
  re-test against your own rendering pipeline before relying on them for
  a new run.

**Not available anywhere, at any point:** token usage and API cost per
response. This was never logged during the original data collection.
`raw_outputs/combined_final_strict_prompt.csv` has no token/cost columns
because none were ever recorded -- this is not a redaction.

**Not included in this bundle:** raw free-format responses for GPT-5.6 and
Claude Opus 4.8 (the supplementary tests referenced in `prompts/`). These
exist in the original study but were only available in this later session
as content pasted into a chat transcript, not as an exact file re-upload
the way `combined_final_strict_prompt.csv` was -- ask if you want a
best-effort reconstruction of those flagged the same way.

## Directory map

```
taxonomy/
  taxonomy.csv, taxonomy.json    20-error taxonomy: code, plot type, severity,
                                  governing standard, injection method,
                                  injector/checker function names

injectors/
  injectors_kaplan_meier.R       R functions/patches that inject each error
  injectors_waterfall.R          into a reference plot's data or script.
  injectors_forest.R             data_injection errors are R functions;
  injectors_swimmer.R            code_patch errors are documented as a
  injectors_consort.R            diff against the production script (the
                                  real scripts have no configurable "flag"
                                  hooks, so a patched copy is required).

checkers/
  checkers.py                    Ground-truth checker functions -- one per
                                  error, programmatically verifying
                                  presence/absence independent of any
                                  model judgment. 3 are verbatim from the
                                  original session; the rest are
                                  reconstructed to the same specification.

prompts/
  prompt_strict_closed_list.txt  Exact prompt for the main test (fixed
                                  VERDICT/ERROR_CATEGORY/DESCRIPTION/
                                  CONFIDENCE structure, closed-list code).
  prompt_freeformat.txt          Exact prompt for the supplementary
                                  open-ended test (no fixed structure).

parsing/
  parse_response.py              parse_strict() -- exact field extraction
                                  from the structured prompt's output.
                                  parse_freeformat() -- heuristic verdict
                                  inference for the unstructured prompt's
                                  output (documented as less reliable).

raw_outputs/
  combined_final_strict_prompt.csv   The real 228-response file: both
                                      models, all 25 samples, 3 repeats,
                                      up to 2 input conditions each.

semantic_scoring/
  semantic_score.py              Bug-fixed keyword-based semantic scorer.
                                  Includes the polarity-check fix and a
                                  documented account of the bug it fixed
                                  (an earlier version over-credited GPT-5.6
                                  by matching keywords regardless of
                                  negation).

scoring/
  score.py                       Primary metrics: sensitivity, specificity,
                                  critical-error FNR, category-match
                                  accuracy. Repeat-level (each of 3 repeats
                                  is an independent observation).
  score_unique_sample.py         Same metrics at the unique-sample level
                                  (majority vote across repeats), with the
                                  majority-vote rule explicitly justified
                                  in the file's docstring.
```

## How the pieces connect (pipeline order)

1. `taxonomy/` defines what error goes with what plot type, severity, and
   injection method.
2. `injectors/` implements each taxonomy entry's injection against the
   reference catalog's real production R scripts (external to this
   directory -- the catalog itself is a separate repository component).
3. `checkers/` independently verifies each injected error actually
   produced the intended defect (used to validate injectors during
   construction, not used at model-evaluation time).
4. `prompts/` were shown to each model alongside each rendered image (and
   optionally its underlying code, per `input_condition`).
5. `parsing/` converts each raw model response into structured fields.
6. The parsed results, pooled across both models/all samples/all repeats/
   both conditions, are `raw_outputs/combined_final_strict_prompt.csv`.
7. `scoring/score.py` and `scoring/score_unique_sample.py` compute
   detection metrics from that file.
8. `semantic_scoring/semantic_score.py` computes a second, independent
   metric (semantic recognition) from the same file's free-text
   `raw_description` column, to separate "did the model understand the
   error" from "did the model output the exact required code."

## Known open items (carried over honestly, not resolved here)

- Statistical tests referenced in early planning documents (weighted
  Cohen's kappa, Fleiss' kappa, McNemar's test) were never actually
  computed against this data -- only confusion-matrix metrics and Wilson
  95% CIs were. If a future version of the manuscript references
  kappa/McNemar's as performed, verify against this fact before trusting it.
- The forest-plot reference image (`FP_reference`) had a real,
  undisclosed-at-the-time rendering defect (CI-clipping bug) during the
  strict-prompt test. Re-scoring with this corrected in the ground truth
  changes specificity for the 4 remaining clean reference images to 100%
  (Claude) and 50% (GPT-5.6), and reveals one additional false negative
  for Claude (it did not catch the defect; GPT-5.6 did, under the
  free-format condition where the bug was actually discovered).
- Only 20 of the fully-planned 46-error taxonomy exist as itemized,
  implemented entries. The remaining 26 (across 7 additional plot types)
  are unscoped beyond being named in the plot-type catalog.
