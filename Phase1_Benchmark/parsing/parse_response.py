"""
ONCVIZ-001 Phase 1 -- structured-output parsing rubric.

Converts a raw model response into three fields: a binary flag
(CORRECT/FLAWED), a matched error category (closed-list code, OTHER, or
NONE), and a confidence level. Two parsers are provided:

  parse_strict()      -- for responses to prompt_strict_closed_list.txt,
                          which follow a fixed VERDICT/ERROR_CATEGORY/
                          DESCRIPTION/CONFIDENCE structure. Parsing is
                          exact-field extraction (regex), not inference.

  parse_freeformat()  -- for responses to prompt_freeformat.txt, which have
                          no fixed structure. This is necessarily a
                          heuristic: it infers the verdict from keyword
                          cues in the first ~300 characters of the
                          response, since there is no structured VERDICT
                          field to read. This heuristic was validated by
                          manual spot-check in the original study but is
                          materially less precise than parse_strict() and
                          should not be treated as equally reliable --
                          the study explicitly flags this in its
                          Limitations.

Reconstructed from the documented build session for this benchmark.
"""

import re

VALID_CATEGORIES = {
    "KM01", "KM02", "KM03", "KM04",
    "FP01", "FP02", "FP03", "FP04",
    "WF01", "WF02", "WF03", "WF04",
    "SW01", "SW02", "SW03", "SW04",
    "CN01", "CN02", "CN03", "CN04",
    "NONE", "OTHER",
}


def parse_strict(raw_response: str) -> dict:
    """Extract VERDICT / ERROR_CATEGORY / DESCRIPTION / CONFIDENCE from a
    response to the strict closed-list prompt. Returns None for fields it
    cannot find (do not silently default -- an unparseable field must be
    visible downstream as a scoring exclusion, not a false negative)."""
    out = {"predicted_verdict": None, "predicted_error_id": None,
           "description": None, "confidence": None}

    m = re.search(r"VERDICT:\s*\[?(\w+)\]?", raw_response, re.IGNORECASE)
    if m:
        v = m.group(1).upper()
        if v in ("CORRECT", "FLAWED"):
            out["predicted_verdict"] = v

    m = re.search(r"ERROR_CATEGORY:\s*\[?([A-Za-z0-9]+)", raw_response, re.IGNORECASE)
    if m:
        code = m.group(1).upper()
        # Accept codes even if the model wrote a trailing description after
        # the code (e.g. "KM01 - non-monotonic survival curve"), which
        # happened under strict format for GPT-5.6 -- this is a real,
        # documented format-compliance failure mode, not a parsing choice
        # to paper over: log it, don't silently "fix" it into an exact
        # match if you are trying to reproduce the original exact-match
        # scoring exactly.
        if code in VALID_CATEGORIES:
            out["predicted_error_id"] = code
        else:
            out["predicted_error_id"] = "OTHER"

    m = re.search(r"DESCRIPTION:\s*(.+?)(?:\nCONFIDENCE:|\Z)", raw_response,
                   re.IGNORECASE | re.DOTALL)
    if m:
        out["description"] = m.group(1).strip()

    m = re.search(r"CONFIDENCE:\s*\[?(\w+)\]?", raw_response, re.IGNORECASE)
    if m:
        c = m.group(1).upper()
        if c in ("LOW", "MEDIUM", "HIGH"):
            out["confidence"] = c

    return out


# ---------------------------------------------------------------------------
# Free-format heuristic parser
# ---------------------------------------------------------------------------

_FLAWED_HINTS = [
    "error", "flaw", "incorrect", "wrong", "issue", "inconsist",
    "does not", "doesn'?t", "missing", "mismatch", "violat",
]
_CORRECT_HINTS = [
    "no error", "no issue", "looks correct", "is correct",
    "appears correct", "no problem", "consistent with", "no flaw",
]


def infer_verdict_freeformat(text: str) -> str:
    """Heuristically infer CORRECT / FLAWED / UNCLEAR from the first ~300
    characters of a free-format response. This is a documented, imperfect
    heuristic -- validate the UNCLEAR rate on any new batch before trusting
    the aggregate numbers (the original study reported this breakdown
    explicitly rather than treating it as fully reliable)."""
    head = text.lower()[:300]
    if any(re.search(p, head) for p in _CORRECT_HINTS):
        return "CORRECT"
    if any(re.search(p, head) for p in _FLAWED_HINTS):
        return "FLAWED"
    return "UNCLEAR"


def parse_freeformat(raw_response: str) -> dict:
    """Returns {'predicted_verdict': ..., 'raw_description': ...}. There is
    no structured category or confidence field to extract in free-format
    mode -- category matching for these responses is done separately via
    the semantic keyword taxonomy (see semantic_scoring/), not by this
    parser."""
    return {
        "predicted_verdict": infer_verdict_freeformat(raw_response),
        "raw_description": raw_response,
    }
