# Expert Review: gpa-disease-risk-query v0.2

**Reviewer persona:** `gpa-genomic-phenotype` / `dgra-prefilter` maintainer  
**Review date:** 2026-06-13  
**Scope:** LiftOver refactor, VEP 115 compatibility shim, end-to-end integration.

## 1. Overall impression

The skill cleanly separates concerns:
- `liftover.py` owns genome-build normalization.
- `gpa_compat.py` isolates the VEP-version mismatch workaround.
- `pipeline.py` orchestrates the 6-step PRD flow without embedding GPA internals.

The decision to **pre-annotate the filtered VCF with VEP 115-compatible flags and then pass the annotated VCF to GPA** is the least-invasive fix for the GPA/VEP parameter mismatch. It avoids forking GPA code and keeps all tier-classification logic inside `gpa-genomic-phenotype`.

## 2. What was done well

1. **No LLM-based variant interpretation.** The pipeline consistently delegates variant grading to GPA; the skill only adds disease-context filtering and scoring.
2. **Local-first resource usage.** OMIM SQLite, HPO gene-phenotype file, dgra-prefilter BEDs, GRCh38 FASTA, and VEP Docker cache are all reused rather than re-downloaded.
3. **Defensive error handling.**
   - REF validation after liftover aborts if mismatch rate > 5%.
   - `pre_annotate_with_vep` fails loudly when the VEP Docker image is missing.
   - GPA result is propagated as an error dict rather than raising unhandled exceptions.
4. **Parser monkey-patch is narrowly scoped.** Only `VCFParser._csq_to_variant` is patched, and only for the missing `gnomAD_AF` fallback.
5. **Preflight is comprehensive.** It now checks `pyliftover`, VEP Docker image, VEP cache, and chain file availability.

## 3. Concerns and recommendations

### 3.1 LiftOver REF/ALT correction for indels is conservative

**Concern:** The current implementation drops indels when the target reference base does not match the source REF. This is safe but may lose real variants in patch regions.

**Recommendation (P2):** For future versions, consider running `bcftools norm -f target.fa` after coordinate remapping, or implement a prefix-alignment routine for small indels. For the MVP, document this limitation clearly.

### 3.2 VEP pre-annotation writes to `filtered_vcf.parent`

**Concern:** `gpa_runner.py` writes `*.vep115.vcf.gz` next to the filtered VCF. In a multi-sample or long-running workflow this could clutter the work directory.

**Recommendation (P2):** Place the pre-annotated VCF in a temporary subdirectory by default, with an optional `--keep-annotated-vcf` flag for debugging.

### 3.3 GPA `GeneListSynchronizer` warning is harmless but noisy

**Observation:** GPA prints `RuntimeWarning: coroutine 'GeneListSynchronizer.get_merged_gene_lists' was never awaited`. This is a known GPA v0.10.5 warning when called inside an existing event loop; it does not affect results.

**Recommendation (P3):** Suppress or filter this specific warning in `gpa_runner.py` once the behavior is confirmed stable.

### 3.4 Tier classification still defaults to Tier 3 for some ClinVar pathogenic missense variants

**Observation:** In the end-to-end test, real ClinVar Pathogenic missense variants in APP/PSEN1 were classified as Tier 3 by GPA because of low quality confidence (missing VAF before AD was added) and gnomAD query failures. After adding AD/VAF, classification remained Tier 3 due to GPA's internal weighting.

**Recommendation (P2):** This is GPA behavior, not a bug in this skill. However, the report should explicitly surface ClinVar pathogenic status and the `upgrade_conditions` so users do not misinterpret Tier 3 as benign. The current report already includes ClinVar info, which is good.

### 3.5 Test coverage for chromosome naming edge cases

**Concern:** The liftover test uses chr-prefixed contigs. Non-chr input (e.g., Ensembl-style `1`) is handled by the code but not exercised in tests.

**Recommendation (P3):** Add a unit test with Ensembl-style chromosome names to prevent regressions.

## 4. Action items taken during review

- [x] Replaced broken CrossMap path with `pyliftover` + `vcfpy` + `samtools faidx`.
- [x] Added FASTA chromosome-naming normalization (`chr` vs no-`chr`).
- [x] Added `gpa_compat.py` with VEP 115-compatible Docker flags and `gnomADe_AF`/`gnomADg_AF` parser fallback.
- [x] Updated `preflight_check` to verify `pyliftover`, VEP Docker image, VEP cache, and liftover chain.
- [x] Added unit tests for liftover and VEP 115 compat annotation.
- [x] Documented the VEP 115 compatibility layer in `SKILL.md`.

## 5. Verdict

**Approved with P2/P3 suggestions.** The skill is functional with the local toolchain and correctly delegates interpretation to GPA. The remaining recommendations are improvements rather than blockers.
