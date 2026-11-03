# Source data and controlled inputs

Figure-specific aggregate numerical tables are exported separately with `tools/build_source_data.py` after the analysis completes. They are not bundled in the code repository. Before release, the institutional data custodian must approve aggregate outputs and reference-model parameters derived from restricted study data.

The full workflow uses six locally staged files under the untracked `rebuild/data/` directory:

| File | Contents |
|---|---|
| `nested.csv` | Discovery sample metadata |
| `m_tr.csv` | Discovery methylation measurements |
| `longitudinal.csv` | Paired longitudinal sample metadata |
| `m_te.csv` | Longitudinal methylation measurements |
| `ph-ext.csv` | External cohort metadata |
| `m-ext.csv` | External methylation measurements |

Individual-level data, sample identifiers, exact collection dates, posterior memberships and clinical metadata are not distributed through this repository. Access to restricted data depends on the originating institutions and study agreements. The external PBMC methylation dataset is publicly available as NCBI GEO accession GSE275371. It should be obtained from the original record under the applicable terms.

The full pipeline additionally requires the recurrent-program discovery checkpoint under `rebuild/cache/discovery/`, or the option `--recompute-discovery` with sufficient time and computing resources.
