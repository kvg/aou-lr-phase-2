# FLARE model fixtures

Synthetic full-format (2A+5) models used when GCS download is unavailable.

Intended real URIs (AoU storage workspace):

- `gs://fc-secure-8f7d6a20-.../lai_exp.em_all_pops.AFR.model`
- `gs://fc-secure-8f7d6a20-.../lai_exp.pin_all_pops_gen_by_pop.AFR.model`

Replace these files with real downloads when credentials are available;
`test_flare_model.py` round-trips whatever is present.
