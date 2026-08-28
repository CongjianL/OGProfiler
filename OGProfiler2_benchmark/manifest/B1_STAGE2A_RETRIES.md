# B1 Stage 2A retry ledger

| Original Slurm job | Method | Failure class | Diagnosis | Corrective action |
|---|---|---|---|---|
| `1404265` | FastOMA 0.5.1 smoke | environment/network | Compute node `cu04` had no DNS access while fetching the official `LUCA.h5` OMAmer database; FastOMA itself did not start. | Stage the same official database outside the compute step, record its SHA256, then submit a new smoke directory linked to this job. |
| `1404267` | Proteinortho 6.3.6 formal | input compatibility | Proteinortho rejected `*` stop markers in the immutable Open Orthobench FASTA files before similarity search. | Preserve and digest-check the originals; create a run-local copy that removes only `*`, record per-file source/prepared hashes and removal counts, and submit `rep1_retry1`. No scientific parameter changes. |

Failed directories remain preserved. They are excluded from accuracy and resource tables.
