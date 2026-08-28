# B1 Stage 2A retry ledger

| Original Slurm job | Method | Failure class | Diagnosis | Corrective action |
|---|---|---|---|---|
| `1404265` | FastOMA 0.5.1 smoke | environment/network | Compute node `cu04` had no DNS access while fetching the official `LUCA.h5` OMAmer database; FastOMA itself did not start. The cluster maps `zenodo.org` to `0.0.0.0`. | Stage the same official Zenodo record under Slurm with TLS hostname verification, enforce its published size (`9918199250`) and MD5 (`87a322a0a6381a8de0b539d62e277054`), record SHA256, then submit a new smoke directory linked to this job. |
| `1404267` | Proteinortho 6.3.6 formal | input compatibility | Proteinortho rejected `*` stop markers in the immutable Open Orthobench FASTA files before similarity search. | Preserve and digest-check the originals; create a run-local copy that removes only `*`, record per-file source/prepared hashes and removal counts, and submit `rep1_retry1`. No scientific parameter changes. |
| `1404268` | Proteinortho 6.3.6 formal retry1 | official-format conversion | Proteinortho completed successfully in 8:10, converter validation passed for all 251,378 proteins, but CRLF line endings left `\r` after each protein ID in the awk aggregation. Python universal-newline parsing therefore read 251,378 singleton prediction lines and official scoring reached a zero denominator. | Keep the completed tool output; strip the trailing `\r` while regenerating only the official-format file, rerun official and extended postprocessing, and use no additional Proteinortho inference run. |
| `1404269` | FastOMA database staging | environment/tool version | Compute-node curl predates `--retry-all-errors`; all eight transfers exited before downloading bytes. | Remove only the unsupported retry flag and resubmit the same published-size/checksum-verified staging workflow. |

Failed directories remain preserved. They are excluded from accuracy and resource tables.
