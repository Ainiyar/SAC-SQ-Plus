# R382 pairwise-only ablation results

The public runner was re-executed from the shared R392 frozen inputs. The six statuses and created-label counts match the archived R382 experiment. For the three strict-reference arms that complete, the R354 exact-support outputs match the R392 reference exactly.

| Arm | Pairwise-only status | Created labels | Sink candidates | Exact possible paths after R354 | R392 support match |
|---|---:|---:|---:|---:|---|
| Occupancy lookahead | LABEL_LIMIT | 2,000,000 | 0 | — | — |
| Occupancy CROPS+PELT | COMPLETE | 56,422 | 5 | 5 | Yes |
| EEG lookahead | LABEL_LIMIT | 2,000,000 | 0 | — | — |
| EEG CROPS+PELT | COMPLETE | 184,502 | 15 | 15 | Yes |
| SECOM lookahead | LABEL_LIMIT | 2,000,000 | 0 | — | — |
| SECOM CROPS+PELT | COMPLETE | 213,211 | 21 | 18 | Yes |

The three strict arms demonstrate why sink candidates are not interchangeable with possible complete paths: SECOM has 21 sink candidates but only 18 survive exact support filtering. The three pairwise-only SAC arms exhaust the 2M-label budget without completing. The corresponding multi-witness comparison and historical CPU figures are retained in `../R392-R393/figure_data/R382_SIX_ARM_PAIRWISE_VS_MULTI.csv` and `results/original_r382_20261109/`.

The current rerun CPU values are host-specific and are recorded in `results/summary.json`; they are not substituted for the original same-host values in the manuscript.
