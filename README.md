# FWB-with-CF: Combining Capacity and Reliability via Inter-Frequency Handovers

This repository contains the code for the paper:

> S. Datta, R. Budhiraja, P. Liu, and S. S. Panwar, "Combining Capacity and Reliability via Inter-Frequency Handovers," *IEEE Journal on Selected Areas in Communications*, vol. 44, pp. 3168–3184, 2026, doi: [10.1109/JSAC.2025.3650542](https://doi.org/10.1109/JSAC.2025.3650542).
>
> ([IEEE Xplore](https://ieeexplore.ieee.org/document/11328784))

We propose offloading high data rate user equipments (UEs), such as those running eXtended Reality (XR) applications, from high-band (HB) frequencies (FR2 and upper FR3) to mid-band (MB) frequencies (FR1 and lower FR3) during HB link outages caused by blockages. Using a cell-free massive MIMO (CF mMIMO) architecture in MB, the proposed inter-frequency handover scheme ensures stable throughput for both low data rate UEs and offloaded high data rate UEs, achieving a hundred-fold reduction in mean outage probability with minimal impact on existing MB users.

## Repository structure

### `Simulation/` — discrete-event simulation of the inter-frequency network

Entry-point scripts (MATLAB, tested with R2023a):

| Script | Purpose |
|---|---|
| `SimulationMain_MIMO.m` | Main discrete-event simulation of HB blockages with offloading to MB CF mMIMO (Algorithm 1, its energy-efficient variant 1E, and exhaustive search). Produces outage duration and outage probability results. |
| `AnalysisMain_MIMO.m` | Analytical model for the outage statistics of the inter-frequency system. |
| `cf_check_limit_density.m` | Motivation experiment: share of Class 2 (high data rate) UEs that MB alone fails to serve; also supports hardware-impairment coefficients. |
| `gen_mmse_cf_stats.m` | Generates CF mMIMO rate statistics with LP-MMSE beamforming (including imperfect CSI and hardware impairments). |
| `try_pei_response_new.m` | Rate comparison of Class 1/Class 2 UEs before and after offloading (impact of power factor and offloading algorithm). |
| `Validate_Result.m` | Validation with ray-tracing channels from the [DeepMIMO](https://deepmimo.net/) 'O1-Blockage' scenario (see [DeepMIMO datasets](#deepmimo-datasets) below). |

Supporting functions include blockage/mobility generation (`Generate_Mobility.m`, `computeBlockageEvents.m`, `lineSegmentIntersect.m`), the discrete-event simulator (`discreteSimulator_MIMO.m` and the `try*/get*/update*` helpers), CF mMIMO channel and rate computation (`generateSetup.m`, `computePhysicalChannels_sub6_MIMO.m`, `compute_link_rates_MIMO*.m`), and truncated-normal UE placement (`TruncatedGaussian.m`).

The `submit*.sbatch` scripts run the corresponding entry points on a SLURM cluster; outputs are written to `Simulation/resultData/`. `dataProcess/combineData.m` aggregates the per-run CSV outputs into summary statistics. Generated figures are collected in `Simulation/plots/`.

#### DeepMIMO datasets

The ray-tracing channel datasets for `Validate_Result.m`, generated from the DeepMIMO 'O1-Blockage' scenario (FR1 at 3.5 GHz and FR2 at 28 GHz), are too large for the git tree and are provided as assets of the [`deepmimo-data` release](https://github.com/sdatta97/FWB-with-CF/releases/tag/deepmimo-data). Download `data-sub6.mat` and `data-mmW.mat` into a `data/` folder at the repository root:

```
FWB-with-CF/
├── data/
│   ├── data-sub6.mat
│   └── data-mmW.mat
├── Simulation/
└── Theory/
```

### `Theory/` — analytical model

- `NumericalChainComputation.m` — numerical computation of the Markov-chain model for HB multi-connectivity and blockage dynamics (uses `State.m`, `bl_rate_calc.m`).
- `CF_outage_analytical.m` — analytical outage probability for the MB CF mMIMO system (uses `hfunc.m`, `jfunc.m`, `term1.m`, `term2.m`).

## Citation

If you use this code in your research, please cite:

```bibtex
@ARTICLE{datta2026combining,
  author={Datta, Soumyadeep and Budhiraja, Rohit and Liu, Pei and Panwar, Shivendra S.},
  journal={IEEE Journal on Selected Areas in Communications},
  title={Combining Capacity and Reliability via Inter-Frequency Handovers},
  year={2026},
  volume={44},
  pages={3168-3184},
  doi={10.1109/JSAC.2025.3650542}
}
```
