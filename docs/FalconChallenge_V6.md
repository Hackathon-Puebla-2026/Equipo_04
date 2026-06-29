# Guided Quantum Computing Challenges for Transboundary Water Systems

**Reservoir Operations and Hydrological Monitoring in the Mexico–United States Rio Grande / Río Bravo Basin**

Supporting SDG 6.4, SDG 6.5, and SDG 13.1

*June 24, 2026*

> Markdown transcription of `FalconChallenge/FalconChallenge_V6.pdf` (the authoritative spec). Provided so it can be loaded into context cheaply and grepped. If anything here disagrees with the PDF, the PDF wins.

---

# Guided Challenge A: Resilient Release Scheduling for the International Falcon Reservoir

## 1. Specific SDG Sub-Target

This challenge is primarily aligned with **SDG 6: Clean Water and Sanitation**, with emphasis on the following targets:

- **Target 6.4:** By 2030, substantially increase water-use efficiency across all sectors and ensure sustainable withdrawals and supply of freshwater to address water scarcity.
- **Target 6.5:** By 2030, implement integrated water resources management at all levels, including through transboundary cooperation as appropriate.

The challenge also contributes to **SDG 13: Climate Action**, particularly:

- **Target 13.1:** Strengthen resilience and adaptive capacity to climate-related hazards and natural disasters in all countries.

The connection to **SDG 6.4** comes from the need to manage limited freshwater availability under scarcity conditions. The connection to **SDG 6.5** is direct because Falcon Reservoir is part of the Mexico–United States Rio Grande / Río Bravo binational water-management system.

The connection to **SDG 13.1** arises because reservoir operation is affected by climate-related hydrological stressors, including droughts, prolonged low-inflow periods, increased evaporation, and extreme precipitation events. Improving release scheduling can therefore support resilience and adaptive capacity under climate-related water-availability variability.

The concrete question is whether data-driven optimization can help identify release schedules that improve reservoir resilience during low-storage periods. The measurable quantity of interest is the **Storage Resilience Score (SRS)**, defined as follows:

$$\text{SRS} = -\left(w_1 C_{\text{crit}} + w_2 C_{\text{dev}} + w_3 C_{\text{smooth}}\right) \tag{1}$$

The weight $w_1$ gives priority to avoiding critically low storage. The weights $w_2$ and $w_3$ penalize unrealistic deviations from historical operation, including large release changes and abrupt week-to-week variations.

Let $S_{\text{opt}}(t)$ be the optimized reservoir storage, $R^{\text{obs}}(t)$ the observed historical release, and $u(t)$ the release adjustment. The optimized release is

$$R(t) = R^{\text{obs}}(t) + u(t), \tag{2}$$

$u(t) = 0$ reproduces historical operation. Note that this formulation searches for $u(t)$. Using the observed storage change $\Delta S_{\text{obs}}(t)$ from `Discharge.Total.Change-in-Storage@08461200`, the simplified storage dynamics are

$$S_{\text{opt}}(t+1) = S_{\text{opt}}(t) + \Delta S_{\text{obs}}(t) - u(t). \tag{3}$$

Participants should use the following definitions:

$$C_{\text{crit}} = \sum_{t=0}^{T} \left[\max\left(0,\, S_{\text{min}} - S_{\text{opt}}(t)\right)\right]^2,$$

$$C_{\text{dev}} = \sum_{t=0}^{T-1} u(t)^2,$$

$$C_{\text{smooth}} = \sum_{t=1}^{T-1} \left[u(t) - u(t-1)\right]^2.$$

Alternative normalizations or weights may be explored, but they must be reported separately from the official benchmark score.

## 2. Identify the Computational Bottleneck

Falcon Reservoir operation involves time-dependent decisions: at each time step, a release level or release adjustment must be selected. Over many time steps, the number of possible release schedules grows combinatorially. The computational bottleneck is therefore the search for a feasible release-adjustment sequence $u(t)$ that improves storage resilience while respecting operational feasibility.

The optimization variable is the release adjustment $u(t)$, while the release itself is derived as

$$R(t) = R^{\text{obs}}(t) + u(t). \tag{4}$$

The feasible release-adjustment sequence must satisfy the following restrictions:

$$R(t) \geq 0, \tag{5}$$

$$|u(t)| \leq u_{\text{max}}, \tag{6}$$

$$0 \leq S_{\text{opt}}(t) \leq S_{\text{max}}, \tag{7}$$

and

$$\left|\sum_{t=0}^{T-1} u(t)\right| \leq \eta \sum_{t=0}^{T-1} R^{\text{obs}}(t). \tag{8}$$

The last constraint prevents solutions that improve storage only by systematically reducing total releases over the full time window. It forces the optimizer to redistribute releases in time rather than simply withholding water.

If $L$ release-adjustment levels are allowed at each of $T$ time steps, the number of candidate schedules scales as $L^T$.

## 3. Baseline and Official Benchmark Setup

Classical reservoir-operation problems are commonly formulated as time-dependent control problems, where release decisions are selected over a planning horizon subject to storage, feasibility, and operational constraints [1] [2] [3]. This challenge uses a simplified benchmark version of that setting. The official benchmark setup below defines the common parameters used to evaluate both classical and quantum or hybrid solutions.

For the official benchmark, participants should use $L = 5$ release-adjustment levels:

$$u(t) \in \{-2\Delta u,\, -\Delta u,\, 0,\, \Delta u,\, 2\Delta u\}. \tag{9}$$

The adjustment scale is defined from the observed release record as

$$\Delta u = 0.25\, \widetilde{R}^{\text{obs}}_{\text{week}}, \tag{10}$$

where $\widetilde{R}^{\text{obs}}_{\text{week}}$ is the median weekly observed release volume over the selected time window. Therefore,

$$u_{\text{max}} = 2\Delta u. \tag{11}$$

The official SRS weights are

$$w_1 = \frac{1}{(T+1) S_{\text{scale}}^2}, \quad w_2 = \frac{0.1}{T u_{\text{max}}^2}, \quad w_3 = \frac{0.1}{(T-1)(2 u_{\text{max}})^2}, \tag{12}$$

with

$$S_{\text{scale}} = S_{\text{min}}. \tag{13}$$

The cumulative release-balance constraint should use

$$\eta = 0.10. \tag{14}$$

For the official benchmark, $S_{\text{max}}$ is taken as the Falcon total conservation storage capacity reported in the reservoir overview, and the critical storage threshold is set to

$$S_{\text{min}} = 0.25\, S_{\text{max}}. \tag{15}$$

Alternative values of $L$, $\Delta u$, $w_i$, or $\eta$ may be explored as sensitivity analyses, but should be reported separately from the official benchmark score.

**Historical replay reference.** The historical replay reference uses the observed release record `Discharge.Best Available@08461300` without optimization. In this case,

$$u_{\text{hist}}(t) = 0, \tag{16}$$

and therefore

$$R_{\text{hist}}(t) = R^{\text{obs}}(t). \tag{17}$$

The corresponding storage-resilience score is denoted by $\text{SRS}_{\text{hist}}$.

**Classical baseline.** As a simple classical baseline, participants should evaluate a threshold-based conservation rule using the same official value of $\Delta u$:

$$u_{\text{rule}}(t) = \begin{cases} -\Delta u, & S_{\text{rule}}(t) < S_{\text{min}}, \\ 0, & S_{\text{rule}}(t) \geq S_{\text{min}}. \end{cases} \tag{18}$$

The resulting release is

$$R_{\text{rule}}(t) = R^{\text{obs}}(t) + u_{\text{rule}}(t). \tag{19}$$

This provides a reproducible classical comparison that reduces releases during low-storage periods while remaining close to historical operation.

**Optional literature-informed baseline.** Participants may also implement a stronger classical optimization baseline, such as dynamic programming, mixed-integer optimization, or an evolutionary algorithm, using the same SRS, constraints, dataset, and release-adjustment levels. The quantum or hybrid implementation should be compared against the selected classical baseline using the same official benchmark setup.

## 4. Quantum Implementation

The recommended quantum approaches are:

- Quantum Approximate Optimization Algorithm (QAOA),
- quantum annealing,
- hybrid quantum-classical optimization.

Participants should:

- discretize the possible release adjustments,
- encode the release adjustments them as binary variables, and
- justify why the resulting release-scheduling problem has a structure appropriate for quantum or quantum-inspired optimization.

The implementation should be self-contained, reproducible, and benchmarked against the selected classical baseline.

**What remains open to participants.** The benchmark setup defines the data, variables, constraints, and scoring metric needed for fair comparison. Participants are still responsible for choosing the discretization strategy for larger instances, deriving a QUBO or Ising formulation, selecting the quantum or hybrid optimizer, enforcing constraints, tuning algorithmic parameters, and analyzing scaling. Thus, the benchmark is specified, but the computational solution is not provided.

## 5. Benchmarked Solution

The benchmarked solution is the release-adjustment policy obtained with a quantum, quantum-inspired, or hybrid optimization method. Participants should compute an optimized adjustment sequence $u_{\text{opt}}(t)$ and compare it against the classical baseline using the same dataset, time window, storage model, and Storage Resilience Score (SRS).

The optimized release is

$$R_{\text{opt}}(t) = R^{\text{obs}}(t) + u_{\text{opt}}(t). \tag{20}$$

The main benchmark quantity is

$$\Delta \text{SRS} = \text{SRS}_{\text{opt}} - \text{SRS}_{\text{baseline}}. \tag{21}$$

Participants should report:

- $\Delta \text{SRS}$ relative to the selected classical baseline;
- runtime of the classical and quantum or hybrid approaches;
- scaling behaviour as the number of time steps and adjustment levels increases.

## 6. Scaling instances

Participants should evaluate scalability by increasing the number of time steps $T$ and, when feasible, the number of release-adjustment levels $L$.

- **Small instance:** $T = 12$ weeks, $L = 3$ adjustment levels.
- **Medium instance:** $T = 26$ weeks, $L = 5$ adjustment levels.
- **Large instance:** $T = 52$ weeks, $L = 5$ or $L = 7$ adjustment levels.

For $L$ adjustment levels and $T$ time steps, the number of candidate release schedules scales as

$$N_{\text{schedules}} = L^T. \tag{22}$$

The small instance may be used for validation and debugging. The medium instance defines the official benchmark comparison. The large instance is used for scaling analysis and optional advanced comparisons (as sensitivity analyses, but should be reported separately from the official benchmark score).

## 7. Current Quantum Limitations and Scope

Current quantum hardware is not expected to outperform mature classical reservoir-operation or water-resource optimization methods for realistic large-scale instances. The purpose of this challenge is not to demonstrate quantum advantage, but to evaluate whether a simplified reservoir release-scheduling problem can be meaningfully encoded, benchmarked, and scaled using quantum or hybrid optimization methods.

The reservoir model used in this challenge is intentionally simplified. It does not attempt to produce an official operational policy for Falcon Reservoir, nor does it represent the full legal, treaty-based, hydrological, agricultural, municipal, or political complexity of the Rio Grande / Río Bravo system. The goal is to provide a tractable benchmark instance where classical and quantum approaches can be compared on the same storage-resilience metric.

## 8. Data Availability

Recommended public data sources from the IBWC (<https://www.ibwc.gov/water-data/>, particularly <https://waterdata.ibwc.gov/AQWebportal/Data/Dashboard/5>) data portal include:

- **`08461200` — International Falcon Reservoir**
  - `Total Storage.Web-Daily-tcm@08461200`
  - `Total Storage.Web-Telemetry-tcm@08461200`
  - `Reservoir Elevation.Web-Daily-m@08461200`
  - `Lake Area.Best Available@08461200`
  - `Evaporation,accumltd.Daily Evaporation - mm@08461200`
  - `Discharge.Total.Change-in-Storage@08461200`
  - `Percentage.Conservation-Web-Telemetry@08461200`
- **`08461300` — Rio Grande Below Falcon Dam**
  - `Discharge.Best Available@08461300`

The last dataset is used as the observed historical release from Falcon Dam. The weekly IBWC ownership reports may also be used as contextual information because they include Falcon storage, requested releases, total releases, inflows, losses, and binational ownership information.

## 9. Optional Scaling Data

For the official benchmark, participants will be provided with the full downloaded dataset in the shared challenge folder:

```
https://cicesemx0-my.sharepoint.com/:f:/g/personal/fadomin_cicese_mx/IgDamhhYBeWfSb73hdPj-EmXAWlRxnqXUyBlI0u9GUiJvrM?e=YDLUL0
```

This dataset contains the time series required to reproduce the historical baseline, compute the Storage Resilience Score (SRS), and compare quantum, quantum-inspired, hybrid, and classical optimization methods on the same benchmark instance.

Participants are not required to download data directly from the IBWC portal during the hackathon. However, teams interested in extending the scaling analysis may access additional public data from the USIBWC Water Data Portal:

```
https://www.ibwc.gov/water-data/
```

A useful starting point is the section titled *Daily Rio Grande Flow Conditions*. This section includes links to:

- *Rio Grande Flow / Reservoir Conditions*;
- *Stored Water in Large Reservoirs of the Rio Grande Basin*;
- *Stored Water in Large Reservoirs of the Colorado Basin*;
- *Some Mexican Reservoirs*.

For this challenge, the relevant reservoir entry is `FALCON -- 08461200`, located under *Stored Water in Large Reservoirs of the Rio Grande Basin*. Within the USIBWC dashboard, the *Summary* tab lists the available datasets for the selected station, including storage, storage change, deduced inflow, and related reservoir variables. The *Dashboard* and *Welcome* pages provide the navigation tools needed to explore additional parameters and time windows.

For release data, participants should use station `08461300 -- Rio Grande Below Falcon Dam`, specifically the dataset `Discharge.Best Available@08461300`. If discharge is reported in m³/s, it should be converted to a volume over the selected time step before being used in the storage-balance model.

For optional scaling studies, participants may extend the official benchmark in three ways:

- increasing the number of time steps $T$ by using longer Falcon Reservoir time windows;
- increasing the number of release-adjustment levels $L$;
- incorporating additional operation-related information from the Rio Grande reservoir overview or weekly ownership reports, such as ownership, requested releases, total releases, inflows, and losses.

These optional extensions should not replace the official benchmark score. They extend the proposed optimization method to handle more realistic operational constraints.

## References

[1] Giuliani, M., Lamontagne, J. R., Reed, P. M., & Castelletti, A., "A State-of-the-Art Review of Optimal Reservoir Control for Managing Conflicting Demands in a Changing World," *Water Resources Research*, vol. 57, e2021WR029927, 2021. doi: 10.1029/2021WR029927.

[2] Sun, X., Luo, J., & Xie, J., "Multi-Objective Optimization for Reservoir Operation Considering Water Diversion and Power Generation Objectives," *Water*, vol. 10, no. 11, p. 1540, 2018. doi: 10.3390/w10111540.

[3] Zhao, T., & Zhao, J., "Optimizing Operation of Water Supply Reservoir: The Role of Constraints," *Mathematical Problems in Engineering*, vol. 2014, Article ID 853186, 15 pages, 2014. doi: 10.1155/2014/853186.
