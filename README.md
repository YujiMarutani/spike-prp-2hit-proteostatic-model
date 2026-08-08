# Mathematical Modeling and Parameter-Space Analysis of the Spike–PrP Neurodegenerative Signaling Hypothesis: A 2-Hit Proteostatic Framework

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-blue.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Epistemic Status: Hypothesis Model](https://img.shields.io/badge/Epistemic_Status-Hypothesis--Driven_Model-orange.svg)](#epistemic-status)

This repository contains the computational framework, ordinary differential equation (ODE) solvers, global sensitivity analysis (LHS-SRCC), two-dimensional parameter-space routines, and full LaTeX source for evaluating the **Spike–PrP 2-Hit Proteostatic Framework**.

---

## 🛡️ Epistemic Status

> **CRITICAL DISCLAIMER**
> 
> * This repository presents a **hypothesis-driven mathematical model**, not an empirical demonstration or clinical prediction.
> * **No Causal Claim**: The model **does not establish** that exposure to the SARS-CoV-2 spike protein causes cellular prion protein ($\text{PrP}^\text{C}$) misfolding, Creutzfeldt–Jakob disease (CJD), or any human prion disease.
> * **Unvalidated Parameter ($k_{\text{conv}}$)**: The key conversion parameter $k_{\text{conv}}$ represents an **experimentally unvalidated theoretical construct** within the modeled spike-exposure context. If $k_{\text{conv}} = 0$, the modeled misfolded compartment remains identically zero ($S_{\text{c}}(t) \equiv 0$) for all time.
> * **Model-Level Results**: Sensitivity metrics and numerical transition boundaries reflect mathematical properties of the formulated ODE system under assumed parameter ranges, not empirical biological effect sizes or verified physiological tipping points.

---

## 📋 Executive Summary

* **Objective**: To translate qualitative, literature-based signaling hypotheses linking spike protein exposure to cellular stress and prion homeostasis into a quantitative 7-dimensional ODE system.
* **Core Finding**: In a healthy proteostatic environment ($A_0 = 1.0$), transient spike exposure generates only a temporary, self-limiting aggregate pulse ($S_{\text{c}} \to 0$). Persistent accumulation within the model requires a **2-Hit state**: transient signaling initiation coupled with a pre-existing, severe degradation of basal autophagic capacity ($A_0 \le 0.25$).
* **Sensitivity Insights**: Aggregate peak height ($\max S_{\text{c}}$) is governed primarily by initial signal surge ($v_1$) and conversion kinetics ($k_{\text{conv}}$), whereas total integrated aggregate residence time ($\text{AUC of } S_{\text{c}}$) is overwhelmingly controlled by host autophagic clearance capacity ($k_{\text{clear}}$, $A_0$).

---

## 📄 Abstract

Recent narrative reviews and theoretical models have hypothesized pathways through which the SARS-CoV-2 spike protein might interact with host signaling cascades—specifically TLR4–MAPK–p53 signaling, autophagic suppression, and molecular mimicry—potentially impacting cellular prion protein ($\text{PrP}^\text{C}$) homeostasis. However, these conceptual frameworks often lack formal dynamic evaluation, leaving unresolved whether proposed interactions represent transient stress responses or self-sustaining pathological cascades.

Here, we construct a 7-dimensional system of ordinary differential equations (ODEs) to formalize this hypothesis and analyze its nonlinear dynamics. Using Latin Hypercube Sampling paired with Spearman’s Rank Correlation Analysis (LHS-SRCC), we evaluate parameter sensitivity across peak aggregate exposure ($\max \text{PrP}^{\text{Sc}}$) and integrated exposure ($\text{AUC of } \text{PrP}^{\text{Sc}}$). Furthermore, two-dimensional phase space analysis demonstrates that transient spike protein exposure alone is insufficient to induce persistent misfolding in a healthy proteostatic environment. Instead, persistent aggregation requires a "2-Hit" mechanism: transient signaling initiation coupled with a pre-existing degradation of basal autophagic capacity ($A_0 \le 0.25$). This model provides a quantitative framework for distinguishing non-causal transient stress from critical parameter-space thresholds in neurodegenerative signaling hypotheses.

**Keywords**: SARS-CoV-2 Spike Protein, Prion Homeostasis, $\text{PrP}^\text{C} \to \text{PrP}^{\text{Sc}}$ Misfolding, Ordinary Differential Equations, Parameter-Space Analysis, Autophagy, 2-Hit Model, Global Sensitivity Analysis, Systems Biology.

---

## ❓ Research Questions

1. **Dynamic Separation**: Which upstream signaling nodes govern the transient peak height versus the integrated residence time of modeled aggregates?
2. **Buffer Capacity**: Can basal host autophagic capacity ($A_0$) completely buffer and eliminate transient aggregate accumulation under hypothetical conversion kinetics?
3. **Transition Thresholds**: Under what parameter combinations does the model transition from complete physiological clearance to persistent aggregate accumulation?

---

## ⚠️ Epistemic Scope & Limitations

* **Hypothetical Coupling**: The model assumes a mathematical coupling between p53 activation, $PRNP$ transcription, and molecular conversion. These links represent theoretical hypotheses derived from literature reviews (Kyriakopoulos et al., 2022; Seneff et al., 2023) and lack direct experimental validation in vivo.
* **Operational Threshold**: The contour at $S_{\text{c}}(t=150) = 0.05$ in Experiment B is an operational numerical boundary used for parameter-space mapping. It must **not** be interpreted as mathematical proof of a formal transcritical bifurcation or a validated clinical threshold.
* **Non-Negativity Constraints**: Numerical state variables are protected via $\max(X, 0)$ operations to maintain stability during integration, preventing unphysical negative concentration values.

---

## {🧮 Mathematical Model

### State Variables ($N = 7$)

| Variable | Biological Definition |
| :--- | :--- |
| $M(t)$ | Activated p38 MAPK / JNK kinase level |
| $D(t)$ | DUSP / Wip1-like phosphatase negative feedback activity |
| $P(t)$ | Nuclear p53 activity index |
| $C(t)$ | Cellular prion protein compartment ($\text{PrP}^\text{C}$) |
| $I(t)$ | Neuroinflammatory / NF-$\kappa$B-associated signaling index |
| $A(t)$ | Autophagic / lysosomal clearance capacity |
| $S_{\text{c}}(t)$ | Modeled misfolded PrP compartment ($\text{PrP}^{\text{Sc}}$) |

### Dynamic Equations

$$\frac{dM}{dt} = v_1 \left( \frac{S(t)}{K_S + S(t)} \right) - k_{\text{dusp}} D M - d_M M$$

$$\frac{dD}{dt} = g_0 - k_{\text{supp}} M D - d_D D$$

$$\frac{dP}{dt} = k_{\text{p53}} M - \gamma_P P$$

$$\frac{dI}{dt} = k_{\text{NF}} M - d_I I$$

$$\frac{dA}{dt} = \frac{A_0}{1 + \left( \frac{\max(I, 0)}{K_I} \right)^m} - d_A A$$

$$\frac{dC}{dt} = V_0 \left( 1 + \beta \frac{P^n}{K_P^n + P^n} \right) - k_{\text{conv}} C \left(\max(S_{\text{c}}, 0)\right)^a - \left( d_C + k_{\text{auto}} A \right) C$$

$$\frac{dS_{\text{c}}}{dt} = k_{\text{conv}} C \left(\max(S_{\text{c}}, 0)\right)^a - k_{\text{clear}} A \max(S_{\text{c}}, 0) - d_{\text{Sc}} S_{\text{c}}$$

### Input Transient Pulse

$$S(t) = \begin{cases} S_{\max} e^{-k_s t}, & t \le t_{\text{spike\_end}} \\ 0, & t > t_{\text{spike\_end}} \end{cases}$$

---

## 📊 Baseline Parameters

| Parameter | Baseline Value | Description |
| :--- | :---: | :--- |
| $S_{\max}$ | $2.0$ | Peak transient input amplitude |
| $k_s$ | $0.1$ | Input decay rate constant |
| $t_{\text{spike\_end}}$ | $20.0$ | Input duration cutoff |
| $v_1$ | $1.0$ | MAPK activation rate |
| $K_S$ | $0.5$ | Half-saturation constant for input $S$ |
| $k_{\text{dusp}}$ | $0.3$ | DUSP inactivation rate on MAPK |
| $d_M$ | $0.2$ | Basal MAPK decay rate |
| $g_0$ | $0.5$ | Basal DUSP synthesis rate |
| $k_{\text{supp}}$ | $0.4$ | MAPK-mediated DUSP suppression rate |
| $d_D$ | $0.2$ | Basal DUSP degradation rate |
| $k_{\text{p53}}$ | $0.5$ | p53 activation rate by MAPK |
| $\gamma_P$ | $0.2$ | Basal p53 decay rate |
| $V_0$ | $0.2$ | Basal $\text{PrP}^\text{C}$ transcription rate |
| $\beta$ | $3.0$ | p53-mediated $\text{PrP}^\text{C}$ induction gain |
| $K_P$ | $0.5$ | Half-saturation constant for p53 induction |
| $n$ | $2.0$ | Hill coefficient for p53 transcription |
| $d_C$ | $0.1$ | Basal $\text{PrP}^\text{C}$ decay rate |
| $k_{\text{auto}}$ | $0.1$ | Autophagic turnover rate of $\text{PrP}^\text{C}$ |
| $k_{\text{NF}}$ | $0.4$ | Inflammatory activation rate by MAPK |
| $d_I$ | $0.2$ | Inflammatory signal decay rate |
| $A_0$ | $1.0$ | Basal autophagic capacity |
| $K_I$ | $0.4$ | Half-saturation constant for autophagic suppression |
| $m$ | $2.0$ | Hill coefficient for inflammatory suppression |
| $d_A$ | $0.2$ | Basal autophagic decay rate |
| $k_{\text{clear}}$ | $0.3$ | Autophagic clearance rate of $S_{\text{c}}$ |
| $d_{\text{Sc}}$ | $0.1$ | Basal clearance/decay rate of $S_{\text{c}}$ |
| $a$ | $1.0$ | Conversion reaction order |
| $k_{\text{conv}}$ | $0.08$ | **Theoretical conversion parameter (experimentally unvalidated)** |

---

## 💻 Environment & Quick Start

### Prerequisites
* Python 3.8 or higher
* Core libraries: `numpy`, `scipy`, `matplotlib`

```bash
pip install numpy scipy matplotlib
