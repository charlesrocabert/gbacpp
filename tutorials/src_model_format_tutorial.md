<p align="center">
  <img src="https://github.com/user-attachments/assets/cb2ea68d-7cb8-4094-8e69-d0c84af8bbbf" width=200 />

</p>

<p align="center">
  <h3 align="center">Self-replicating cell (SRC) model format tutorial</h3>
</p>

-----------------

Self-replicating cell (SRC) models must comply to a standard format to be compatible with growth balance analysis (GBA) formalism. Please refer to <a href="https://doi.org/10.1371/journal.pcbi.1011156" target="_blank">Dourado et al. (2023)</a> and the tutorial available on https://cellgrowthsim.com/ for other detailed sources.

Two formats are used to distribute SRC models: the text format CSV and the OpenDocument Spreadsheet format ODS. Data organization is strictly identical in both formats. As for now, <strong>gbacpp</strong> reads SRC models in CSV format only.

# Table of contents
- [1) Files organization](#files_organization)
- [2) Files content](#files_content)
  - [2.1) Model information (<code>Info.csv</code>)](#info)
  - [2.2) Mass fraction matrix (<code>M.csv</code>)](#M)
  - [2.3) Forward and backward turnover rates vectors (<code>kcat.csv</code>)](#kcat)
  - [2.4) Michaelis constants matrix (<code>K.csv</code>)](#K)
  - [2.5) External conditions matrix (<code>conditions.csv</code>)](#conditions)
  - [2.6) Initial solution (<code>f0.csv</code>)](#f0)
  - [2.7) Activation constants matrix (<code>KA.csv</code>)](#KA)
  - [2.8) Inhibition constants matrix (<code>KI.csv</code>)](#KI)
  - [2.9) List of constant reactions (<code>constant_reactions.csv</code>)](#constant_reactions)
  - [2.10) Enzyme to protein mass concentration mapping (<code>protein_contributions.csv</code>)](#protein_contributions)

# 1) Files organization <a name="files_organization"></a>

The SRC model is organized as a set of CSV files located in a folder having the name of the model. They usually contain GBA variables (such as matrices, or vectors), but also additional variables and information.

      └── model_name
           ├── Info.csv
           ├── M.csv
           ├── kcat.csv
           ├── K.csv
           ├── conditions.csv
           ├── f0.csv
           ├── KA.csv
           ├── KI.csv
           ├── constant_rhs.csv
           ├── constant_reactions.csv
           └── protein_contributions.csv

Some files are mandatory to set a minimal SRC model. These are:
- `M.csv`
- `kcat.csv`
- `K.csv`
- `conditions.csv`
- `f0.csv`

Other files are optional.

⚠️ All CSV files have `;` separators, and `\n` line breaks.

# 2) Files content <a name="files_content"></a>

### 2.1) Model information (<code>Info.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="info"></a>

The optional file `Infos.csv` contains various information about the model (units, file description, ...).

The content is free but must follow a standard structure:

      └── Category
           ├── Key 1: Descriptor 1
           ├── Key 2: Descriptor 2
           └── ...

Usually, categories are "General" (model name, short description, ...), "Units", "Sheets", etc...

For example, the toy model A has the following information:

      └── General
           ├── Name: A
           └── Description: Simplest model with two reactions
      └── Units
           ├── KM: g/L
           ├── kcat: 1/h ([mass of products]/[mass of protein]/h)
           └── rho: g/L
      └── Sheets
           ├── M: Mass fraction matrix
           ├── K: Forward Michaelis constant matrix
           ├── kcat: Turnover numbers
           └── conditions: Value of rho and external concentrations at different growth conditions

### 2.2) Mass fraction matrix $\mathbf{M}$ (<code>M.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="M"></a>

The file `M.csv` contains the mass fraction matrix $\mathbf{M}$, which is the pendant of the stoichiometric matrix in normalized mass units (see <a href="https://doi.org/10.1371/journal.pcbi.1011156" target="_blank">Dourado et al., 2023</a>). This file is mandatory to have minimal kinetics:
- Metabolites are in row, reactions in columns.
- The last row corresponds to total protein amount (`Protein`).
- The last column corresponds to the ribosome reaction (`Ribosome`), producing the total protein amount.
- All metabolites starting with `x_` are external metabolites with constant concentration.
- :warning: Stoichiometric coefficients must be converted following GBA formalism (see the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">units conversion tutorial</a>).

For example, the toy model B has the following mass fraction matrix:

|             | **rxn1** | **rnx2** | **Ribosome** |
|:-----------:|:--------:|:--------:|:------------:|
|   **x_G**   |    -1    |     0    |       0      |
|    **G**    |     1    |    -1    |       0      |
|    **AA**   |     0    |     1    |      -1      |
| **Protein** |     0    |     0    |       1      |

### 2.3) Forward and backward $k_\text{cat}$ vectors (<code>kcat.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="kcat"></a>

The file `kcat.csv` contains the vectors of forward (`kcat_f`) and backward (`kcat_b`) turnover rates $k_\text{cat}$ (usually, in h<sup>-1</sup>). This file is mandatory to have minimal kinetics:
- Reactions are in column.
- `kcat_f` and `kcat_b` vectors are in row.
- For forward irreversible reactions, backward values will be zero.
- :warning: $k_\text{cat}$ values must be converted following GBA formalism (see the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">units conversion tutorial</a>).

For example, the toy model B has the following $k_\text{cat}$ vectors:

|            | **rxn1** | **rnx2** | **Ribosome** |
|:----------:|:--------:|:--------:|:------------:|
| **kcat_f** |    150   |    50    |     4.55     |
| **kcat_b** |     0    |     0    |       0      |

### 2.4) Michaelis constants matrix $\mathbf{K}$ (<code>K.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="K"></a>

The file `K.csv` contains the matrix of Michaelis constants $\mathbf{K}$ (usually, in g.L<sup>-1</sup>). This file is mandatory to have minimal kinetics:
- Metabolites are in row, reactions in columns (as in the matrix $\mathbf{M}$).
- The matrix maps Michaelis constants from reactions to substrates and products, therefore including forward and backward $K_\text{M}$ values.
- :warning: $K_\text{M}$ values must be converted following GBA formalism (see the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">units conversion tutorial</a>).

For example, the toy model B has the following Michaelis constant matrix:

|             | **rxn1** | **rnx2** | **Ribosome** |
|:-----------:|:--------:|:--------:|:------------:|
|   **x_G**   |    10    |     0    |       0      |
|    **G**    |     0    |    10    |       0      |
|    **AA**   |     0    |     0    |      8.3     |
| **Protein** |     0    |     0    |       0      |

### 2.5) External conditions matrix (<code>conditions.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="conditions"></a>

The file `conditions.csv` contains the list of external conditions. This file is mandatory to have minimal kinetics. Each condition contains:
- Condition identifiers are in column (usually numbered from 1 to N),
- External metabolite concentrations are in row (usually, in g.L<sup>-1</sup>).
- The total density $\rho$ is located in the first row (usually, in g.L<sup>-1</sup>).

For example, the toy model B has 25 conditions with a glucose gradient. Here are the first 5 conditions. $\rho = 340 g/L$ represents the typical <em>E. coli</em> dry mass:

|            | **1** | **2** | **3** | **4** | **5** |
|:----------:|:-----:|:-----:|:-----:|-------|-------|
| **$\rho$** |  340  |  340  |  340  | 340   | 340   |
|   **x_G**  |  100  | 66.67 | 44.44 | 29.63 | 19.75 |

### 2.6) Initial solution $\mathbf{f}_0$ (<code>f0.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="f0"></a>

This file `f0.csv` contains a flux fraction vector $f_0$ being an initial valid solution for the SRC model. This solution must be generated by the user (see the Python module <a href="https://github.com/charlesrocabert/gbapy" target="_blank">gbapy</a>). This solution is used as a starting point for the gradient ascent algorithm.

For example, the toy model A has the following initial solution $f_0$:

| **Reaction** | **f0** |
|:------------:|:------:|
|   **rxn1**   |   1.0  |
|   **rxn2**   |  0.97  |
| **Ribosome** | 0.93   |

### 2.7) Activation constants matrix $\mathbf{K_\text{A}}$ (<code>KA.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="KA"></a>

The optional file `KA.csv` contains activation constants $K_\text{A}$ (usually in g.L<sup>-1</sup>), where some metabolites acts as activators of one or more reactions. The structure is the same than the Michaelis constants matrix.

### 2.8) Inhibition constants matrix $\mathbf{K_\text{I}}$ (<code>KI.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="KI"></a>

The optional file `KI.csv` contains inhibition constants $K_\text{I}$ (usually in g.L<sup>-1</sup>), where some metabolites acts as inhibitors of one or more reactions. The structure is the same than the Michaelis constants matrix.

### 2.9) List of constant reactions (<code>constant_reactions.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="constant_reactions"></a>

The optional file `constant_reactions.csv` contains a list of reactions and their associated values in the flux fraction vector $f$, that are kept constant during optimization. For example:

| **Reaction** | **Value** |
|:------------:|:------:|
|   **r1**   |   0.2  |
|   **r2**   |  0.01  |

### 2.10) Enzyme to protein mass concentration mapping (<code>protein_contributions.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="protein_contributions"></a>

The optional file `protein_contributions.csv` contains a mapping linking enzyme mass concentrations (the vector $p$ in GBA formalism) and protein mass concentrations. This file is useful to build predicted proteomics of a SRC model and compare it to experimental datasets.

:warning: This file could only be obtained with biological knowledge of the constructed SRC model, see the Python module <a href="https://github.com/charlesrocabert/gbapy" target="_blank">gbapy</a>.

Here is an example:

| **Reaction** | **Protein**  | **Contribution**   |
|:------------:|:------------:|:------------------:|
| DADK         | protein_0651 | 1.0                |
| DADNK        | protein_0330 | 0.979131688089874  |
| DADNK        | protein_0382 | 1.0208683119101258 |
| DADNabc      | protein_0008 | 0.1358904650088565 |
| DADNabc      | protein_0009 | 0.3862941601280572 |
