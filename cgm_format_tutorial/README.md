<p align="center">
  <img src="https://github.com/user-attachments/assets/4f06bbdf-ef2f-4500-b775-fd8cfe9dd590" width=250 />
</p>

<p align="center">
  <h3 align="center">Cell growth model (CGM) format tutorial</h3>
</p>

-----------------

Cell growth models (CGMs) must comply to a standard format to be compatible with growth balance analysis (GBA) formalism. Please refer to <a href="https://doi.org/10.1371/journal.pcbi.1011156" target="_blank">Dourado et al. (2023)</a> and the tutorial available on https://cellgrowthsim.com/ for other detailed sources.

Two formats are used to distribute CGMs: the text format CSV and the OpenDocument Spreadsheet format ODS. Data organization is strictly identical in both formats. As for now, <strong>gbacpp</strong> reads CGMs in CSV format only.

# 1) Files organization

The CGM is organized as a set of CSV files located in a folder having the name of the model. They usually contain GBA variables (such as matrices, or vectors), but also additional variables and information.

      └── CGM_name
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

Some files are mandatory to set a minimal CGM. These are:
- M.csv
- kcat.csv
- K.csv
- KA.csv
- conditions.csv
- f0.csv

Other files are optional.

⚠️ All CSV files have `;` separators, and `\n` line breaks.

# 2) Files content

### 2.1) CGM information (<code>Info.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="info"></a>

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
           └── Description: Simplest CGM with two reactions
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

The mandatory file `M.csv` contains the mass fraction matrix $\mathbf{M}$, which is the pendant of the stoichiometric matrix in normalized mass units (see <a href="https://doi.org/10.1371/journal.pcbi.1011156" target="_blank">Dourado et al., 2023</a>):
- Metabolites are in row, reactions in columns.
- The last row corresponds to total protein amount (`Protein`).
- The last column corresponds to the ribosome reaction (`Ribosome`), producing the total protein amount.
- All metabolites starting with `x_` are external metabolites with constant concentration.

For example, the toy model B has the following mass fraction matrix:

|             | **rxn1** | **rnx2** | **Ribosome** |
|:-----------:|:--------:|:--------:|:------------:|
|   **x_G**   |    -1    |     0    |       0      |
|    **G**    |     1    |    -1    |       0      |
|    **AA**   |     0    |     1    |      -1      |
| **Protein** |     0    |     0    |       1      |

### 2.3) Forward and backward $k_\text{cat}$ vectors (<code>kcat.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="kcat"></a>
### 2.4) Michaelis constants matrix $\mathbf{K}$ (<code>K.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="K"></a>
### 2.4) External conditions matrix (<code>conditions.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="conditions"></a>
### 2.4) Initial solution $\mathbf{f}_0$ (<code>f0.csv</code>) <img src="https://img.shields.io/badge/mandatory-red" /> <a name="f0"></a>
### 2.5) Activation constants matrix $\mathbf{K_A}$ (<code>KA.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="KA"></a>
### 2.6) Inhibition constants matrix $\mathbf{K_I}$ (<code>KI.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="KI"></a>
### 2.6) List of constant metabolite fractions in the initial solution (<code>constant_rhs.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="constant_rhs"></a>
### 2.6) List of constant reactions (<code>constant_reactions.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="constant_reactions"></a>
### 2.6) Enzyme to protein mass concentration mapping (<code>protein_contributions.csv</code>) <img src="https://img.shields.io/badge/optional-grey" /> <a name="protein_contributions"></a>



