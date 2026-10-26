# EMSC-NLTHG

MATLAB code for **Error-robust multi-view subspace clustering with nonconvex low-rank tensor approximation and hyper-Laplacian graph embedding**, Baicheng Pan, Chuandong Li and Hangjun Che, *Engineering Applications of Artificial Intelligence*, 133 (2024), 108274.

[Paper / DOI](https://doi.org/10.1016/j.engappai.2024.108274)

This compact release runs one complete dataset. It contains the original optimization updates, the required helper functions, and the prokaryotic dataset; parameter sweeps and other experiments are omitted.

## Requirements

- MATLAB with **Statistics and Machine Learning Toolbox** (`kmeans`, `KDTreeSearcher`, `knnsearch`).
- No compilation, FLANN binaries, or full GSPbox installation is needed.

## Quick start

Clone or download this repository, open MATLAB, change to the repository folder, and run:

```matlab
demo_prokaryotic
```

The script resolves data and helper paths relative to its own location. It sets the random seed to 0, normalizes each sample in each view, runs optimization, and performs spectral clustering 10 times (20 k-means replicates each).

For a command-line run from the repository folder:

```sh
matlab -batch "demo_prokaryotic"
```

## Included dataset

`data/prokaryotic.mat` is the dataset supplied with the author's local release:

| Variable | Shape | Meaning |
| --- | --- | --- |
| `text` | 551 x 438 | Text view |
| `proteome_comp` | 551 x 3 | Proteome composition view |
| `gene_repert` | 551 x 393 | Gene repertoire view |
| `truth` | 551 x 1 | Reference labels |

The demo transposes the three views so that each column is a sample. Reference labels are used for evaluation and the number of clusters, as in the original experiment.

## Parameters and output

The original `ER_pro.m` settings are retained in `demo_prokaryotic.m`:

| Parameter | Value |
| --- | --- |
| Error penalty `lambda1` | 0.1 |
| Graph penalty `lambda2` | 2 |
| Tensor penalty `alfa`, `theta` | 100, 10 |
| Neighbors | 3 |
| Initial ADMM `mu` | 1e-4 |
| Multiplier growth | 2 |
| Maximum `mu` | 1e11 |
| Residual tolerance | 1e-7 |

The original stopping rule allows up to 202 iterations. The script reports residuals; reaching the iteration limit does not by itself establish convergence.

Results are printed and saved to `results/prokaryotic_result.mat` (generated files are ignored by Git):

- `result`: mean and standard deviation of ACC, NMI, adjusted Rand index (`AR`), F-score, precision and recall.
- `Result`: one column per clustering run; rows are ACC, NMI, AR, F-score, precision, recall.
- `C`: predicted labels for the 10 clustering runs.
- `S`: learned affinity matrix.
- `Recons`, `Match`: reconstruction and tensor matching residual histories.
- `iter`, `time_cost`: iteration count and optimization time in seconds (excluding clustering).

This is a single-dataset demonstration. Scores may vary with MATLAB/toolbox versions and numerical libraries; it does not reproduce all paper experiments.

## Validation

The complete 551-sample, 4-class dataset was run successfully using MATLAB R2024b Update 3 and Statistics and Machine Learning Toolbox, starting from MATLAB's default path. The run took 34 optimization iterations and produced ACC 0.7241, NMI 0.4254 and adjusted Rand index 0.4433 (means over 10 clustering runs). Output dimensions and finite metric values were checked.

## Files and provenance

```text
demo_prokaryotic.m       Single-dataset entry script
core/                   Logarithmic error shrinkage and Geman tensor shrinkage
data/prokaryotic.mat     Included dataset
third_party/clustering/ Normalization, graph initialization and evaluation
third_party/gspbox/      Required GSPbox 0.7.0 helper subset and license
citation.bib            Paper citation
```

The demo is extracted from the author's `ER_pro.m` in the original experiment directory. It preserves the update order (including the error update inside the per-view loop), tensor shrinkage, hypergraph construction, parameters and evaluation routines. Packaging changes add relative paths, a fixed seed, dependency checks and saved output and actual residual logging below the stopping tolerance; commented experiments and unused initialization are removed. The helper that was named `wshrinkObj_log.m` in the local publishing folder uses the same Geman weights as `wshrinkObj_Geman.m`; this release uses the latter name from the original experiment code.

Third-party copyright notices are retained. The bundled GSPbox subset is distributed under GPL-3.0-or-later; see [its license](third_party/gspbox/LICENSE). This subset includes the author's locally modified hypergraph Laplacian implementation. Spectral clustering credits Ehsan Elhamifar, and label matching credits Deng Cai. MATLAB source encodings and trailing whitespace are normalized for publication. The only behavioral change to the bundled clustering helpers removes the global `warning off` command. No additional blanket license is assigned to the remaining author-supplied code or dataset.

## Citation

If you use this code, please cite:

```bibtex
@article{pan2024error,
  title={Error-robust multi-view subspace clustering with nonconvex low-rank tensor approximation and hyper-Laplacian graph embedding},
  author={Pan, Baicheng and Li, Chuandong and Che, Hangjun},
  journal={Engineering Applications of Artificial Intelligence},
  volume={133},
  pages={108274},
  year={2024},
  publisher={Elsevier},
  doi={10.1016/j.engappai.2024.108274}
}
```

## 中文运行说明

安装 MATLAB 和 Statistics and Machine Learning Toolbox，将当前目录切换到仓库根目录，运行 `demo_prokaryotic` 即可。数据和所需辅助函数已包含，无需手动配置路径。程序使用完整的 prokaryotic 数据集，结果保存在 `results/prokaryotic_result.mat`。此版本仅提供一个数据集的核心示例。
