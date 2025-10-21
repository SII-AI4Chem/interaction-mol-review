# interactions_view
=======
# 配体相互作用可视化流程

此文件夹包含基于 PLIP 与辅助脚本的一站式配体分析与可视化工具链。

## 目录说明

- `run_plip_ligand_map.sh`：主流程脚本（结构转换 → PLIP 分析 → 配体抽取 → 2D 相互作用图）。
- `extract_ligand.py`：从 PDB 中抽取指定配体残基。
- `mol_map6.py`：根据 PLIP 生成的 XML 和配体文件绘制 SVG 交互图。
- `plip/`：本地 PLIP 源码（脚本内部调用 `python -m plip.plipcmd`）。
- 示例结构文件：`9hvx.cif`、`b3.cif`。
- 示例运行结果位于 `example_output/`（见下文）。

## 快速上手

确保已激活 `plip-env` 环境，在当前目录下运行：

```bash
scripts/run_plip_ligand_map.sh <structure.cif|pdb> [hetid] [chain] [resid] [scale] [output_dir]
```

- 如果 `hetid` / `chain` / `resid` 传入 `auto` 或省略，脚本会读取 PLIP 的 `report.xml`，自动选取第一个存在相互作用的配体。
- `scale` 默认为 `55`，将传给 `mol_map6.py` 控制图像缩放。
- `output_dir` 默认为当前目录；脚本会在其中创建 `plip_output_<HETID>_<CHAIN>_<RESID>/` 子目录。

示例（自动识别配体并将结果输出到 `my_results/`）：

```bash
scripts/run_plip_ligand_map.sh 9hvx.cif auto auto auto 60 ./my_results
```

## 示例说明

我们已使用 `b3.cif` 自动识别配体并跑完流程，结果位于 `example_output/`：

- `b3.pdb`：由 CIF 转换得到的 PDB。
- `ligand_LIG_B_0.pdb`：抽取出的配体残基。
- `ligand_LIG_B_0_2d.sdf`：通过 Open Babel 生成的 2D 坐标。
- `plip_output_LIG_B_0/`：PLIP 报告及可视化产物：
  - `report.xml` / `report.txt`
  - `b3_protonated.pdb`
  - `B3_PROTEIN_LIG_B_0.pse`（PyMOL 会话文件）
  - `mol_LIG_B_0.svg`（`mol_map6.py` 生成的交互图）

打开 `example_output/plip_output_LIG_B_0/mol_LIG_B_0.svg` 即可预览最终图像。

## 环境依赖

- 需先激活已安装 PLIP、Open Babel 等依赖的 conda 环境，例如 `conda activate plip-env`。
- 辅助脚本仅依赖 Python 3 标准库。
