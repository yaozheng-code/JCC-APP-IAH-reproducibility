# APP–IAH研究｜MIMIC-IV v3.1本地数据库复现代码

**研究题目：** *Abdominal perfusion pressure at confirmation of intra-abdominal hypertension and subsequent mortality in critically ill patients: a retrospective cohort study*

本代码包提供一套可独立执行的复现流程，面向已经获得MIMIC-IV访问资格、并在本地PostgreSQL部署MIMIC-IV v3.1及相关派生概念表的研究者。

**数据来源仅为本地MIMIC-IV v3.1数据库。** 本包不含患者数据、既往数据导出或预先构建的分析队列。复现者无需提供原研究的患者级CSV文件。

## 一、复现内容

1. 按研究定义从MIMIC-IV筛选成人ICU住院，依据两次膀胱压测量确认IAH，建立T0和分析队列。
2. 从数据库原生表及派生表提取APP、IAP、血压、实验室结果、输入输出、尿量、器官功能和病情调整变量。
3. 按研究采用的六药VIS定义及原始数据单位计算血管活性药物评分。
4. 保留原SQL规定的来源记录缺失机制，不将未记录的测量误作零。
5. 使用随包的原统计代码生成主要/敏感性分析、表格、图形和汇总基准核查。

**复现类型：独立近似复现。** 不同机构的派生表构建版本、数据库部署、计算环境等可能造成小幅数值差异。评价时应同时核查队列、主要效应方向、效应量、置信区间及统计结果。代码不会为了匹配预期指标而调整患者记录或修改模型。

## 二、数据库前提

- 完整的MIMIC-IV **v3.1** PostgreSQL数据库；
- `mimiciv_icu`、`mimiciv_hosp`、`mimiciv_derived`数据库模式；
- 已建立`icustay_detail`、`vitalsign`、`bg`、`chemistry`、`urine_output`、`vasoactive_agent`、`sofa`、`kdigo_stages`、`charlson`、`weight_durations`、`sepsis3`、`ventilation`等所需派生表；
- Python 3.11及以上、能够以只读方式查询上述PostgreSQL数据。

本项目不提供数据库安装包、MIMIC原始数据或派生表创建所需的全部外部脚本。应根据PhysioNet官方授权及MIMIC衍生概念文档完成本地准备。

## 三、直接运行

下载并解压代码包，在包含`run_all.py`的根目录中执行：

```bash
python -m pip install -r requirements.txt
python run_all.py --database mimiciv --bootstrap 500 --seed 42
```

如需先测试流程，可使用`--bootstrap 200`，但正式复现时推荐按照原统计程序的500次Bootstrap配置运行。

可指定独立输出位置：

```bash
python run_all.py --database mimiciv \
  --outdir "$HOME/mimic-project/results/jcc_app_iah_reproduction" \
  --bootstrap 500 --seed 42
```

Jupyter使用者可参考`Run_One_Click_Jupyter.ipynb`和`JUPYTER_一格启动_示例.txt`。

## 四、代码结构

| 文件/文件夹 | 作用 |
|---|---|
| `run_all.py` | 自动完成从本地数据库到论文统计输出的完整流程 |
| `stages/01_cohort_*.sql`、`02_cohort_*.sql` | 队列筛选SQL |
| `stages/native_extract_*.sql` | 十类临床数据轨迹查询 |
| `stages/native_stage1.py` | IAP确认T0及T0相关临床指标 |
| `stages/native_stage2.py` | SOFA、AKI、Charlson、体重等协变量 |
| `stages/compute_vis.py` | 六药VIS计算 |
| `stages/source_eligibility.py` | 输入、输出、尿量的来源记录与缺失处理 |
| `stages/run_app_iah.py` | 原论文统计分析实现（保持原样） |
| `METHODS_VARIABLES.md` | 变量、时间窗及研究方法 |
| `README_EN.md` | 英文使用说明 |
| `benchmarks/REFERENCE_AGGREGATES.json` | 汇总参考指标 |
| `tests/test_core_logic.py` | 不使用患者数据的算法单元测试 |

## 五、输出和使用规范

默认结果目录：

`~/mimic-project/results/jcc_app_iah_native_approx_reproduction/`

`01_cohort`、`02_native_source_trajectories`、`03_native_covariates`及`04_local_clinical_dataset`包含患者级编号或临床记录，**只能保留在授权本地环境中，不应公开上传或随论文作为开放数据附件分发**。

`05_statistical_results`中包含表格、模型结果和图形，项目根目录的`REPRODUCTION_SUMMARY_AGGREGATE_ONLY.json`包含总体运行摘要。所有对外输出均应先检查是否符合MIMIC数据使用协议以及期刊要求。

## 六、验证和报告

研究者已在本地部署的MIMIC-IV v3.1数据库上成功完成整合后的一键流程端到端运行，配套的离线逻辑测试和统计回归测试也已通过。**该结果说明代码在研究者本地环境下能够执行，但不保证所有计算机和数据库部署均可直接运行。**实际运行取决于数据库及派生概念表版本、Python依赖、PostgreSQL配置、访问权限和计算资源等条件。复现者应记录所用环境、检查运行日志，并核对实际输出结果。

自动生成的汇总基准核查用`PASS`/`REVIEW`标记结果。复现者应保留完整运行日志、数据库及派生表版本、Python环境和汇总质量报告，并如实报告实际结果。

## 七、来源与引用

请依据PhysioNet对MIMIC-IV v3.1数据库及其派生概念的引用要求标注数据来源。此代码包不授予MIMIC-IV数据访问资格。
