# Data analysis script applicable to Shimadzu GCMS

基于 FastAPI + Vue3 + ECharts 的 GC-MS 文字报告交互式分析工具。

> 说明：仓库默认**不包含**本地数据库、质谱库文件和 GC-MS 原始文本报告。`chemical_data.db`、`.msp`、`.txt` 等文件会被 `.gitignore` 排除，需由使用者在本地运行时自行生成或导入。

## 功能

- 上传岛津 GCMSsolution 导出的 `.txt` 文字报告
- 自动解析 TIC 色谱图、115 个峰的质谱数据、2875 条候选化合物命中
- 自动爬取 ChemicalBook 获取化合物中文名、FEMA、气味、香型、沸点等性质（SQLite 缓存）
- 交互式 TIC 色谱图，点击峰区域查看局部色谱图 + 质谱图 + 候选化合物列表
- 点击候选行展开化合物性质面板（气味、香型、类别、外观、沸点、用途），辅助判断正确命中物
- 标准谱对比：支持导入 MoNA / MassBank 开源质谱库（MSP 格式），点击候选时显示实验谱与标准谱的镜像对比图
- 逐峰选择或丢弃化合物，内标法（2-辛醇，4 mg/ml）自动定量计算
- 结果可视化（饼图 + 柱状图）+ 导出 Excel
<img width="1920" height="869" alt="8545538f092d579a6c0e1e7c6169cd57" src="https://github.com/user-attachments/assets/5ada40ac-2229-4c0a-b322-ac1bfbba97c0" />
<img width="1909" height="866" alt="d50a88c73373fc1a6b0774fabeedce68" src="https://github.com/user-attachments/assets/7acc596b-d4be-4de4-b448-719462668c46" />




## 快速开始

### 环境准备

```bash
conda create -n gcms_web python=3.12 -y
conda activate gcms_web
pip install fastapi uvicorn pandas openpyxl requests beautifulsoup4 tqdm jinja2 python-multipart
```

### 启动

```bash
python app.py
```

浏览器打开 http://127.0.0.1:8765

首次运行后，程序会在本地创建 `chemical_data.db`，用于缓存 ChemicalBook 查询结果。

### 导入标准质谱库（可选，推荐）

从以下数据源下载 MSP 格式的质谱文件：

- **MoNA**: <https://mona.fiehnlab.ucdavis.edu/downloads> → GC-MS Spectra → MSP
- **MassBank**: <https://github.com/MassBank/MassBank-data/releases> → `MassBank_NIST.msp`

```bash
# 导入 MoNA（推荐优先）
python import_msp.py MoNA-export-GC-MS_Spectra.msp --source mona

# 导入 MassBank（补充）
python import_msp.py MassBank_NIST.msp --source massbank
```

导入后，分析界面中点击候选化合物行时将自动显示标准谱与实验谱的镜像对比图。

> 提示：导入 MSP 后，标准谱数据也会写入本地 `chemical_data.db`。该数据库文件默认不会提交到 Git 仓库。
> 注意：这两个开源的谱数据库大多以化合物的InChiKey作为索引而不是CAS
> 推荐获取到两个谱库之后，按InChiKey为索引，利用PubChem官方的API查询对应的CAS

### 使用流程

1. 拖拽上传 GC-MS 文字报告 (.txt)
2. 等待化合物数据查询完成（有缓存，重复 CAS 秒查）
3. 在 TIC 色谱图上点击峰，查看候选化合物，点击行展开性质详情
4. 为每个峰选择正确的化合物或丢弃
5. 点击「完成定量计算」，查看结果并导出 Excel

## 项目结构

```
Data-analysis-script-applicable-to-Shimadzu-GCMS/
├── app.py              # FastAPI 后端（解析、爬虫、API、定量计算）
├── enrich_cas.py       # 补充/更新化合物 CAS 信息的辅助脚本
├── import_msp.py       # 导入 MoNA / MassBank 的 MSP 标准谱
├── templates/
│   └── index.html      # Vue3 + ECharts 单页前端
└── README.md
```


## 技术栈

- 后端: Python 3.12 + FastAPI + uvicorn
- 前端: Vue3 (CDN) + ECharts 5 (CDN)，无需 npm 构建
- 数据源: ChemicalBook 爬虫 + SQLite 缓存
- 定量方法: 内标法，内标物 2-辛醇 (CAS 123-96-6)，浓度 4 mg/ml，响应因子 = 1

## 支持的报告格式

岛津 GCMSsolution 导出的文字报告 (.txt)，需包含以下段落：

- `[MC Peak Table]` — 峰表
- `[MS Similarity Search Results for Spectrum Process Table]` — 候选命中
- `[MS Chromatogram]` — TIC 色谱图数据
- `[MS Spectrum]` — 各峰质谱数据

## License

MIT
