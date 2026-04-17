"""
GCMS Web 应用 — FastAPI 后端
启动: python app.py
"""
import sys, os, re, json, time, random, uuid, logging, threading, sqlite3, io
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from bs4 import BeautifulSoup
from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.responses import HTMLResponse, StreamingResponse, JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from starlette.requests import Request

app = FastAPI(title="GCMS 分析工具")
app.mount("/static", StaticFiles(directory=str(Path(__file__).parent / "templates")), name="static")
HTML_PATH = str(Path(__file__).parent / "templates" / "index.html")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# ======================== 常量 ========================
ISTD_CAS = "123-96-6"    # 内标物CAS号，需要请自改
ISTD_CONC = 4.0    # 内标物浓度，请根据实际情况填写
DB_PATH = str(Path(__file__).parent / "chemical_data.db")
HEADERS_WEB = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Accept-Language": "en-US,en;q=0.9",
    "Referer": "https://www.chemicalbook.com/"
}
DELAY_RANGE = (1, 3)
MAX_RETRIES = 3

# ======================== 内存存储 ========================
sessions: dict = {}

# ======================== txt 解析 ========================

def sanitize_for_excel(df: pd.DataFrame) -> pd.DataFrame:
    """清除 Excel 不支持的控制字符"""
    illegal = re.compile(r'[\x00-\x08\x0B\x0C\x0E-\x1F]')
    def clean(c):
        return illegal.sub('', c) if isinstance(c, str) else c
    df_clean = df.copy()
    for col in df_clean.select_dtypes(include='object'):
        mask = df_clean[col].notna()
        df_clean.loc[mask, col] = df_clean.loc[mask, col].apply(clean)
    return df_clean


def normalize_cas(cas_str):
    if not isinstance(cas_str, str):
        return ""
    return re.sub(r"\s*-\s*", "-", cas_str.strip())


def is_placeholder_cas(cas_str) -> bool:
    cas = normalize_cas(cas_str)
    if not cas:
        return True
    digits = re.sub(r"\D", "", cas)
    return bool(digits) and set(digits) == {"0"}


def normalize_compound_name(name) -> str:
    if not isinstance(name, str):
        return ""
    return re.sub(r"\s+", " ", name.strip()).lower()


def build_compound_group_key(cas, name, peak_no=None, hit=None) -> str:
    cas = normalize_cas(cas)
    if cas and not is_placeholder_cas(cas):
        return f"cas::{cas}"
    norm_name = normalize_compound_name(name)
    if norm_name:
        return f"name::{norm_name}"
    if peak_no is not None and hit is not None:
        return f"peak::{peak_no}::hit::{hit}"
    if peak_no is not None:
        return f"peak::{peak_no}"
    return "unknown"

def read_sections(text: str) -> dict:
    sections = {}
    current = None
    for line in text.split("\n"):
        line = line.rstrip("\r")
        m = re.match(r"^\[(.+)\]$", line)
        if m:
            current = m.group(1)
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    return sections
# PLACEHOLDER_PARSE

def parse_mc_peak_table(lines):
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("Peak#\t"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("MC Peak Table 中未找到表头行")
    headers = lines[header_idx].split("\t")
    rows = []
    for line in lines[header_idx + 1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        padded = parts + [""] * max(0, len(headers) - len(parts))
        rows.append(dict(zip(headers, padded[:len(headers)])))
    df = pd.DataFrame(rows)
    for col in ["Peak#", "Area", "Height", "SI"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["Ret.Time", "A/H", "Conc.", "Area%", "Height%"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "CAS #" in df.columns:
        df["CAS #"] = df["CAS #"].apply(normalize_cas)
    return df


def parse_spectrum_process_table(lines):
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("Spectrum#\t"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Spectrum Process Table 中未找到表头行")
    headers = lines[header_idx].split("\t")
    rows = []
    for line in lines[header_idx + 1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        padded = parts + [""] * max(0, len(headers) - len(parts))
        rows.append(dict(zip(headers, padded[:len(headers)])))
    df = pd.DataFrame(rows)
    for col in ["Spectrum#", "Hit #", "SI"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "CAS #" in df.columns:
        df["CAS #"] = df["CAS #"].apply(normalize_cas)
    if "Name" in df.columns:
        df["Name"] = df["Name"].apply(lambda x: x.split("$$")[0].strip() if isinstance(x, str) else x)
    return df


def parse_chromatogram(lines):
    """解析 TIC 色谱图数据"""
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("Ret.Time\t"):
            data_start = i + 1
            break
    if data_start is None:
        return [], []
    times, intensities = [], []
    for line in lines[data_start:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            try:
                times.append(float(parts[0]))
                intensities.append(float(parts[1]))
            except ValueError:
                continue
    return times, intensities


def parse_spectra(sections):
    """解析所有峰的质谱数据，返回 {peak_no: {mz: [...], intensity: [...]}}"""
    spectra = {}
    # 找到所有 MS Spectrum section（它们在 sections 里可能重复 key）
    # 需要重新从原始文本解析
    return spectra


def parse_spectra_from_text(text: str, peak_count: int):
    """从原始文本解析所有质谱"""
    spectra = {}
    lines = text.split("\n")
    i = 0
    peak_idx = 0
    while i < len(lines):
        line = lines[i].rstrip("\r")
        if line == "[MS Spectrum]":
            peak_idx += 1
            if peak_idx > peak_count:
                break
            mz_list, int_list = [], []
            # 跳到数据行
            i += 1
            data_started = False
            while i < len(lines):
                l = lines[i].rstrip("\r")
                if l.startswith("["):
                    break
                if l.startswith("m/z\t"):
                    data_started = True
                    i += 1
                    continue
                if data_started and l.strip():
                    parts = l.split("\t")
                    if len(parts) >= 2:
                        try:
                            mz_list.append(float(parts[0]))
                            int_list.append(float(parts[1]))
                        except ValueError:
                            pass
                i += 1
            spectra[peak_idx] = {"mz": mz_list, "intensity": int_list}
            continue
        i += 1
    return spectra

# PLACEHOLDER_CHEM

# ======================== ChemicalBook 爬虫 ========================

def get_db_conn():
    db = DB_PATH
    if not os.path.exists(db):
        os.makedirs(os.path.dirname(db), exist_ok=True)
        conn = sqlite3.connect(db, check_same_thread=False)
        conn.execute("""CREATE TABLE IF NOT EXISTS chemical_data (
            cas TEXT PRIMARY KEY, data_json TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP)""")
        conn.commit()
        return conn
    return sqlite3.connect(db, check_same_thread=False)


def db_get(conn, cas):
    row = conn.execute("SELECT data_json FROM chemical_data WHERE cas = ?", (cas,)).fetchone()
    return json.loads(row[0]) if row else None


def db_put(conn, cas, data):
    conn.execute("INSERT OR REPLACE INTO chemical_data (cas, data_json, updated_at) VALUES (?, ?, CURRENT_TIMESTAMP)",
                 (cas, json.dumps(data, ensure_ascii=False)))
    conn.commit()


def parse_chemicalbook_html(html, cas):
    if not html:
        return None
    soup = BeautifulSoup(html, "html.parser")
    result = {"CAS": cas}
    basics_div = soup.find("div", class_="Basicsl")
    if basics_div:
        for item in basics_div.find_all("div"):
            if item.span and item.span.text.strip():
                key = item.span.text.strip().replace("：", "")
                value = item.get_text().replace(key, "").strip()
                result[key] = value
    subclass_div = soup.find("div", id="SubClass")
    if subclass_div:
        for section in subclass_div.find_all("div", class_="sxlist"):
            title = section.find("h2")
            if not title:
                continue
            section_name = title.text.strip()
            if section_name == "基本信息":
                for item in section.select(".cwb > div.tbt"):
                    key = item.text.strip()
                    sib = item.find_next_sibling("span")
                    value = sib.text.strip() if sib else ""
                    result[f"基本信息_{key}"] = value.replace("\n", " | ")
            elif section_name == "物理化学性质":
                for item in section.select(".xztable .xztr"):
                    if item.span:
                        key = item.span.text.strip()
                        value = item.get_text().replace(key, "").strip()
                        result[f"物化性质_{key}"] = value
    return result


def fetch_chemical_data(cas, conn, retry=0):
    cached = db_get(conn, cas)
    if cached:
        return cached
    try:
        time.sleep(random.uniform(*DELAY_RANGE))
        url = f"https://www.chemicalbook.com/CAS_{cas}.htm"
        resp = requests.get(url, headers=HEADERS_WEB, timeout=10)
        resp.raise_for_status()
        data = parse_chemicalbook_html(resp.text, cas)
        if data:
            db_put(conn, cas, data)
        return data
    except Exception as e:
        if retry < MAX_RETRIES:
            time.sleep(2)
            return fetch_chemical_data(cas, conn, retry + 1)
        logging.warning(f"CAS {cas} 爬取失败: {e}")
        return None


def extract_chinese_name(data):
    if not data:
        return ""
    for key in ["中文名称", "中文名", "中文别名"]:
        if key in data and data[key]:
            return data[key]
    return ""


def extract_fema(data):
    if not data:
        return ""
    for key in data:
        if "FEMA" in key.upper():
            val = str(data[key]).strip()
            nums = re.findall(r"\d{3,5}", val)
            if nums:
                return nums[0]
    return ""


def extract_properties(data):
    """提取用于 hover 展示的关键性质"""
    if not data:
        return {}
    props = {}
    mapping = {
        "odor": "物化性质_气味 (Odor)",
        "aroma_type": "物化性质_香型",
        "appearance": "物化性质_外观性状",
        "boiling_point": "物化性质_沸点",
        "category": "基本信息_所属类别",
        "usage": "应用领域_用途1",
        "mol_weight": "分子量",
        "mol_formula": "分子式",
    }
    for key, src in mapping.items():
        if src in data and data[src]:
            props[key] = str(data[src]).strip()
    # FEMA 单独提取
    fema = extract_fema(data)
    if fema:
        props["fema"] = fema
    return props

# ======================== 标准谱库查询 ========================

def ensure_spectra_table():
    """确保 spectra_library 表存在"""
    conn = get_db_conn()
    conn.execute("""CREATE TABLE IF NOT EXISTS spectra_library (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        cas TEXT NOT NULL,
        name TEXT,
        source TEXT,
        mz_json TEXT NOT NULL,
        intensity_json TEXT NOT NULL,
        comment TEXT,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )""")
    # 为查询加速建索引（IF NOT EXISTS 防止重复）
    conn.execute("CREATE INDEX IF NOT EXISTS idx_spectra_cas ON spectra_library(cas)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_spectra_name ON spectra_library(name)")
    conn.commit()
    conn.close()


def get_reference_spectrum(cas: str):
    """从 spectra_library 查标准谱，MoNA 优先"""
    cas = normalize_cas(cas)
    if not cas or is_placeholder_cas(cas):
        return None
    conn = get_db_conn()
    row = conn.execute(
        "SELECT mz_json, intensity_json, source, name FROM spectra_library "
        "WHERE cas = ? ORDER BY CASE source WHEN 'mona' THEN 0 ELSE 1 END LIMIT 1",
        (cas,)
    ).fetchone()
    conn.close()
    if not row:
        return None
    return {
        "mz": json.loads(row[0]),
        "intensity": json.loads(row[1]),
        "source": row[2],
        "name": row[3],
    }

# 启动时确保表结构存在
ensure_spectra_table()

# PLACEHOLDER_CRAWL_TASK

# ======================== 后台爬虫任务 ========================

def crawl_task(sid: str):
    """后台线程: 爬取所有 CAS 数据"""
    session = sessions[sid]
    cas_list = session["cas_to_fetch"]
    session["crawl_total"] = len(cas_list)
    session["crawl_done"] = 0
    session["crawl_running"] = True

    conn = get_db_conn()
    chem_data = {}
    for cas in cas_list:
        data = fetch_chemical_data(cas, conn)
        chem_data[cas] = data
        session["crawl_done"] += 1

    conn.close()
    session["chem_data"] = chem_data

    # 给候选表补充中文名和 FEMA
    spec_df = session["spec_df"]
    spec_df["中文名称"] = spec_df["CAS #"].map(lambda c: extract_chinese_name(chem_data.get(c, None)))
    spec_df["FEMA"] = spec_df["CAS #"].map(lambda c: extract_fema(chem_data.get(c, None)))
    session["spec_df"] = spec_df
    session["crawl_running"] = False
    logging.info(f"[{sid}] 爬虫完成, {len(chem_data)} 个 CAS")

# PLACEHOLDER_ROUTES

# ======================== API 路由 ========================

@app.get("/", response_class=HTMLResponse)
async def index():
    with open(HTML_PATH, "r", encoding="utf-8") as f:
        return HTMLResponse(f.read())


@app.post("/api/upload")
async def upload(file: UploadFile = File(...)):
    content = (await file.read()).decode("utf-8")
    sid = str(uuid.uuid4())[:8]

    sections = read_sections(content)
    if "MC Peak Table" not in sections:
        raise HTTPException(400, "txt 中未找到 [MC Peak Table]")

    spt_key = "MS Similarity Search Results for Spectrum Process Table"
    if spt_key not in sections:
        raise HTTPException(400, f"txt 中未找到 [{spt_key}]")

    peak_df = parse_mc_peak_table(sections["MC Peak Table"])
    spec_df = parse_spectrum_process_table(sections[spt_key])

    # 色谱图 — 取第二个 MS Chromatogram (TIC)，如果只有一个就取第一个
    chrom_sections = [k for k in sections if k == "MS Chromatogram"]
    chrom_times, chrom_intensities = [], []
    if "MS Chromatogram" in sections:
        chrom_times, chrom_intensities = parse_chromatogram(sections["MS Chromatogram"])

    # 如果上面没拿到数据（因为重复 key 只保留了最后一个），从原始文本重新解析
    if not chrom_times:
        lines = content.split("\n")
        in_chrom = False
        data_started = False
        for line in lines:
            line = line.rstrip("\r")
            if line == "[MS Chromatogram]":
                in_chrom = True
                data_started = False
                chrom_times, chrom_intensities = [], []
                continue
            if in_chrom and line.startswith("[") and line != "[MS Chromatogram]":
                if chrom_times:
                    break
                in_chrom = False
                continue
            if in_chrom and line.startswith("Ret.Time\t"):
                data_started = True
                continue
            if in_chrom and data_started and line.strip():
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        chrom_times.append(float(parts[0]))
                        chrom_intensities.append(float(parts[1]))
                    except ValueError:
                        pass

    # 质谱
    spectra = parse_spectra_from_text(content, len(peak_df))

    # 收集需要爬取的 CAS
    all_cas = spec_df["CAS #"].dropna().unique().tolist()
    cas_to_fetch = [c for c in set(all_cas) if c and not is_placeholder_cas(c)]

    # 初始化选择状态
    selections = {}
    for _, row in peak_df.iterrows():
        pn = int(row["Peak#"]) if pd.notna(row["Peak#"]) else 0
        selections[pn] = None  # None = 未选择

    sessions[sid] = {
        "peak_df": peak_df,
        "spec_df": spec_df,
        "chrom_times": chrom_times,
        "chrom_intensities": chrom_intensities,
        "spectra": spectra,
        "selections": selections,
        "cas_to_fetch": cas_to_fetch,
        "chem_data": {},
        "crawl_total": len(cas_to_fetch),
        "crawl_done": 0,
        "crawl_running": True,
    }

    # 启动后台爬虫
    t = threading.Thread(target=crawl_task, args=(sid,), daemon=True)
    t.start()

    return {"sid": sid, "peaks": len(peak_df), "candidates": len(spec_df), "cas_count": len(cas_to_fetch)}


@app.get("/api/status/{sid}")
async def status(sid: str):
    s = sessions.get(sid)
    if not s:
        raise HTTPException(404, "Session not found")
    return {"total": s["crawl_total"], "done": s["crawl_done"], "running": s["crawl_running"]}


@app.get("/api/chromatogram/{sid}")
async def chromatogram(sid: str):
    s = sessions.get(sid)
    if not s:
        raise HTTPException(404)
    times = s["chrom_times"]
    intensities = s["chrom_intensities"]
    # 降采样到 ~3000 点
    step = max(1, len(times) // 3000)
    return {"times": times[::step], "intensities": intensities[::step]}


@app.get("/api/peaks/{sid}")
async def peaks(sid: str):
    s = sessions.get(sid)
    if not s:
        raise HTTPException(404)
    df = s["peak_df"]
    cols = ["Peak#", "Ret.Time", "Area", "Area%", "Name", "SI", "CAS #", "Mark"]
    cols = [c for c in cols if c in df.columns]
    result = []
    for _, row in df.iterrows():
        pn = int(row["Peak#"]) if pd.notna(row["Peak#"]) else 0
        sel = s["selections"].get(pn)
        r = {c: (None if pd.isna(row[c]) else row[c]) for c in cols}
        r["selected"] = sel
        result.append(r)
    return result

# PLACEHOLDER_ROUTES2

@app.get("/api/peak/{sid}/{peak_no}")
async def peak_detail(sid: str, peak_no: int):
    s = sessions.get(sid)
    if not s:
        raise HTTPException(404)
    df = s["peak_df"]
    row = df[df["Peak#"] == peak_no]
    if row.empty:
        raise HTTPException(404, f"Peak {peak_no} not found")
    row = row.iloc[0]

    ret_time = float(row["Ret.Time"]) if pd.notna(row["Ret.Time"]) else 0
    proc_from = float(row.get("Proc.From", ret_time - 0.5)) if pd.notna(row.get("Proc.From")) else ret_time - 0.5
    proc_to = float(row.get("Proc.To", ret_time + 0.5)) if pd.notna(row.get("Proc.To")) else ret_time + 0.5

    # 局部色谱图: 峰前后 ±1min
    margin = 1.0
    t_start, t_end = proc_from - margin, proc_to + margin
    local_t, local_i = [], []
    for t, intensity in zip(s["chrom_times"], s["chrom_intensities"]):
        if t_start <= t <= t_end:
            local_t.append(t)
            local_i.append(intensity)

    # 质谱
    spectrum = s["spectra"].get(peak_no, {"mz": [], "intensity": []})

    # 候选化合物
    spec_df = s["spec_df"]
    chem_data = s.get("chem_data", {})
    candidates = spec_df[spec_df["Spectrum#"] == peak_no]
    cand_list = []
    for _, c in candidates.iterrows():
        cas = c.get("CAS #", "")
        props = extract_properties(chem_data.get(cas))
        cand_list.append({
            "hit": int(c["Hit #"]) if pd.notna(c["Hit #"]) else 0,
            "si": int(c["SI"]) if pd.notna(c["SI"]) else 0,
            "cas": cas,
            "name": c.get("Name", ""),
            "chinese_name": c.get("中文名称", "") if pd.notna(c.get("中文名称", "")) else "",
            "fema": c.get("FEMA", "") if pd.notna(c.get("FEMA", "")) else "",
            "mol_weight": c.get("Mol.Weight", ""),
            "mol_form": c.get("Mol.Form", ""),
            "properties": props,
        })

    return {
        "peak_no": peak_no,
        "ret_time": ret_time,
        "area": float(row["Area"]) if pd.notna(row["Area"]) else 0,
        "area_pct": float(row["Area%"]) if pd.notna(row["Area%"]) else 0,
        "name": row.get("Name", ""),
        "si": int(row["SI"]) if pd.notna(row["SI"]) else 0,
        "cas": row.get("CAS #", ""),
        "local_chrom": {"times": local_t, "intensities": local_i},
        "spectrum": spectrum,
        "candidates": cand_list,
        "selected": s["selections"].get(peak_no),
        "proc_from": proc_from,
        "proc_to": proc_to,
    }


@app.get("/api/ref_spectrum/{cas}")
async def ref_spectrum(cas: str):
    spec = get_reference_spectrum(cas)
    if not spec:
        return {"found": False}
    return {"found": True, **spec}


@app.post("/api/select/{sid}/{peak_no}")
async def select_compound(sid: str, peak_no: int, request: Request):
    s = sessions.get(sid)
    if not s:
        raise HTTPException(404)
    body = await request.json()
    # body: {"hit": 1, "cas": "xxx", "name": "xxx", "chinese_name": "xxx"} 或 null
    if not isinstance(body, dict) or body.get("cas") is None:
        s["selections"][peak_no] = None
    else:
        s["selections"][peak_no] = {
            "hit": body.get("hit"),
            "cas": body.get("cas", ""),
            "name": body.get("name", ""),
            "chinese_name": body.get("chinese_name", ""),
        }
    return {"ok": True}


@app.post("/api/quantify/{sid}")
async def quantify(sid: str):
    s = sessions.get(sid)
    if not s:
        raise HTTPException(404)

    peak_df = s["peak_df"]
    selections = s["selections"]
    chem_data = s["chem_data"]

    # 找内标物 Area
    istd_row = peak_df[peak_df["CAS #"] == ISTD_CAS]
    if istd_row.empty:
        raise HTTPException(400, f"未找到内标物 (CAS {ISTD_CAS})")
    istd_area = float(istd_row["Area"].sum())
    if istd_area == 0:
        raise HTTPException(400, "内标物峰面积为 0")

    # 收集已选择的峰
    rows = []
    for pn, sel in selections.items():
        if sel is None:
            continue
        pr = peak_df[peak_df["Peak#"] == pn]
        if pr.empty:
            continue
        pr = pr.iloc[0]
        rows.append({
            "peak_no": pn,
            "hit": sel.get("hit"),
            "cas": sel["cas"],
            "name": sel["name"],
            "chinese_name": sel.get("chinese_name", ""),
            "area": float(pr["Area"]) if pd.notna(pr["Area"]) else 0,
            "area_pct": float(pr["Area%"]) if pd.notna(pr["Area%"]) else 0,
        })

    if not rows:
        return {"results": [], "istd_area": istd_area}

    rdf = pd.DataFrame(rows)
    rdf["group_key"] = rdf.apply(
        lambda r: build_compound_group_key(r["cas"], r["name"], r["peak_no"], r["hit"]),
        axis=1,
    )

    # 真实 CAS 继续按 CAS 合并，占位 CAS 用名称/峰号兜底，避免 00-00-0 被错误合并
    grouped = rdf.groupby("group_key", sort=False).agg({
        "cas": "first",
        "area": "sum", "area_pct": "sum",
        "name": "first", "chinese_name": "first",
        "peak_no": lambda x: ",".join(str(v) for v in dict.fromkeys(x)),
    }).reset_index()

    grouped["concentration"] = (grouped["area"] / istd_area) * ISTD_CONC
    grouped["fema"] = grouped["cas"].map(
        lambda c: extract_fema(chem_data.get(c)) if c and not is_placeholder_cas(c) else ""
    )

    # 排除内标物
    grouped = grouped[grouped["cas"] != ISTD_CAS]
    grouped = grouped.sort_values("concentration", ascending=False).reset_index(drop=True)
    grouped.index += 1

    results = []
    for _, r in grouped.iterrows():
        results.append({
            "cas": r["cas"],
            "chinese_name": r["chinese_name"],
            "name": r["name"],
            "fema": r["fema"],
            "area_pct": round(r["area_pct"], 2),
            "concentration": round(r["concentration"], 4),
            "peaks": r["peak_no"],
        })

    s["quant_results"] = results
    return {"results": results, "istd_area": istd_area, "istd_conc": ISTD_CONC}


@app.get("/api/export/{sid}")
async def export(sid: str):
    s = sessions.get(sid)
    if not s or "quant_results" not in s:
        raise HTTPException(400, "请先执行定量计算")

    results = s["quant_results"]
    df = pd.DataFrame(results)
    df.index += 1
    df.index.name = "序号"
    df.columns = ["CAS #", "中文名称", "英文名称", "FEMA", "Area%", "浓度(mg/ml)", "对应峰号"]

    df = sanitize_for_excel(df)
    buf = io.BytesIO()
    df.to_excel(buf, engine="openpyxl")
    buf.seek(0)
    return StreamingResponse(
        buf,
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": "attachment; filename=quantitative_results.xlsx"}
    )


# ======================== 启动 ========================

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8765)
