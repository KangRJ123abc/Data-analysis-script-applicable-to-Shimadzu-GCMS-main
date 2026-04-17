"""
用 PubChem API 为 MoNA / MassBank MSP 文件中缺少 CAS 的记录补全 CAS 号，
并将补全后的谱图写入 spectra_library 表。

断点续传：已查询的 InChIKey 结果缓存在 inchikey_cas_cache 表中，
重新运行时跳过已查过的 InChIKey。

用法:
    python enrich_cas.py --mona MoNA-export-GC-MS_Spectra.msp
    python enrich_cas.py --massbank MassBank_NISTformat.msp
    python enrich_cas.py --mona MoNA-export-GC-MS_Spectra.msp --massbank MassBank_NISTformat.msp
"""

import argparse
import json
import os
import re
import sqlite3
import sys
import time
from pathlib import Path

import requests

DB_PATH = str(Path(__file__).parent / "chemical_data.db")

# PubChem synonyms API，返回该化合物的所有同义名（含 CAS）
PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{key}/synonyms/JSON"

# 保守频率：每次请求后等待 0.4s，约 2.5 req/s，远低于 PubChem 5 req/s 限制
REQUEST_DELAY = 0.4
MAX_RETRIES = 3
RETRY_DELAY = 10  # 遇到 429/5xx 时等待秒数


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


def is_valid_cas(s: str) -> bool:
    """简单校验 CAS 格式: digits-digits-digit"""
    return bool(re.match(r"^\d{2,7}-\d{2}-\d$", s.strip()))


def get_db_conn():
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
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
    conn.execute("CREATE INDEX IF NOT EXISTS idx_spectra_cas ON spectra_library(cas)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_spectra_name ON spectra_library(name)")
    # 断点续传缓存表
    conn.execute("""CREATE TABLE IF NOT EXISTS inchikey_cas_cache (
        inchikey TEXT PRIMARY KEY,
        cas TEXT,           -- NULL 表示查过但没找到
        queried_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )""")
    conn.commit()
    return conn


def normalize_intensity(intensity_list, target_max=999):
    if not intensity_list:
        return []
    max_val = max(intensity_list)
    if max_val == 0:
        return [0.0] * len(intensity_list)
    factor = target_max / max_val
    return [round(v * factor, 1) for v in intensity_list]


def lookup_cas_from_pubchem(inchikey: str) -> str | None:
    """
    查询 PubChem synonyms，从中提取 CAS 号。
    返回 CAS 字符串，或 None（未找到 / 请求失败）。
    """
    url = PUBCHEM_URL.format(key=inchikey)
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code == 404:
                return None  # 化合物不在 PubChem 中
            if resp.status_code == 429 or resp.status_code >= 500:
                wait = RETRY_DELAY * (attempt + 1)
                print(f"  [{resp.status_code}] 等待 {wait}s 后重试 ({inchikey})")
                time.sleep(wait)
                continue
            resp.raise_for_status()
            data = resp.json()
            synonyms = data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
            for syn in synonyms:
                if is_valid_cas(syn):
                    return normalize_cas(syn)
            return None
        except requests.exceptions.Timeout:
            print(f"  超时，重试 {attempt+1}/{MAX_RETRIES} ({inchikey})")
            time.sleep(RETRY_DELAY)
        except Exception as e:
            print(f"  请求异常: {e}，重试 {attempt+1}/{MAX_RETRIES}")
            time.sleep(RETRY_DELAY)
    return None


def parse_msp_no_cas(filepath: str, source: str):
    """
    解析 MSP 文件，只 yield 没有 CAS 但有 InChIKey 的记录。
    yield dict: {inchikey, name, source, mz, intensity, comment}
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    record = {}
    in_peaks = False
    mz_list = []
    int_list = []

    def flush():
        nonlocal record, in_peaks, mz_list, int_list
        inchikey = record.get("inchikey", "")
        cas = record.get("cas", "")
        has_peaks = bool(mz_list)
        # 只处理：有 InChIKey、没有有效 CAS、有峰数据
        if inchikey and not cas and has_peaks:
            yield {
                "inchikey": inchikey,
                "name": record.get("name", ""),
                "source": source,
                "mz": list(mz_list),
                "intensity": list(int_list),
                "comment": record.get("comment", ""),
            }
        record = {}
        in_peaks = False
        mz_list.clear()
        int_list.clear()

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            yield from flush()
            continue
        if in_peaks:
            for mz_s, int_s in re.findall(r"([\d.]+)[\s:]+([\d.]+)", line):
                try:
                    mz_list.append(float(mz_s))
                    int_list.append(float(int_s))
                except ValueError:
                    pass
            continue
        upper = line.upper()
        if upper.startswith("NAME:"):
            record["name"] = line.split(":", 1)[1].strip()
        elif upper.startswith("INCHIKEY:"):
            record["inchikey"] = line.split(":", 1)[1].strip()
        elif re.match(r"^CAS[#:NO\s]", line, re.IGNORECASE):
            val = re.split(r"[:#]\s*", line, maxsplit=1)[-1].strip()
            cas = normalize_cas(val)
            if cas and not is_placeholder_cas(cas):
                record["cas"] = cas
        elif upper.startswith("COMMENT"):
            comment_text = line.split(":", 1)[1].strip() if ":" in line else ""
            record["comment"] = record.get("comment", "") + (" " if record.get("comment") else "") + comment_text
            # MoNA 风格: "cas=xxx"
            if not record.get("cas"):
                m = re.search(r'"cas=([^"]+)"', comment_text, re.IGNORECASE)
                if m:
                    cas = normalize_cas(m.group(1))
                    if cas and not is_placeholder_cas(cas):
                        record["cas"] = cas
        elif re.match(r"^NUM\s*PEAKS?\s*:", line, re.IGNORECASE):
            in_peaks = True
        elif re.match(r"^\d", line) and record.get("name"):
            pairs = re.findall(r"([\d.]+)[\s:]+([\d.]+)", line)
            if pairs:
                in_peaks = True
                for mz_s, int_s in pairs:
                    try:
                        mz_list.append(float(mz_s))
                        int_list.append(float(int_s))
                    except ValueError:
                        pass

    yield from flush()


def process_file(filepath: str, source: str, conn: sqlite3.Connection):
    print(f"\n处理文件: {filepath} (source={source})")
    records = list(parse_msp_no_cas(filepath, source))
    print(f"  无 CAS 记录数: {len(records)}")

    # 按 InChIKey 去重（同一 InChIKey 可能有多条谱图）
    by_inchikey: dict[str, list] = {}
    for rec in records:
        by_inchikey.setdefault(rec["inchikey"], []).append(rec)

    unique_keys = list(by_inchikey.keys())
    print(f"  唯一 InChIKey 数: {len(unique_keys)}")

    # 读取已缓存的查询结果
    cached = {}
    for row in conn.execute("SELECT inchikey, cas FROM inchikey_cas_cache"):
        cached[row[0]] = row[1]  # None 表示查过但没找到

    todo = [k for k in unique_keys if k not in cached]
    print(f"  已缓存: {len(cached)} 个，待查询: {len(todo)} 个")
    print(f"  预计耗时: ~{len(todo) * REQUEST_DELAY / 3600:.1f} 小时\n")

    found = 0
    inserted = 0

    for i, inchikey in enumerate(todo, 1):
        cas = lookup_cas_from_pubchem(inchikey)
        # 写入缓存（无论是否找到）
        conn.execute(
            "INSERT OR REPLACE INTO inchikey_cas_cache (inchikey, cas) VALUES (?, ?)",
            (inchikey, cas)
        )
        if (i % 100) == 0:
            conn.commit()
            print(f"  进度: {i}/{len(todo)}，已找到 CAS: {found}，已插入谱图: {inserted}")

        if not cas or is_placeholder_cas(cas):
            time.sleep(REQUEST_DELAY)
            continue

        found += 1
        cached[inchikey] = cas

        # 插入该 InChIKey 对应的所有谱图
        batch = []
        for rec in by_inchikey[inchikey]:
            norm_int = normalize_intensity(rec["intensity"])
            batch.append((cas, rec["name"], rec["source"],
                          json.dumps(rec["mz"]), json.dumps(norm_int), rec["comment"]))
        conn.executemany(
            "INSERT INTO spectra_library (cas, name, source, mz_json, intensity_json, comment) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            batch
        )
        inserted += len(batch)
        time.sleep(REQUEST_DELAY)

    conn.commit()

    # 处理之前已缓存但还没插入的（断点续传后重跑时）
    already_in_db = {row[0] for row in conn.execute(
        "SELECT DISTINCT cas FROM spectra_library WHERE source = ?", (source,)
    )}
    late_batch = []
    for inchikey, cas in cached.items():
        if not cas or is_placeholder_cas(cas):
            continue
        if cas in already_in_db:
            continue
        if inchikey not in by_inchikey:
            continue
        for rec in by_inchikey[inchikey]:
            norm_int = normalize_intensity(rec["intensity"])
            late_batch.append((cas, rec["name"], rec["source"],
                               json.dumps(rec["mz"]), json.dumps(norm_int), rec["comment"]))
    if late_batch:
        conn.executemany(
            "INSERT INTO spectra_library (cas, name, source, mz_json, intensity_json, comment) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            late_batch
        )
        conn.commit()
        inserted += len(late_batch)
        print(f"  补充插入（断点续传）: {len(late_batch)} 条")

    print(f"\n  完成: 查询 {len(todo)} 个 InChIKey，找到 CAS {found} 个，插入谱图 {inserted} 条")
    return found, inserted


def main():
    parser = argparse.ArgumentParser(description="用 PubChem API 为 MSP 文件补全 CAS 并导入谱图库")
    parser.add_argument("--mona", help="MoNA MSP 文件路径")
    parser.add_argument("--massbank", help="MassBank MSP 文件路径")
    parser.add_argument("--db", default=DB_PATH, help="SQLite 数据库路径")
    args = parser.parse_args()

    if not args.mona and not args.massbank:
        parser.print_help()
        sys.exit(1)

    conn = get_db_conn()

    total_found = total_inserted = 0

    if args.mona:
        if not os.path.exists(args.mona):
            print(f"文件不存在: {args.mona}")
            sys.exit(1)
        f, i = process_file(args.mona, "mona", conn)
        total_found += f; total_inserted += i

    if args.massbank:
        if not os.path.exists(args.massbank):
            print(f"文件不存在: {args.massbank}")
            sys.exit(1)
        f, i = process_file(args.massbank, "massbank", conn)
        total_found += f; total_inserted += i

    conn.close()

    print(f"\n{'='*50}")
    print(f"全部完成: 找到 CAS {total_found} 个，新增谱图 {total_inserted} 条")

    # 最终统计
    conn = sqlite3.connect(DB_PATH)
    total = conn.execute("SELECT COUNT(*) FROM spectra_library").fetchone()[0]
    ucas = conn.execute("SELECT COUNT(DISTINCT cas) FROM spectra_library").fetchone()[0]
    by_src = conn.execute("SELECT source, COUNT(*) FROM spectra_library GROUP BY source").fetchall()
    conn.close()
    print(f"数据库总谱图: {total}，唯一CAS: {ucas}")
    for src, cnt in by_src:
        print(f"  {src}: {cnt} 条")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
