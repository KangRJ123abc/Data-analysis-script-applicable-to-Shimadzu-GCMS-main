"""
MSP 质谱库导入工具 — 将 MoNA / MassBank 的 MSP 文件导入 SQLite

用法:
    python import_msp.py <msp文件路径> --source mona
    python import_msp.py <msp文件路径> --source massbank

示例:
    python import_msp.py MoNA-export-GC-MS_Spectra.msp --source mona
    python import_msp.py MassBank_NIST.msp --source massbank
"""

import argparse
import json
import os
import re
import sqlite3
import sys
from pathlib import Path

DB_PATH = str(Path(__file__).parent / "chemical_data.db")


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


def ensure_table(conn):
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
    conn.commit()


def extract_cas_from_comments(comments_str):
    """从 MoNA 风格的 Comments 行中提取 CAS 号
    例: '"cas=51794-16-2" "pubchem cid=521349"' -> '51794-16-2'
    """
    m = re.search(r'"cas=([^"]+)"', comments_str, re.IGNORECASE)
    if m:
        return normalize_cas(m.group(1))
    m = re.search(r'\bcas[=:]\s*([\d]+-[\d]+-[\d]+)', comments_str, re.IGNORECASE)
    if m:
        return normalize_cas(m.group(1))
    return ""


def parse_msp_file(filepath: str):
    """
    逐条解析 MSP 文件，yield dict:
    {name, cas, mz: [], intensity: [], comment}

    支持 MoNA 和 NIST/MassBank 两种格式：
    - MoNA: CAS 在 Comments 行的 "cas=xxx" 中，峰数据用空格分隔
    - NIST: CAS 在 CAS#: 行中，峰数据用空格/tab/分号分隔
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    record = {}
    in_peaks = False
    mz_list = []
    int_list = []

    def make_record():
        """将当前累积的数据打包为一条记录，如果有效则返回 dict，否则返回 None"""
        if not mz_list or not record.get("cas"):
            return None
        return {
            "name": record.get("name", ""),
            "cas": record["cas"],
            "mz": list(mz_list),
            "intensity": list(int_list),
            "comment": record.get("comment", ""),
        }

    def reset():
        nonlocal record, in_peaks, mz_list, int_list
        record = {}
        in_peaks = False
        mz_list = []
        int_list = []

    for raw_line in lines:
        line = raw_line.strip()

        # 空行 = 记录分隔
        if not line:
            result = make_record()
            if result:
                yield result
            reset()
            continue

        # 如果在读峰数据
        if in_peaks:
            pairs = re.findall(r"([\d.]+)[\s:]+([\d.]+)", line)
            for mz_s, int_s in pairs:
                try:
                    mz_list.append(float(mz_s))
                    int_list.append(float(int_s))
                except ValueError:
                    pass
            continue

        # 元数据行
        upper = line.upper()
        if upper.startswith("NAME:"):
            record["name"] = line.split(":", 1)[1].strip()
        elif re.match(r"^CAS[#:NO\s]", line, re.IGNORECASE):
            # NIST 格式: CAS#: 100-52-7 / CASNO: 100-52-7 / CAS: 100-52-7
            val = re.split(r"[:#]\s*", line, maxsplit=1)[-1].strip()
            cas = normalize_cas(val)
            if cas:
                record["cas"] = cas
        elif upper.startswith("COMMENT"):
            # Comments: 或 Comment:
            comment_text = line.split(":", 1)[1].strip() if ":" in line else ""
            record["comment"] = record.get("comment", "") + (" " if record.get("comment") else "") + comment_text
            # MoNA 把 CAS 藏在 Comments 里: "cas=51794-16-2"
            if not record.get("cas"):
                cas = extract_cas_from_comments(comment_text)
                if cas:
                    record["cas"] = cas
        elif re.match(r"^NUM\s*PEAKS?\s*:", line, re.IGNORECASE):
            in_peaks = True
        elif re.match(r"^\d", line) and record.get("name"):
            # 有些 MSP 没有 Num Peaks 行，直接跟峰数据
            pairs = re.findall(r"([\d.]+)[\s:]+([\d.]+)", line)
            if pairs:
                in_peaks = True
                for mz_s, int_s in pairs:
                    try:
                        mz_list.append(float(mz_s))
                        int_list.append(float(int_s))
                    except ValueError:
                        pass
        # 其他元数据行 (DB#, InChIKey, Formula, MW, etc.) 跳过

    # 文件末尾最后一条
    result = make_record()
    if result:
        yield result


def normalize_intensity(intensity_list, target_max=999):
    """将 intensity 归一化到 [0, target_max]"""
    if not intensity_list:
        return []
    max_val = max(intensity_list)
    if max_val == 0:
        return [0.0] * len(intensity_list)
    factor = target_max / max_val
    return [round(v * factor, 1) for v in intensity_list]


def main():
    parser = argparse.ArgumentParser(description="将 MSP 质谱文件导入 SQLite 数据库")
    parser.add_argument("msp_file", help="MSP 文件路径")
    parser.add_argument("--source", required=True, choices=["mona", "massbank"],
                        help="数据来源标识")
    parser.add_argument("--db", default=DB_PATH, help="SQLite 数据库路径")
    args = parser.parse_args()

    if not os.path.exists(args.msp_file):
        print(f"错误: 文件不存在 — {args.msp_file}")
        sys.exit(1)

    file_size = os.path.getsize(args.msp_file)
    print(f"文件: {args.msp_file} ({file_size / 1024 / 1024:.1f} MB)")
    print(f"来源: {args.source}")
    print(f"数据库: {args.db}")
    print()

    conn = sqlite3.connect(args.db, check_same_thread=False)
    ensure_table(conn)

    imported = 0
    skipped_no_cas = 0
    skipped_placeholder = 0
    skipped_no_peaks = 0
    total = 0

    try:
        has_tqdm = True
        from tqdm import tqdm
    except ImportError:
        has_tqdm = False
        print("提示: 安装 tqdm 可以显示进度条 (pip install tqdm)")

    records = parse_msp_file(args.msp_file)

    batch = []
    batch_size = 500

    for rec in records:
        total += 1

        cas = rec.get("cas", "")
        if not cas:
            skipped_no_cas += 1
            continue
        if is_placeholder_cas(cas):
            skipped_placeholder += 1
            continue
        if not rec.get("mz"):
            skipped_no_peaks += 1
            continue

        norm_int = normalize_intensity(rec["intensity"])

        batch.append((
            cas,
            rec.get("name", ""),
            args.source,
            json.dumps(rec["mz"]),
            json.dumps(norm_int),
            rec.get("comment", ""),
        ))

        if len(batch) >= batch_size:
            conn.executemany(
                "INSERT INTO spectra_library (cas, name, source, mz_json, intensity_json, comment) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                batch
            )
            conn.commit()
            imported += len(batch)
            batch = []

            if not has_tqdm:
                print(f"\r已导入 {imported} 条 (扫描 {total} 条) ...", end="", flush=True)

    # 最后一批
    if batch:
        conn.executemany(
            "INSERT INTO spectra_library (cas, name, source, mz_json, intensity_json, comment) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            batch
        )
        conn.commit()
        imported += len(batch)

    conn.close()

    print(f"\n\n{'='*50}")
    print(f"导入完成!")
    print(f"  扫描记录:          {total}")
    print(f"  成功导入:          {imported}")
    print(f"  跳过(无CAS):       {skipped_no_cas}")
    print(f"  跳过(占位CAS):     {skipped_placeholder}")
    print(f"  跳过(无峰数据):    {skipped_no_peaks}")
    print(f"{'='*50}")

    # 统计信息
    conn = sqlite3.connect(args.db)
    total_rows = conn.execute("SELECT COUNT(*) FROM spectra_library").fetchone()[0]
    unique_cas = conn.execute("SELECT COUNT(DISTINCT cas) FROM spectra_library").fetchone()[0]
    by_source = conn.execute("SELECT source, COUNT(*) FROM spectra_library GROUP BY source").fetchall()
    conn.close()

    print(f"\n数据库统计:")
    print(f"  总谱图数:    {total_rows}")
    print(f"  唯一CAS数:   {unique_cas}")
    for src, cnt in by_source:
        print(f"  {src}: {cnt} 条")


if __name__ == "__main__":
    main()
