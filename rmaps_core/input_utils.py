from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

def _excel_col_to_idx(col_ref: str) -> int:
    value = 0
    for ch in col_ref:
        value = value * 26 + (ord(ch) - ord("A") + 1)
    return value - 1


def _read_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []

    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for si in root.findall("a:si", _NS):
        text = "".join((node.text or "") for node in si.findall(".//a:t", _NS))
        strings.append(text)
    return strings


def _xlsx_sheet1_to_tsv(xlsx_path: Path, tsv_path: Path) -> None:
    with zipfile.ZipFile(xlsx_path) as zf:
        if "xl/worksheets/sheet1.xml" not in zf.namelist():
            raise RuntimeError(f"Could not find sheet1 in workbook: {xlsx_path}")

        shared = _read_shared_strings(zf)
        sheet_root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))

        rows: list[dict[int, str]] = []
        max_col = -1
        for row in sheet_root.findall("a:sheetData/a:row", _NS):
            cell_map: dict[int, str] = {}
            for cell in row.findall("a:c", _NS):
                ref = cell.attrib.get("r", "")
                match = re.match(r"([A-Z]+)(\d+)$", ref)
                if not match:
                    continue

                col_idx = _excel_col_to_idx(match.group(1))
                max_col = max(max_col, col_idx)

                cell_type = cell.attrib.get("t")
                value = ""
                if cell_type == "inlineStr":
                    node = cell.find("a:is/a:t", _NS)
                    value = node.text if node is not None and node.text is not None else ""
                else:
                    node = cell.find("a:v", _NS)
                    if node is not None and node.text is not None:
                        raw = node.text
                        if cell_type == "s":
                            try:
                                value = shared[int(raw)]
                            except (ValueError, IndexError):
                                value = raw
                        else:
                            value = raw

                cell_map[col_idx] = value
            rows.append(cell_map)

    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with tsv_path.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        width = max_col + 1
        for row in rows:
            writer.writerow([row.get(idx, "") for idx in range(width)])


def maybe_prepare_rmats_input(rmats: str, output_dir: Path) -> str:
    """
    Return an rMATS path that legacy scripts can parse.

    Legacy scripts read plain tab-delimited text. If an `.xlsx` workbook
    is provided, convert `sheet1` to UTF-8 TSV in `<output>/temp/` and return it.
    """
    if rmats == "NA":
        return rmats

    rmats_path = Path(rmats)
    if rmats_path.suffix.lower() != ".xlsx":
        return rmats

    converted_path = output_dir / "temp" / f"{rmats_path.stem}.from_xlsx.tsv"
    _xlsx_sheet1_to_tsv(rmats_path, converted_path)
    return str(converted_path)
