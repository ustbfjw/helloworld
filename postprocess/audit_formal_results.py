"""独立审计 MATLAB 检查点和 summary.csv；只写工作区指定的审计 JSON。

不重算 bootstrap，不判断方法是否更优；不修改 MATLAB、检查点或正式结果。
依赖 NumPy / SciPy。示例见 --help。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat


RELATIVE_TOLERANCE = 1e-9
ZERO_SCALE = 1e-12  # 两数均接近零时的比较尺度，绝对容差为 1e-21。
SCHEMAS = {
    1: (
        ["static", "global_offbin", "near_inside", "near_outside", "same_frequency", "local_left_half"],
        ["median", "pn", "full", "no_spatial_gate", "no_magnetic_guard", "no_refinement"],
    ),
    2: (
        ["static", "global_offbin", "same_frequency", "gaussian_patch", "tilted_smooth", "oblique_step"],
        ["median", "pn", "full", "spatial_validated"],
    ),
}
GAIN_FIELDS = (
    "gain_rms", "gain_rms_left", "gain_rms_right", "baseline_gain_rms",
    "baseline_gain_rms_left", "baseline_gain_rms_right",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def native(value: Any) -> Any:
    # simplify_cells 在二维 MATLAB struct 数组内仍可能保留 mat_struct。
    if hasattr(value, "_fieldnames"):
        return {name: native(getattr(value, name)) for name in value._fieldnames}
    if isinstance(value, dict):
        return {name: native(item) for name, item in value.items() if not name.startswith("__")}
    if isinstance(value, np.ndarray) and value.dtype == object:
        result = np.empty(value.shape, dtype=object)
        for index, item in np.ndenumerate(value):
            result[index] = native(item)
        return result
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    return value


def read_mat(path: Path) -> dict:
    return native(loadmat(path, simplify_cells=True))


def text(value: Any) -> str:
    if isinstance(value, str):
        return value
    array = np.asarray(value)
    if array.size == 0:
        return ""
    if array.dtype.kind in "US":
        return "".join(array.reshape(-1).tolist())
    raise ValueError(f"不是 MATLAB 文本字段: {type(value).__name__}")


def number(value: Any) -> float:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError("数值字段不是标量")
    return float(array.reshape(-1)[0])


def exact_equal(left: Any, right: Any) -> bool:
    """数值完全相同（NaN 必须位于同一位置）；不使用近似比较缓存。"""
    if isinstance(left, dict) or isinstance(right, dict):
        return isinstance(left, dict) and isinstance(right, dict) and left.keys() == right.keys() and all(
            exact_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, (list, tuple)) or isinstance(right, (list, tuple)):
        return isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)) and len(left) == len(right) and all(
            exact_equal(a, b) for a, b in zip(left, right)
        )
    a, b = np.asarray(left), np.asarray(right)
    if a.shape != b.shape:
        return False
    if a.dtype.kind in "biufc" and b.dtype.kind in "biufc":
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


def cache_errors(items: np.ndarray, methods: list[str]) -> tuple[list[str], int]:
    errors, reused = [], 0
    for (condition_index, method_index), item in np.ndenumerate(items):
        reference = text(item["reused_from"])
        if not reference:
            continue
        reused += 1
        if reference not in methods or methods.index(reference) >= method_index:
            errors.append(f"condition {condition_index + 1}, {methods[method_index]}: 缓存来源不是先前方法 {reference!r}")
        elif not exact_equal(item["stats"], items[condition_index, methods.index(reference)]["stats"]):
            errors.append(f"condition {condition_index + 1}, {methods[method_index]}: 缓存统计与 {reference} 不完全一致")
    return errors, reused


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def aggregate(items: list[dict], p: dict, diameter: float, schema: int) -> dict[str, float]:
    """只使用充分统计；用中心化内积实现独立 OLS，不调用 MATLAB 拟合。"""
    stats = [item["stats"] for item in items]
    nl, nf = int(p["ana"]["n_lag"]), int(p["ana"]["n_fit"])
    counts, moments = [], []
    for stat in stats:
        c = np.asarray(stat["counts"], dtype=float).reshape(-1)
        m = np.asarray(stat["moments"], dtype=float)
        require(c.shape == (nl,) and m.shape == (nl, 6), "充分统计维度错误")
        require(bool(np.all(np.isfinite(c)) and np.all(c >= 0) and np.all(np.isfinite(m))), "位移统计非有限或计数为负")
        require(bool(np.all(np.isfinite(stat["driftA"])) and np.all(np.isfinite(stat["driftB"]))), "漂移统计非有限")
        counts.append(c)
        moments.append(m)
    counts_array = np.stack(counts)
    total_count = counts_array.sum(axis=0)
    total_moment = np.stack(moments).sum(axis=0)
    y = np.full(nl, np.nan)
    np.divide(total_moment[:, 1], total_count, out=y, where=total_count > 0)
    tau = np.arange(1, nl + 1, dtype=float) / float(p["sim"]["fps"])
    valid = np.isfinite(y) & (np.arange(nl) < nf)
    diffusion = math.nan
    if valid.sum() >= 2:
        centered_t = tau[valid] - tau[valid].mean()
        centered_y = y[valid] - y[valid].mean()
        diffusion = float(np.dot(centered_t, centered_y) / np.dot(centered_t, centered_t) / 2)
    phys = p["phys"]
    se = float(phys["kB"]) * float(phys["T"]) / (3 * np.pi * float(phys["eta"]))
    diameter_estimate = se / (diffusion * 1e-12) * 1e9 if math.isfinite(diffusion) and diffusion > 0 else math.nan
    expected = {
        "diameter_true_nm": diameter,
        "D_true_um2_s": se / (diameter * 1e-9) * 1e12,
        "n_videos": len(items),
        "active_videos": int(np.any(counts_array > 0, axis=1).sum()),
        "D_um2_s": diffusion,
        "diameter_nm": diameter_estimate,
    }
    for name in ("n_tracks", "n_loc", "unmatched_selected", "matched_localizations"):
        values = [number(stat[name]) for stat in stats]
        require(all(math.isfinite(x) and x >= 0 and x.is_integer() for x in values), f"无效计数: {name}")
        expected[name] = sum(values)
    squared = sum(number(stat["localization_squared_sum_nm2"]) for stat in stats)
    expected["localization_rms_nm"] = math.sqrt(squared / expected["matched_localizations"]) if expected["matched_localizations"] > 0 else math.nan
    metrics = [item["metric"] for item in items]
    for name in GAIN_FIELDS:
        values = np.array([number(metric[name]) for metric in metrics])
        require(bool(np.all(np.isfinite(values)) and np.all(values >= 0)), f"照明 RMS 非有限或为负: {name}")
        expected[name + "_pct"] = float(np.sqrt(np.mean(values ** 2)) * 100)
    baseline = expected["baseline_gain_rms_pct"]
    expected["gain_rms_reduction_pct"] = 100 * (1 - expected["gain_rms_pct"] / baseline) if baseline > 1e-12 else math.nan
    applied = [number(metric["applied"]) for metric in metrics]
    require(all(x in (0, 1) for x in applied), "metric.applied 非逻辑标量")
    expected["applied_videos"] = sum(applied)
    expected["saturated_pixels"] = sum(number(metric["saturated_pixels"]) for metric in metrics)
    expected["adequate_evidence"] = int(expected["active_videos"] >= 10 and expected["n_loc"] >= 5000)
    if schema == 2:
        spatial = [item["illumination_info"].get("is_spatial", 0) for item in items]
        require(all(number(x) in (0, 1) for x in spatial), "is_spatial 非逻辑标量")
        expected["spatial_accepted_videos"] = sum(number(x) for x in spatial)
        expected["spatial_decision_videos"] = sum("is_spatial" in item["illumination_info"] for item in items)
    return expected


def audit(data_dir: Path, source_root: Path | None = None) -> dict:
    report = {
        "audit_version": 1, "data_directory": str(data_dir.resolve()), "status": "failed",
        "scope": "checkpoint_integrity_and_independent_point_estimates_only",
        "scientific_performance_judgement": "not_assessed",
        "bootstrap_or_interval_recomputed": False,
        "relative_tolerance": RELATIVE_TOLERANCE, "near_zero_comparison_scale": ZERO_SCALE,
        "errors": [], "missing": [], "groups": [], "rows": [],
        "source_hash_check": "not_requested", "parameters_hash_check": "not_in_protocol",
        "cache_reuse_checked": 0, "numeric_comparisons": 0, "max_scaled_difference": 0.0,
    }
    errors, missing = report["errors"], report["missing"]
    try:
        protocol_path = data_dir / "protocol.json"
        stored_hash = (data_dir / "protocol.sha256").read_text(encoding="utf-8").strip()
        require(sha256(protocol_path) == stored_hash, "protocol.json 实际 SHA256 与 protocol.sha256 不符")
        protocol = json.loads(protocol_path.read_text(encoding="utf-8-sig"))
        schema = int(protocol["schema_version"])
        require(schema in SCHEMAS, "不支持的协议版本")
        conditions, methods = SCHEMAS[schema]
        require([c["name"] for c in protocol["conditions"]] == conditions and protocol["methods"] == methods, "条件或方法顺序与协议schema不符")
        groups = protocol["groups"]
        require([g["name"] for g in groups] == ["mono30", "mono100", "mono200"], "三组名称错误")
        require([g["n_videos"] for g in groups] == [20, 60, 60], "检查点计划必须为20/60/60")
        require([g["diameter_nm"] for g in groups] == [30, 100, 200], "组名义粒径错误")
        report.update(schema_version=schema, protocol_hash=stored_hash, expected_checkpoints=140,
                      synthetic_statistics_fixture=bool(protocol.get("synthetic_statistics_fixture", False)))
        frozen_parameters = None
        if "parameters_sha256" in protocol:
            parameter_path = data_dir / "parameters.mat"
            require(sha256(parameter_path) == protocol["parameters_sha256"], "parameters.mat SHA256 不符")
            frozen_parameters = read_mat(parameter_path)["parameters"]
            require(len(frozen_parameters) == 3, "冻结参数组数不为3")
            report["parameters_hash_check"] = "matched"
        if source_root is not None:
            source_root = source_root.resolve()
            files, hashes = protocol["source_files"], protocol["source_hashes"]
            require(len(files) == len(hashes) and len(set(files)) == len(files) and len(files) > 0, "源码清单长度/唯一性错误")
            for name, expected_hash in zip(files, hashes):
                path = (source_root / name).resolve()
                require(path.is_relative_to(source_root), f"源码清单路径越界: {name}")
                require(sha256(path) == expected_hash, f"冻结源码 SHA256 不符: {name}")
            report["source_hash_check"] = "matched"
            report["source_files_checked"] = len(files)

        expected_files = {}
        for group in groups:
            wanted = {f"video_{index:03d}.mat" for index in range(1, group["n_videos"] + 1)}
            present = {path.name for path in (data_dir / group["name"]).glob("video_*.mat")}
            extras = sorted(present - wanted)
            require(not extras, f"{group['name']} 有多余检查点: {extras}")
            missing.extend(f"{group['name']}/{name}" for name in sorted(wanted - present))
            report["groups"].append({"group": group["name"], "expected": len(wanted), "present": len(present)})
            expected_files[group["name"]] = sorted(wanted)
        summary_path = data_dir / "summary.csv"
        if not summary_path.is_file():
            missing.append("summary.csv")
        report["present_checkpoints"] = sum(group["present"] for group in report["groups"])
        if missing:
            report["status"] = "incomplete"
            return report  # 未全齐时不将子集的统计结果包装成完整验收。
        with summary_path.open(newline="", encoding="utf-8-sig") as handle:
            csv_rows = list(csv.DictReader(handle))
        expected_keys = {(g["name"], c, m) for g in groups for c in conditions for m in methods}
        rows = {(row["group"], row["condition"], row["method"]): row for row in csv_rows}
        require(len(rows) == len(csv_rows) and set(rows) == expected_keys, "summary.csv 行数、唯一键或条件/方法集合错误")
        report["summary_rows"] = len(rows)

        for group_index, group in enumerate(groups):
            checkpoints, reference_p = [], None
            for video_index, filename in enumerate(expected_files[group["name"]], start=1):
                path = data_dir / group["name"] / filename
                data = read_mat(path)
                require(number(data["video_index"]) == video_index, f"{path.name}: video_index不符")
                expected_seed = group["seed_base"] + protocol["seed_stride"] * (video_index - 1)
                require(number(data["seed"]) == expected_seed, f"{group['name']}/{path.name}: seed不符")
                require(text(data["protocol_hash"]) == stored_hash, f"{path.name}: 协议hash不符")
                items = data["items"]
                require(isinstance(items, np.ndarray) and items.shape == (len(conditions), len(methods)), f"{path.name}: items矩阵维度错误")
                p = data["p"]
                require(text(p["ana"]["fit_model"]) == "perp", "主估计器并非perp")
                require(number(p["sim"]["n_frames"]) == protocol["frames"] and number(p["sim"]["n_videos"]) == group["n_videos"], "检查点帧数/视频计划与协议不符")
                if reference_p is None:
                    reference_p = p
                require(exact_equal(p, reference_p), f"{path.name}: 组内参数不完全一致")
                if frozen_parameters is not None:
                    require(exact_equal(p, frozen_parameters[group_index]), f"{path.name}: 参数与冻结parameters.mat不完全一致")
                problems, reused = cache_errors(items, methods)
                report["cache_reuse_checked"] += reused
                errors.extend(f"{group['name']}/{path.name}: {message}" for message in problems)
                checkpoints.append(items)
            for condition_index, condition in enumerate(conditions):
                for method_index, method in enumerate(methods):
                    key = (group["name"], condition, method)
                    expected = aggregate([cp[condition_index, method_index] for cp in checkpoints], reference_p,
                                         float(group["diameter_nm"]), schema)
                    report["rows"].append({"group": key[0], "condition": condition, "method": method, **expected})
                    for field, estimate in expected.items():
                        actual = float(rows[key][field])
                        report["numeric_comparisons"] += 1
                        if math.isnan(estimate) and math.isnan(actual):
                            continue
                        if not math.isfinite(estimate) or not math.isfinite(actual):
                            errors.append(f"{'/'.join(key)} {field}: NaN/非有限状态不符: independent={estimate}, csv={actual}")
                            continue
                        scaled = abs(actual - estimate) / max(abs(estimate), abs(actual), ZERO_SCALE)
                        report["max_scaled_difference"] = max(report["max_scaled_difference"], scaled)
                        if scaled > RELATIVE_TOLERANCE:
                            errors.append(f"{'/'.join(key)} {field}: independent={estimate:.17g}, csv={actual:.17g}, scaled_difference={scaled:.3g}")
        report["status"] = "failed" if errors else "integrity_and_point_estimates_match"
    except Exception as exc:
        # MAT解析损坏、无效元数据或数值异常均须明确失败并留下审计报告。
        errors.append(f"{type(exc).__name__}: {exc}")
        report["status"] = "failed"
    return report


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", type=Path, help="包含protocol及mono30/100/200检查点的目录")
    parser.add_argument("--output", type=Path, required=True, help="work内的审计JSON文件，不得位于被审计数据目录")
    parser.add_argument("--source-root", type=Path, help="可选：按协议中源码清单核对全部源码SHA256")
    args = parser.parse_args()
    work_root, output, data_dir = Path(__file__).resolve().parent, args.output.resolve(), args.data_dir.resolve()
    if not output.is_relative_to(work_root) or output.is_relative_to(data_dir) or output.suffix.lower() != ".json":
        parser.error("--output必须是脚本所在work目录内、被审计数据目录外的.json文件")
    report = audit(data_dir, args.source_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(json_safe(report), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    print(f"{report['status']}: checkpoints={report.get('present_checkpoints', 0)}/140, "
          f"rows={report.get('summary_rows', 0)}, errors={len(report['errors'])}, "
          f"missing={len(report['missing'])}; {output}")
    return 0 if report["status"] == "integrity_and_point_estimates_match" else 2 if report["status"] == "incomplete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
