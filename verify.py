#!/usr/bin/env python3
"""Run all scientific checks and integrity checks into a new external directory.

Scientific assertions, canonical-report byte identity, and input integrity are
reported separately. All three must pass for an overall passing receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source"
MANIFEST = SOURCE / "MANIFEST.json"
REFERENCES = ROOT / "verification/reference"
SOURCE_FILES = {
    "readout.py", "check_resources.py", "finite_repulsion.py",
    "check_finite_repulsion.py", "exact_short_time.py", "check_hard_core.py",
}
REFERENCE_HASHES = {
    "readout-resources.json": "156a8a349331c02d6c67a89805576ba1ef0ae001f8bb86994983830ac55af0d7",
    "finite-repulsion.json": "09c63d61e989e9a0567cf6b3da481f8a095c65306f6d1986917ea755ce8edf4a",
    "short-time.json": "c7b21ca5f259538dc175ff37799b4b1ec664475fd9c6ce3284494ba8f1dedf0e",
    "hard-core.json": "23ad86e618f28b70b048eed0e0ea5a675fbe9fe947c250228bafaacbdaf63416",
}
# The order preserves certificate -> parameter robustness -> finite calibration.
SUITES = (
    ("readout_resources", "source/check_resources.py", "readout-resources.json", 5, "test_groups"),
    ("finite_repulsion", "source/check_finite_repulsion.py", "finite-repulsion.json", 5, "test_groups"),
    ("short_time", "source/exact_short_time.py", "short-time.json", 1, "test_groups"),
    ("hard_core", "source/check_hard_core.py", "hard-core.json", 5, "test_groups"),
    ("complete_parity_diagnostics", "research/check_complete_parity.py", "complete-parity-diagnostics.json", 4, "scientific_groups"),
    ("complete_parity_certificate", "research/parity_certificate.py", "complete-parity-certificate.json", 1, None),
    ("channel_calibration", "research/check_channel_calibration.py", "channel-calibration.json", 3, "scientific_groups"),
    ("calibration_cost", "research/check_calibration_cost.py", "calibration-cost.json", 3, "scientific_groups"),
    ("passive_parity", "research/check_passive_parity.py", "passive-parity.json", 5, "scientific_groups"),
    ("parameter_robustness", "research/check_parameter_robustness.py", "parameter-robustness.json", 5, "scientific_groups"),
    ("finite_calibration", "research/check_finite_calibration.py", "finite-calibration.json", 2, "scientific_groups"),
)


def digest(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected a regular input file: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory():
    """Inventory every source file except the manifest and the reader index."""
    if SOURCE.is_symlink() or not SOURCE.is_dir():
        raise ValueError("The source directory must be a regular directory.")
    files = {}
    for path in sorted(SOURCE.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Unexpected source symlink: {path}")
        name = path.relative_to(SOURCE).as_posix()
        if path.is_file() and name not in ("MANIFEST.json", "README.md"):
            data = path.read_bytes()
            files[name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    digest(SOURCE / "README.md")
    return files


def reference_hashes():
    if REFERENCES.is_symlink() or not REFERENCES.is_dir():
        raise ValueError("The reference directory must be a regular directory.")
    hashes = {}
    for path in sorted(REFERENCES.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Unexpected reference symlink: {path}")
        if path.is_file():
            hashes[path.relative_to(REFERENCES).as_posix()] = digest(path)
    return hashes


def code_hashes():
    files = [ROOT / "verify.py", *sorted(SOURCE.glob("*.py")),
             *sorted((ROOT / "research").glob("*.py"))]
    return {path.relative_to(ROOT).as_posix(): digest(path) for path in files}


def differences(expected, actual):
    return {
        "missing": sorted(expected.keys() - actual.keys()),
        "extra": sorted(actual.keys() - expected.keys()),
        "changed": sorted(name for name in expected.keys() & actual.keys()
                          if expected[name] != actual[name]),
    }


def report_differences(expected, actual, path=""):
    rows = []
    if isinstance(expected, dict) and isinstance(actual, dict):
        if expected.keys() != actual.keys():
            rows.append({"path": path, "kind": "keys"})
        for key in sorted(expected.keys() & actual.keys()):
            rows.extend(report_differences(expected[key], actual[key], path + "/" + key))
    elif isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            rows.append({"path": path, "kind": "length"})
        for index, (left, right) in enumerate(zip(expected, actual)):
            rows.extend(report_differences(left, right, path + "/" + str(index)))
    elif expected != actual or type(expected) is not type(actual):
        row = {"path": path, "reference": expected, "actual": actual}
        if type(expected) in (int, float) and type(actual) in (int, float):
            row["absolute_difference"] = abs(expected - actual)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir.expanduser().resolve()
    if not __debug__ or sys.flags.optimize:
        parser.error("Assertions must be enabled; do not use Python optimization.")
    if output.exists() or output == ROOT or ROOT in output.parents:
        parser.error("Use a new output directory outside the repository.")
    output.mkdir(parents=True)
    receipt = {
        "status": "FAIL", "scientific_assertions_pass": False,
        "reference_reports_byte_identical": False, "source_integrity_pass": False,
        "code_inputs_unchanged": False, "reference_inputs_unchanged": False,
        "assertions_enabled": True, "total_scientific_groups": 0,
        "expected_scientific_groups": 39, "source_manifest": "source/MANIFEST.json",
        "reference_directory": "verification/reference", "suites": {},
        "source_checks": {}, "errors": [],
    }
    checks = receipt["source_checks"]
    manifest_bytes = expected = before = code_before = references_before = index_before = None
    try:
        digest(MANIFEST)
        manifest_bytes = MANIFEST.read_bytes()
        manifest = json.loads(manifest_bytes)
        if manifest["root"] != "source" or set(manifest["files"]) != SOURCE_FILES:
            raise ValueError("Source manifest must name exactly the six scientific scripts.")
        expected = manifest["files"]
        receipt["source_manifest_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()
        receipt["expected_source_files"] = len(expected)
        before = inventory()
        index_before = digest(SOURCE / "README.md")
        checks["before_differences"] = differences(expected, before)
        checks["before_matches_manifest"] = before == expected
        if before != expected:
            raise ValueError("Source membership or bytes differ from the source manifest.")
        references_before = reference_hashes()
        receipt["reference_sha256_before"] = references_before
        if references_before != REFERENCE_HASHES:
            raise ValueError("Canonical reference membership or bytes differ from their fixed hashes.")
        code_before = code_hashes()
        receipt["code_inputs_sha256_before"] = code_before
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONOPTIMIZE="0",
                   OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                   MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
                   NUMEXPR_NUM_THREADS="1")
        for suite, script, report_name, expected_groups, group_field in SUITES:
            row = {
                "script": script, "report": report_name,
                "log": str(Path(report_name).with_suffix(".log")),
                "expected_groups": expected_groups, "scientific_groups": 0,
                "returncode": None, "scientific_assertions_pass": False,
            }
            receipt["suites"][suite] = row
            if report_name in REFERENCE_HASHES:
                row.update(reference=f"verification/reference/{report_name}",
                           reference_report_byte_identical=False)
            try:
                report_path = output / report_name
                command = [sys.executable, "-B", str(ROOT / script),
                           "--output", str(report_path)]
                if suite in ("parameter_robustness", "finite_calibration"):
                    command.extend(["--certificate", str(output / "complete-parity-certificate.json")])
                if suite == "finite_calibration":
                    command.extend(["--robustness-certificate", str(output / "parameter-robustness.json")])
                with (output / row["log"]).open("xb") as log:
                    run = subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, timeout=600)
                row["returncode"] = run.returncode
                report_bytes = report_path.read_bytes()
                report = json.loads(report_bytes)
                if group_field is None:
                    report_pass = report.get("certificate") == (
                        "PASS: exact integer and rational enclosures; no floating arithmetic")
                    groups = int(report_pass)
                else:
                    groups = report.get(group_field, 0)
                    report_pass = report.get("status") == "PASS" and groups == expected_groups
                row["scientific_groups"] = groups
                row["scientific_assertions_pass"] = run.returncode == 0 and report_pass
                row["report_sha256"] = hashlib.sha256(report_bytes).hexdigest()
                if report_name in REFERENCE_HASHES:
                    reference_bytes = (REFERENCES / report_name).read_bytes()
                    row["reference_report_byte_identical"] = report_bytes == reference_bytes
                    row["report_differences"] = report_differences(json.loads(reference_bytes), report)
            except (OSError, ValueError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
                receipt["errors"].append(f"{suite}: {type(error).__name__}: {error}")
    except (OSError, ValueError, KeyError, TypeError) as error:
        receipt["errors"].append(f"{type(error).__name__}: {error}")
    finally:
        try:
            after = inventory()
            checks["manifest_unchanged"] = (
                manifest_bytes is not None and MANIFEST.read_bytes() == manifest_bytes)
            checks["source_unchanged"] = before is not None and before == after
            checks["index_unchanged"] = index_before is not None and digest(SOURCE / "README.md") == index_before
            checks["after_matches_manifest"] = expected is not None and after == expected
            if expected is not None:
                checks["after_differences"] = differences(expected, after)
        except (OSError, ValueError, TypeError) as error:
            receipt["errors"].append(f"Source integrity: {type(error).__name__}: {error}")
        try:
            code_after = code_hashes()
            receipt["code_inputs_sha256_after"] = code_after
            receipt["code_inputs_unchanged"] = code_before is not None and code_before == code_after
            references_after = reference_hashes()
            receipt["reference_sha256_after"] = references_after
            receipt["reference_inputs_unchanged"] = (
                references_before == references_after == REFERENCE_HASHES)
        except (OSError, ValueError) as error:
            receipt["errors"].append(f"Input integrity: {type(error).__name__}: {error}")
    rows = list(receipt["suites"].values())
    receipt["total_scientific_groups"] = sum(row["scientific_groups"] for row in rows)
    receipt["scientific_assertions_pass"] = (
        len(rows) == len(SUITES) and all(row["scientific_assertions_pass"] for row in rows)
        and receipt["total_scientific_groups"] == 39)
    canonical = [row for row in rows if "reference_report_byte_identical" in row]
    receipt["reference_reports_byte_identical"] = (
        len(canonical) == 4 and all(row["reference_report_byte_identical"] for row in canonical))
    receipt["source_integrity_pass"] = all(checks.get(name) is True for name in (
        "before_matches_manifest", "after_matches_manifest", "manifest_unchanged",
        "source_unchanged", "index_unchanged"))
    passed = all(receipt[name] for name in (
        "scientific_assertions_pass", "reference_reports_byte_identical",
        "source_integrity_pass", "code_inputs_unchanged", "reference_inputs_unchanged"
    )) and not receipt["errors"]
    receipt["status"] = "PASS" if passed else "FAIL"
    with (output / "RECEIPT.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
