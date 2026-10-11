#!/usr/bin/env python3
"""Verify the immutable research snapshot into a new directory outside this repo."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SNAPSHOT = ROOT / "source/resource-review-2026-10-10"
MANIFEST = ROOT / "source/IMPORT_MANIFEST.json"


def inventory():
    if SNAPSHOT.is_symlink() or not SNAPSHOT.is_dir():
        raise ValueError("The source snapshot must be a regular directory.")
    files = {}
    for path in sorted(SNAPSHOT.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Unexpected source symlink: {path}")
        if path.is_file():
            data = path.read_bytes()
            files[path.relative_to(SNAPSHOT).as_posix()] = {
                "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()
            }
    return files


def differences(expected, actual):
    return {
        "missing": sorted(expected.keys() - actual.keys()),
        "extra": sorted(actual.keys() - expected.keys()),
        "changed": sorted(name for name in expected.keys() & actual.keys()
                          if expected[name] != actual[name]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir.expanduser().resolve()
    if output.exists() or output == ROOT or ROOT in output.parents:
        parser.error("Use a new output directory outside the repository.")
    output.mkdir(parents=True)
    receipt = {
        "status": "FAIL",
        "snapshot": SNAPSHOT.relative_to(ROOT).as_posix(),
        "scientific_assertions_pass": False,
        "canonical_reports_byte_identical": False,
        "source_integrity_pass": False,
        "total_scientific_groups": 0,
        "source_verifier_returncode": None,
        "source_receipt": "source-run/RECEIPT.json",
        "source_checks": {},
        "errors": [],
    }
    checks = receipt["source_checks"]
    manifest_bytes = expected = before = None
    try:
        manifest_bytes = MANIFEST.read_bytes()
        manifest = json.loads(manifest_bytes)
        if manifest["root"] != SNAPSHOT.name:
            raise ValueError("Import manifest names a different source snapshot.")
        expected = manifest["files"]
        receipt["expected_files"] = len(expected)
        receipt["import_manifest_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()
        before = inventory()
        checks["before_differences"] = differences(expected, before)
        checks["before_matches_manifest"] = before == expected
        if before != expected:
            raise ValueError("Source membership or bytes differ from the import manifest.")
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONOPTIMIZE="0",
                   OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        with (output / "source-verifier.log").open("xb") as log:
            run = subprocess.run(
                [sys.executable, "-B", str(SNAPSHOT / "verify.py"),
                 "--output-dir", str(output / "source-run")],
                cwd=SNAPSHOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                timeout=600,
            )
        receipt["source_verifier_returncode"] = run.returncode
        source = json.loads((output / receipt["source_receipt"]).read_text())
        receipt["total_scientific_groups"] = source.get("total_scientific_groups", 0)
        receipt["scientific_assertions_pass"] = (
            source.get("new_assertions_pass") is True
            and source.get("preserved_assertions_pass") is True
            and receipt["total_scientific_groups"] == 16
        )
        receipt["canonical_reports_byte_identical"] = (
            source.get("all_reports_byte_identical") is True
        )
        checks["archived_source_unchanged"] = source.get("source_unchanged") is True
        checks["archived_verifier_status"] = source.get("status")
    except (OSError, ValueError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        receipt["errors"].append(f"{type(error).__name__}: {error}")
    finally:
        try:
            after = inventory()
            checks["manifest_unchanged"] = (
                manifest_bytes is not None and MANIFEST.read_bytes() == manifest_bytes
            )
            checks["snapshot_unchanged"] = before is not None and before == after
            checks["after_matches_manifest"] = expected is not None and after == expected
            if expected is not None:
                checks["after_differences"] = differences(expected, after)
        except (OSError, ValueError, TypeError) as error:
            receipt["errors"].append(f"Post-run integrity: {type(error).__name__}: {error}")
    receipt["source_integrity_pass"] = all(checks.get(name) is True for name in (
        "before_matches_manifest", "after_matches_manifest", "manifest_unchanged",
        "snapshot_unchanged", "archived_source_unchanged",
    ))
    passed = (receipt["source_verifier_returncode"] == 0
              and checks.get("archived_verifier_status") == "PASS"
              and receipt["scientific_assertions_pass"]
              and receipt["canonical_reports_byte_identical"]
              and receipt["source_integrity_pass"] and not receipt["errors"])
    receipt["status"] = "PASS" if passed else "FAIL"
    with (output / "RECEIPT.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
