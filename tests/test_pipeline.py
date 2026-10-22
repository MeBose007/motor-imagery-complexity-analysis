"""End-to-end run of every script on a small synthetic dataset (see synthetic_data.py)."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import synthetic_data

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def run(script, argv):
    spec = importlib.util.spec_from_file_location(Path(script).stem, SCRIPTS / script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main(argv)


@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("strokepe")
    cfgdir = synthetic_data.make(tmp)
    mp = pytest.MonkeyPatch()
    mp.setenv("STROKEPE_CONFIG_DIR", str(cfgdir))
    yield tmp
    mp.undo()


@pytest.mark.parametrize("mode", ["paper", "published"])
def test_main_pipeline(workspace, mode):
    for s in ["01_compute_pe.py", "02_tables.py", "03_figure2.py", "04_figure3.py", "05_figure4.py"]:
        run(s, ["--mode", mode])
    out = workspace / "results" / "broadband" / mode
    pe = pd.read_csv(out / "pe_table.csv")
    n_seg = 1 if mode == "paper" else 20
    # 8 stroke x 2 hands + 4 EO + 4 EC subjects, 25 channels each
    assert len(pe) == (8 * 2 + 4 + 4) * 25 * n_seg
    assert pe["pe"].between(0, 1).all()
    for f in ["table2.csv", "table3.csv", "figure2.png", "figure3ab.png", "figure3c-h.png", "figure4a.png", "figure4b-e.png"]:
        assert (out / f).exists(), f
    t3 = pd.read_csv(out / "table3.csv")
    assert len(t3) == 15


def test_modes_differ_but_agree_closely(workspace):
    a = pd.read_csv(workspace / "results/broadband/paper/pe_sub.csv")
    b = pd.read_csv(workspace / "results/broadband/published/pe_sub.csv")
    m = a.merge(b, on=["dataset", "subject", "condition"], suffixes=("_paper", "_pub"))
    d = (m["pe_paper"] - m["pe_pub"]).abs()
    assert d.max() > 0 and d.max() < 0.1


def test_check_script_runs(workspace):
    run("06_check_published.py", ["--mode", "published"])
    assert (workspace / "results/broadband/published/check_report.md").exists()


def test_supplementary(workspace):
    run("supplementary/tables_S1_S2_anova.py", [])
    run("supplementary/tables_S1_S2_anova.py", ["--level", "subject"])
    run("supplementary/table_S3_cbpt.py", ["--n-perm", "99", "--k", "2", "3"])
    run("supplementary/figS2_tac.py", [])
    run("supplementary/figS7_S8_psd.py", [])
    sup = workspace / "results" / "broadband" / "supplementary"
    s1 = pd.read_csv(sup / "tableS1_sensor_anova_trial-level.csv")
    assert len(s1) == 25 and s1["ANOVA_p (FDR)"].notna().all()
    s3 = pd.read_csv(sup / "tableS3_cbpt.csv")
    assert len(s3) == 12 and (s3["p"] >= 1 / 100).all()
    tac = pd.read_csv(workspace / "results" / "supplementary" / "figureS2_tac_mu3.csv")
    assert {"Minor L-MIT", "EC-RS", "Pink noise"} <= set(tac.columns)
    for f in ["figureS2.png", "figureS7_psd_eyes_closed.png", "figureS8_psd_eyes_open.png"]:
        assert (workspace / "results" / "supplementary" / f).exists()


def test_alpha_band(workspace):
    for s in ["01_compute_pe.py", "02_tables.py", "03_figure2.py"]:
        run(s, ["--mode", "published", "--band", "alpha"])
    assert (workspace / "results/alpha/published/figureS3.png").exists()
