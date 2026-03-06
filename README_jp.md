# auxiliary for protomidpy


### Weight correction
ALMA の Measurement Set (MS) からnpz を作成し、円盤の幾何を用いて WEIGHT のバイアスを推定・補正するためのコード。

* `ms_to_npz.py` — Measurement Set を npz に変換する。`protomidpy` の入力にも使用可能。
* `compute_bias_weight.py`— 円盤幾何（中心, cosi, PA）を用いて軸対称性を仮定し重みのバイアスを推定する。
* `correct_weight_ms.py`  — 推定したバイアスに基づき MS の WEIGHT（必要なら SIGMA も）を補正する。
* `average_ms.py` - Measurement setをChannel & Time Averageする。データ量を減らしたい時に便利。

### Producing Residual Measurement set
* `ms_to_npz_spw.py` — Measurement Setの(u,v)をspwごとの配列にして出力。residual計算に必須。
* `mode_calc_spw.py` — protomidpyの結果からモデルを計算。加えて、Meaurement setに対応するようにモデルを計算。
* `make_residual_ms.py` - protomidpyで作成したVisibilityのモデルをMeasurement Setのデータから引く。

---

## Quick Start
---

- Weight correction

```bash
# 1) MS → NPZ 変換
casa -c ms_to_npz.py　(ターミナルでCASAを開かず打つ）

npzは`protomidpy`にも使える。

# 2) WEIGHT バイアスの推定（円盤幾何が必要）
python compute_bias_weight.py

# 3) バイアスの適用（
casa -c correct_weight_ms.py　(ターミナルでCASAを開かず打つ）

## Tips
- 自分のmsファイル名、出力ファイル、円盤幾何へ変更すること
- `average_ms.py`を用いて、msをaveragingしてサイズを下げると解析がやりやすい。
```
---

- Producing Residual Measurement set
```bash
# 1) MS → NPZ 変換
casa -c ms_to_npz_for_spw.py

# 2) protomidpyのmodel calc (MSからmodelを引くために必須)
python model_calc_spw.py

# 3) MSからmodelを引く
casa -c make_residual_ms.py　(ターミナルでCASAを開かず打つ）
- modelfile:  
  protomidpy の model file へのパス（`model_calc.py` で作成される `model.npz`）

- msfile:  
  protomidpy で解析した元の Measurement Set へのパス

- cos_dec:  
  cos(declination)。declination は天体 (RA, Dec) の Dec をラジアン単位にしたもの

```
---


