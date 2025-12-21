import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd

# train data
train_df = pd.read_csv("/Users/yahavbar/Documents/METL-Rosetta/METL-Rosetta/output/train_subset.tsv", sep="\t")
test_df = pd.read_csv("/Users/yahavbar/Documents/METL-Rosetta/METL-Rosetta/output/test_subset.tsv", sep="\t")
val_df = pd.read_csv("/Users/yahavbar/Documents/METL-Rosetta/METL-Rosetta/output/val_subset.tsv", sep="\t")
print(train_df.columns)

target_col = "total_score"

models_features = {
    "model_energy_only": [
        "fa_atr", "fa_rep", "fa_sol", "fa_elec", "fa_dun"
    ],
    "model_hbond_env": [
        "hbond_bb_sc", "hbond_lr_bb", "hbond_sc", "env", "vdw"
    ],
    "model_structural": [
        "buried_all", "buried_np", "contact_all",
        "total_sasa", "rg", "pack"
    ],
    "model_all_selected": [
        "fa_atr", "fa_rep", "fa_sol", "fa_elec", "fa_dun",
        "hbond_bb_sc", "hbond_lr_bb", "hbond_sc",
        "buried_all", "buried_np", "contact_all",
        "total_sasa", "rg", "pack"
    ],
    "model_energy_and_structural": [
        "fa_atr", "fa_rep", "fa_sol", "fa_elec", "fa_dun",
        "buried_all", "buried_np", "contact_all", "total_sasa",
        "rg", "pack"
    ]
}
results = {}
for name, feats in models_features.items():
    X_train = train_df[feats]
    y_train = train_df[target_col]

    X_val = val_df[feats]
    y_val = val_df[target_col]

    model = LinearRegression()
    model.fit(X_train, y_train)

    y_val_pred = model.predict(X_val)
    mse = mean_squared_error(y_val, y_val_pred)
    r2 = r2_score(y_val, y_val_pred)

    results[name] = {
        "features": feats,
        "mse_val": mse,
        "r2_val": r2,
        "coef": model.coef_,
        "intercept": model.intercept_,
    }

for name, res in results.items():
    print(name, "MSE:", res["mse_val"], "R2:", res["r2_val"])
# MSE: mean squared error, am liebsten 0
# R^2: bestimmungsmaß, [0,1], 0 = model erklär nichts, 1 = model erklärt gut