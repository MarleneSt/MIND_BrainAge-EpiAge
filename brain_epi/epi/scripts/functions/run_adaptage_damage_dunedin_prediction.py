import pandas as pd
import pyaging as pya
import pyreadr

def run_adaptage_damage_dunedin_prediction(file_path, output_name="pyaging_predictions.csv"):
    """
    Loads R data file with DNA methylation data (unnamed object), imputes missing CpGs,
    predicts epigenetic age using YingDamAge and YingAdaptAge clocks,
    and saves predictions to CSV.
    
    Parameters:
    - file_path (str): Path to R data file (.rds, .RData, etc).
    - output_name (str): Output CSV filename for predictions.
    """
    print(f"Processing: {file_path}")

    # Read R object (assumes unnamed object)
    result = pyreadr.read_r(file_path)
    
    if None not in result:
        raise ValueError("No unnamed object found in the R data file.")
    
    dnam = result[None]

    dnam_transposed = dnam.T
    print("DNAm shape after transpose:", dnam_transposed.shape)
    
    # Convert to AnnData with KNN imputation
    adata = pya.pp.df_to_adata(dnam_transposed, imputer_strategy='knn')
    
    # Predict epigenetic age
    pya.pred.predict_age(adata, ['YingDamAge', 'YingAdaptAge', 'dunedinpace'])
    
    # Print percent missing CpGs
    dam_na = adata.uns.get('yingdamage_percent_na', 'NA')
    adapt_na = adata.uns.get('yingadaptage_percent_na', 'NA')
    dunedin_na = adata.uns.get('dunedinpace_percent_na', 'NA')
    print("DamAge % missing CpGs:", dam_na)
    print("AdaptAge % missing CpGs:", adapt_na)
    print("DunedinPace % missing CpGs:", dunedin_na)
    
    # Save predictions (check columns exist)
    cols = ["yingdamage", "yingadaptage", "dunedinpace"]
    if not all(col in adata.obs.columns for col in cols):
        raise ValueError(f"Expected columns {cols} not found in adata.obs")
    
    adata.obs[cols].to_csv(output_name, index_label='Sample_Name')
    print(f"Saved predictions to: {output_name}\n")