import pandas as pd
import pyaging as pya

def run_bmiq_prediction(file_path, output_name="AltumAge_Horvath_predictions.csv"):
    """
    Loads DNA methylation data, imputes missing CpGs, predicts epigenetic age using AltumAge and Horvath clocks,
    and saves predictions to CSV.
    
    Parameters:
    - file_path (str): Path to CSV file with DNA methylation data (CpGs x samples).
    - output_name (str): Output CSV filename for predictions.
    """
    print(f"Processing: {file_path}")
    
    # Load and transpose data (CpGs as rows, samples as columns)
    dnam = pd.read_csv(file_path, index_col=0)
    dnam_transposed = dnam.T
    print("DNAm shape after transpose:", dnam_transposed.shape)
    
    # Convert to AnnData with KNN imputation
    adata = pya.pp.df_to_adata(dnam_transposed, imputer_strategy='knn')
    
    # Predict epigenetic age using AltumAge and Horvath clocks
    pya.pred.predict_age(adata, ['altumAge', 'horvath2013'])
    
    # Report percent missing CpGs
    altum_na = adata.uns.get('altumage_percent_na', 'NA')
    horvath_na = adata.uns.get('horvath2013_percent_na', 'NA')
    print(f"AltumAge % missing CpGs: {altum_na}")
    print(f"Horvath % missing CpGs: {horvath_na}")
    
    # Save predictions for both clocks
    adata.obs[["altumage", "horvath2013"]].to_csv(output_name, index_label='Sample_Name')
    print(f"Saved predictions to: {output_name}\n")