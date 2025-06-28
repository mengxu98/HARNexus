import miceforest as mf
import pandas as pd
import numpy as np
import datetime as dt
import argparse
from sklearn.impute import SimpleImputer, KNNImputer
import sys
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Impute missing values in UKB Olink data"
    )
    parser.add_argument(
        "--input_file",
        type=str,
        default="data/ukb/olink_instance_0.csv",
        help="Path to input Olink protein data CSV file"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/ukb/",
        help="Directory to save output files"
    )
    parser.add_argument(
        "--gene_lists_dir",
        type=str,
        default="data/ukb/gene_lists/",
        help="Directory to save gene lists (kept and dropped genes)"
    )
    parser.add_argument(
        "--method",
        type=str,
        default="mean",
        choices=["mean", "median", "knn", "mice"],
        help="Imputation method: mean, median, knn, or mice"
    )
    parser.add_argument(
        "--n_neighbors",
        type=int,
        default=5,
        help="Number of neighbors for KNN imputation"
    )
    parser.add_argument(
        "--mice_iterations",
        type=int,
        default=5,
        help="Number of iterations for MICE"
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=2025,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--output_suffix",
        type=str,
        default="",
        help="Suffix to add to output filename"
    )
    parser.add_argument(
        "--drop_threshold",
        type=float,
        default=None,
        help="Drop genes with missing values exceeding this threshold (0.0-1.0). Default is None (keep all genes)."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    random_seed = args.random_seed
    
    print(f"Using imputation method: {args.method}")
    
    # Create output directories if they don't exist
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.gene_lists_dir, exist_ok=True)
    
    # Load protein data
    print(f"\nLoading data from: {args.input_file}")
    data = pd.read_csv(args.input_file)
    print("Original data shape: ", data.shape)
    print("\nOriginal column names (first 20):")
    print(list(data.columns)[:20])
    
    # Convert column names to uppercase except for 'eid'
    data.columns = ["eid"] + [col.upper() for col in data.columns if col != "eid"]
    print("\nColumn names after uppercase conversion (first 20):")
    print(list(data.columns)[:20])
    
    # Save original eid column before any processing
    original_eid = data["eid"].copy()
    
    # Calculate missing values percentage for each column
    missing_percentages = data.isna().mean()
    print("\nMissing values before imputation:")
    missing_counts = data.isna().sum().sort_values(ascending=False)
    print(missing_counts.head())
    print(f"Total missing values: {missing_counts.sum()}")
    
    # Drop genes with high missing values if threshold provided
    if args.drop_threshold is not None:
        threshold = float(args.drop_threshold)
        print(f"\nDropping genes with >{threshold*100:.1f}% missing values...")
        
        # Identify columns to drop and keep, always ensure 'eid' is kept
        cols_to_drop = [col for col in data.columns if col != "eid" and missing_percentages[col] > threshold]
        cols_to_keep = ["eid"] + [col for col in data.columns if col != "eid" and col not in cols_to_drop]
        
        # Save the list of dropped and kept genes
        timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save dropped genes
        pd.DataFrame({'gene': cols_to_drop}).to_csv(
            f"{args.gene_lists_dir}/dropped_genes_{timestamp}.csv", index=False
        )
        print(f"Dropped {len(cols_to_drop)} genes with >{threshold*100:.1f}% missing values")
        print(f"List saved to {args.gene_lists_dir}/dropped_genes_{timestamp}.csv")
        
        # Save kept genes
        kept_genes = [col for col in cols_to_keep if col != "eid"]
        pd.DataFrame({'gene': kept_genes}).to_csv(
            f"{args.gene_lists_dir}/kept_genes_{timestamp}.csv", index=False
        )
        print(f"Kept {len(kept_genes)} genes with <={threshold*100:.1f}% missing values")
        print(f"List saved to {args.gene_lists_dir}/kept_genes_{timestamp}.csv")
        
        # Actually drop the columns
        data = data[cols_to_keep]
        print(f"Data shape after dropping high-missing genes: {data.shape}")
    
    # Create a copy of data for imputation, always exclude eid
    X = data.drop(columns=["eid"]) if "eid" in data.columns else data.copy()
    
    # Perform imputation based on selected method
    suffix = args.output_suffix if args.output_suffix else args.method
    if args.drop_threshold is not None:
        suffix = f"{suffix}_drop{int(args.drop_threshold*100)}"
    output_file = f"{args.output_dir}/olink_data_imputed_{suffix}.csv"
    
    if args.method == "mice":
        # Original MICE method using miceforest
        try:
            print("\nStarting MICE imputation (miceforest)...")
            # Create a dictionary of variables to impute
            exclude = ["eid"]
            dont_impute = ["eid"]
            column_dict = {
                col: [
                    other_col
                    for other_col in data.columns
                    if other_col != col and other_col not in exclude
                ]
                for col in data.columns
                if col not in dont_impute
            }
            
            # run miceforest imputation on multiple cores
            kds = mf.ImputationKernel(
                data,
                num_datasets=1,
                variable_schema=column_dict,
                random_state=random_seed,
            )
            
            # run
            kds.mice(iterations=args.mice_iterations, n_jobs=-1, verbose=True)
            
            # get the completed dataframe from the miceforest object
            olink_data_imputed = kds.complete_data()
            
        except Exception as e:
            print(f"MICE imputation failed: {e}")
            print("Falling back to median imputation...")
            imputer = SimpleImputer(strategy="median")
            X_imputed = imputer.fit_transform(X)
            olink_data_imputed = pd.DataFrame(X_imputed, columns=X.columns)
            
            # Add original eid column back
            if "eid" not in olink_data_imputed.columns:
                olink_data_imputed.insert(0, "eid", original_eid)
                
            suffix = f"median_fallback_{suffix}"
            output_file = f"{args.output_dir}/olink_data_imputed_{suffix}.csv"
            
    elif args.method == "knn":
        print(f"\nStarting KNN imputation with {args.n_neighbors} neighbors...")
        imputer = KNNImputer(n_neighbors=args.n_neighbors, weights="uniform")
        X_imputed = imputer.fit_transform(X)
        olink_data_imputed = pd.DataFrame(X_imputed, columns=X.columns)
        
        # Add original eid column back
        if "eid" not in olink_data_imputed.columns:
            olink_data_imputed.insert(0, "eid", original_eid)
        
    else:  # mean or median
        print(f"\nStarting {args.method} imputation...")
        imputer = SimpleImputer(strategy=args.method)
        X_imputed = imputer.fit_transform(X)
        olink_data_imputed = pd.DataFrame(X_imputed, columns=X.columns)
        
        # Add original eid column back
        if "eid" not in olink_data_imputed.columns:
            olink_data_imputed.insert(0, "eid", original_eid)
    
    print("\nMissing values after imputation:")
    missing_after = olink_data_imputed.isna().sum().sum()
    print(missing_after)
    
    if missing_after > 0:
        print(
            "Warning: Some values could not be imputed. Filling remaining NAs with column medians..."
        )
        olink_data_imputed = olink_data_imputed.fillna(olink_data_imputed.median())
        print(
            f"Missing values after final cleanup: {olink_data_imputed.isna().sum().sum()}"
        )
    
    # Final check to ensure eid is in the output
    if "eid" not in olink_data_imputed.columns:
        olink_data_imputed.insert(0, "eid", original_eid)
        print("Added eid column back to the final output")
    
    # Save imputed data
    olink_data_imputed.to_csv(output_file, index=False)
    print(f"\nImputed data saved to: {output_file}")
    print(f"Final data shape: {olink_data_imputed.shape}")
    print(f"Final data columns include 'eid': {'eid' in olink_data_imputed.columns}")


if __name__ == "__main__":
    main()
