#!/bin/bash

# Step 1: Compute granular age
echo "Step 1: Computing granular age..."
if [ ! -f "data/ukb/granular_age_april_07_2025.csv" ]; then
    echo "  Computing granular age..."
    Rscript code/aging_model/step01-UKB-granular-age.R
else
    echo "Granular age already computed"
fi

# Step 2: Impute missing values in Olink data
echo "Step 2: Imputing missing values in Olink data..."
if [ ! -f "data/ukb/olink_data_imputed_mean_drop30.csv" ]; then
    echo "  Imputing missing values in Olink data..."
    python code/aging_model/step02-UKB-olink-imputation.py \
        --input_file "data/ukb/olink_instance_0.csv" \
        --output_dir "data/ukb/" \
        --gene_lists_dir "data/ukb/gene_lists/" \
        --method mean \
        --drop_threshold 0.3 \
        --output_suffix "mean_drop30"
else
    echo "Olink data imputation already complete"
fi


# if TRUE, re-train models, and existing models will be removed
re_train_models=FALSE

features_num_boruta_list=(100)
features_num_optimized_list=(60)

for features_num_boruta in ${features_num_boruta_list[@]}; do
    for features_num_optimized in ${features_num_optimized_list[@]}; do

        echo "Step 3: Feature selection..."

        # Simplified process:
        # 1. First use LightGBM for preliminary feature selection
        # 2. Then apply Boruta on the selected features
        # 3. Finally use the final selected features to train the model

        n_boruta=100
        nrounds=2000
        train_ratio=0.7
        seed=2025

        if [ $re_train_models = TRUE ]; then
            echo "  Removing existing features selection results..."
            rm -rf results/aging_model_${features_num_optimized}_${features_num_boruta}
        fi

        echo "Parameters set: "
        echo "  features_num_boruta=$features_num_boruta"
        echo "  n_boruta=$n_boruta"
        echo "  features_num_optimized=$features_num_optimized"
        echo "  train_ratio=$train_ratio"
        echo "  seed=$seed"
        echo "Starting feature selection..."

        output_features_dir="results/aging_model_${features_num_optimized}_${features_num_boruta}/features_selection"

        if [ ! -f "$output_features_dir/har/features_optimized.csv" ]; then
            echo "  Selecting features for HAR model..."
            Rscript code/aging_model/step03-feature_selection.R \
                --feature_file "results/networks_pfc/s8_s9_common_targets.csv" \
                --data_file "data/ukb/olink_data_imputed_mean_drop30.csv" \
                --response_file "data/ukb/granular_age_april_07_2025.csv" \
                --output_dir "$output_features_dir/har/" \
                --features_num_boruta $features_num_boruta \
                --n_boruta $n_boruta \
                --nrounds $nrounds \
                --features_num_optimized $features_num_optimized \
                --train_ratio $train_ratio \
                --seed $seed
        else
            echo "  Aging model feature selection already complete for HAR"
        fi

        if [ ! -f "$output_features_dir/hagr/features_optimized.csv" ]; then
            echo "Selecting features for HAGR model..."
            Rscript code/aging_model/step03-feature_selection.R \
                --feature_file "data/genes_set/HAGR_108genes.csv" \
                --data_file "data/ukb/olink_data_imputed_mean_drop30.csv" \
                --response_file "data/ukb/granular_age_april_07_2025.csv" \
                --output_dir "$output_features_dir/hagr/" \
                --features_num_boruta $features_num_boruta \
                --n_boruta $n_boruta \
                --nrounds $nrounds \
                --features_num_optimized $features_num_optimized \
                --train_ratio $train_ratio \
                --seed $seed
        else
            echo "Aging model feature selection already complete for HAGR"
        fi

        if [ ! -f "$output_features_dir/protage/features_optimized.csv" ]; then
            echo "  Selecting features for ProtAge model..."
            Rscript code/aging_model/step03-feature_selection.R \
                --feature_file "data/genes_set/ProtAge_204genes.csv" \
                --data_file "data/ukb/olink_data_imputed_mean_drop30.csv" \
                --response_file "data/ukb/granular_age_april_07_2025.csv" \
                --output_dir "$output_features_dir/protage/" \
                --features_num_boruta $features_num_boruta \
                --n_boruta $n_boruta \
                --nrounds $nrounds \
                --features_num_optimized $features_num_optimized \
                --train_ratio $train_ratio \
                --seed $seed
        else
            echo "  Aging model feature selection already complete for ProtAge"
        fi

        if [ ! -f "$output_features_dir/syngo/features_optimized.csv" ]; then
            echo "  Selecting features for SynGO model..."
            Rscript code/aging_model/step03-feature_selection.R \
                --feature_file "data/genes_set/SynGO_1026genes.csv" \
                --data_file "data/ukb/olink_data_imputed_mean_drop30.csv" \
                --response_file "data/ukb/granular_age_april_07_2025.csv" \
                --output_dir "$output_features_dir/syngo/" \
                --features_num_boruta $features_num_boruta \
                --n_boruta $n_boruta \
                --nrounds $nrounds \
                --features_num_optimized $features_num_optimized \
                --train_ratio $train_ratio \
                --seed $seed
        else
            echo "  Aging model feature selection already complete for SynGO"
        fi

        if [ ! -f "$output_features_dir/brain/features_optimized.csv" ]; then
            echo "  Selecting features for Brain model..."
            Rscript code/aging_model/step03-feature_selection.R \
                --feature_file "data/genes_set/brain_127genes.csv" \
                --data_file "data/ukb/olink_data_imputed_mean_drop30.csv" \
                --response_file "data/ukb/granular_age_april_07_2025.csv" \
                --output_dir "$output_features_dir/brain/" \
                --features_num_boruta $features_num_boruta \
                --n_boruta $n_boruta \
                --nrounds $nrounds \
                --features_num_optimized $features_num_optimized \
                --train_ratio $train_ratio \
                --seed $seed
        else
            echo "  Aging model feature selection already complete for Brain"
        fi

    done
done

re_train_models=FALSE
for features_num_boruta in ${features_num_boruta_list[@]}; do
    for features_num_optimized in ${features_num_optimized_list[@]}; do

        echo "Step 4: Running aging model analysis..."

        output_models_dir="results/aging_model_${features_num_optimized}_${features_num_boruta}/models"

        # Define feature files
        feature_file1="results/aging_model_${features_num_optimized}_${features_num_boruta}/features_selection/har/features_optimized.csv"
        feature_file2="results/aging_model_${features_num_optimized}_${features_num_boruta}/features_selection/hagr/features_optimized.csv"
        feature_file3="results/aging_model_${features_num_optimized}_${features_num_boruta}/features_selection/protage/features_optimized.csv"
        feature_file4="results/aging_model_${features_num_optimized}_${features_num_boruta}/features_selection/syngo/features_optimized.csv"
        feature_file5="results/aging_model_${features_num_optimized}_${features_num_boruta}/features_selection/brain/features_optimized.csv"
        FEATURE_FILES="$feature_file1,$feature_file2,$feature_file3,$feature_file4"

        FEATURE_FILES="$feature_file1,$feature_file2,$feature_file5"

        Rscript code/aging_model/step04-aging_model.R \
            --data_file "data/ukb/olink_data_imputed_mean_drop30.csv" \
            --response_file "data/ukb/granular_age_april_07_2025.csv" \
            --output_dir "$output_models_dir/" \
            --re_train_models $re_train_models \
            --train_ratio $train_ratio \
            --seed $seed \
            --feature_files "$FEATURE_FILES"

        echo "Aging model analysis completed!"

    done
done
