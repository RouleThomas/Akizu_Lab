#!/bin/bash

# Make sure scripts exist
SUMMARY_SCRIPT="scripts/summarize_simulation_results_v3.py"
PLOT_SCRIPT="scripts/plot_signature_summary.py"

# Loop through each results/* folder
for dir in results/*/; do
    # Remove trailing slash and extract base name
    SIGNATURE=$(basename "$dir")
    echo "🔄 Processing $SIGNATURE"

    # Set file paths
    SUMMARY_OUT="${dir}${SIGNATURE}_summary_all_v3.tsv"
    PLOT_OUT="${dir}${SIGNATURE}_summary_plots_v3.pdf"

    # Run summary script
    python "$SUMMARY_SCRIPT" \
        --results-dir "$dir" \
        --output "$SUMMARY_OUT"

    # Run plot script
    python "$PLOT_SCRIPT" \
        --input "$SUMMARY_OUT" \
        --output "$PLOT_OUT"

    echo "✅ Done with $SIGNATURE"
done


