#!/bin/bash

# Make sure scripts exist
SUMMARY_SCRIPT="scripts/summarize_simulation_results.py"
PLOT_SCRIPT="scripts/plot_signature_summary.py"

# Loop through each results_contexts/* folder
for dir in results_contexts/*/; do
    # Remove trailing slash and extract base name
    SIGNATURE=$(basename "$dir")
    echo "🔄 Processing $SIGNATURE"

    # Set file paths
    SUMMARY_OUT="${dir}${SIGNATURE}_summary_all.tsv"
    PLOT_OUT="${dir}${SIGNATURE}_summary_plots.pdf"

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


