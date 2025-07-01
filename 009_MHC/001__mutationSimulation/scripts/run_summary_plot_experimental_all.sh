#!/bin/bash

# Make sure scripts exist
SUMMARY_SCRIPT="scripts/summarize_simulation_results_v2.py"
PLOT_SCRIPT="scripts/plot_signature_summary.py"



# Loop through each folder in results/, excluding those that start with SBS
for dir in results/*/; do
    BASENAME=$(basename "$dir")
    
    # Skip if it starts with SBS
    if [[ $BASENAME == SBS* ]]; then
        continue
    fi

    echo "🔄 Processing $BASENAME"

    # Set file paths
    SUMMARY_OUT="${dir}${BASENAME}_summary_all.tsv"
    PLOT_OUT="${dir}${BASENAME}_summary_plots.pdf"

    # Run summary script
    python "$SUMMARY_SCRIPT" \
        --results-dir "$dir" \
        --output "$SUMMARY_OUT"

    python "$PLOT_SCRIPT" \
        --input "$SUMMARY_OUT" \
        --output "$PLOT_OUT"

    echo "✅ Done with $BASENAME"
done

