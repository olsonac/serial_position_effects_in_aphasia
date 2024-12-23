Table 2 - length and frequency - words correct
- `./output/[ppt]/tables/[ppt]_repetition_byword_len_freq_model_resuls.csv`

Table 3 - length and position  - positions correct
- `./output/[ppt]/reports/[ppt]_repetition_analysis_report.pdf`
- for % change with length, search for `mean change in probability for each additional length:` (multiply by 100 to get %)
- for % change with position, search for `mean change in probability for each additional position (excluding U-shape)` (multiply by 100 to get %)
- for % change with position in final upward portion of the curve, search for `Average upward change after U minimum` (multiply by 100 to get %)
- for range of % preserved, search for `Min/max preserved range:`

Figure 2 - simulated data 
- plots are .tif files in ./src/simulation
- `theory_...` plots show theoretical (not sampled) values
- `simulated_data...` plots show sampled data values generated from theoretical models and corpus characteristics

Table 4 - Relative importance of factors - Analysis of deviance
-  Percentage values are based on Nagelkerke values from `./output/[ppt]/tables/[ppt]_repetition_dominance_analysis.csv` and the percentages themselves are found in `./output/[ppt]/tables/[ppt]_repetition_nagelkerke_contributions_table.csv`.  The % variance accounted for is found in `./output/[ppt]/tables/[ppt]_repetition_see_results_table.csv` in the column `p_accounted_for`

Figure 3 - Plots of p(preserved) predicted by previous errors/previous correct
- previous error plots are in `./output/[ppt]/fig/[ppt]_repetitionPrevErrPlot_ObsPredME.tif` (ME = main effect)
- previous correct plots are in `./output/[ppt]/fig/[ppt]_repetitionPrevCorPlot_ObsPredME.tif`

Figure 4 - plot of % of segments that were coda or satellite syllable components across word position
- plots in `./output/[ppt]/fig/[ppt]_repetition_pos_of_syll_complexities_in_stimuli.png`
- These are based on **stimuli** so are identical in different participant records.  The plot can be obtained from any of the participants.

Supplementary materials D - regression coefficients for logistic regression models for words correct based on length and log(frequency)
- `./output/[ppt]/tables/[ppt]_repetition_byword_len_freq_model_resuls.csv`

Supplementary materials E - Logistic regression coefficients for length X position models
- `./output/[ppt]/tables/[ppt]_repetition_len_pos_model_summary.csv
-`
Supplementary materials F - regression coefficients for models based on previous correct/errors, position, length and log(frequency)
- `./output/[ppt]/tables/[ppt]_repetition_main_effects_plus_one_len_freq_summary.csv`

Supplementary materials F - plots of models based on previous errors, previous correct, full final model
- plot for previous errors is found in `./output/[ppt]/fig/[ppt]_repetitionPrevErrPlot_ObsPredME.tif`
- plot for previous correct is found in `./output/[ppt]/fig/[ppt]_repetitionPrevCorPlot_ObsPredME.tif`
- plot for the full final model is found in `./output/[ppt]/fig/[ppt]_repetition_ObsPred_BestFullModel_prev_correct.tif` or in `./output/[ppt]/fig/[ppt]_reptition_ObsPred_BestFullModel_prev_error.tif`

