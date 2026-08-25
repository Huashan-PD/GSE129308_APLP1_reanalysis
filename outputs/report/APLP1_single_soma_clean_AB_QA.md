# APLP1 single-soma clean AB QA

1. Final figure contains only Panel A and Panel B: YES
2. Panel C table removed from figure: YES
3. Supplementary statistics table exported: YES
4. Panel B displays only Ex2, Ex7, Ex10: YES
5. Panel B excludes non-AD: YES
6. Panel B star source: fraction_FDR only
7. Expr. FDR used for Panel B stars: NO
8. If fraction_FDR >= 0.05 for all displayed subtypes, no stars shown: YES
   Displayed subtype fraction_FDR values: Ex2=0.127, Ex7=0.127, Ex10=0.367
9. Panel B y-axis spans 0 to 1 and starts at zero: YES
10. Caption markdown generated: YES
11. Output files generated:
   - /Users/ywb/Documents/aplp1/GSE129308_APLP1_reanalysis/outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars.pdf
   - /Users/ywb/Documents/aplp1/GSE129308_APLP1_reanalysis/outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars.png
   - /Users/ywb/Documents/aplp1/GSE129308_APLP1_reanalysis/outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars.svg

Star logic implemented in script:
   if fraction_FDR < 0.001: "***"
   else if fraction_FDR < 0.01: "**"
   else if fraction_FDR < 0.05: "*"
   else: no star

Forbidden star sources not used: Expr. FDR, expression_FDR, pseudobulk_FDR, MAST_FDR, MAST adjusted P, pseudobulk expression p value.
Prior expression-FDR-style star behavior is not used in this noExprStars figure; Panel B stars are generated exclusively from donor-level fraction_FDR.
