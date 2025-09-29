"""
📄 2. How to present predictions_pre_post.tsv
This file lists per-query changes in predictions. For the report:
Don’t dump the whole file (it’s huge).
Instead, summarise:
# of queries where prediction changed pre→post (and whether they got better or just shifted errors).
Examples: include a small table (5–10 rows) showing a few interesting cases where the filter fixed a misclassification or where errors persisted.
If the changes are rare (like in your snippet), that’s an important result: filtering didn’t substantially change predictions. That’s a strong, clear point.
👉 Example table for report (manually trimmed):
qseqid	gene_name	true_species	pred_pre	pred_post
133280_0_51	smc_1	Streptococcus_agalactiae	Staphylococcus_aureus	Staphylococcus_aureus
…	…	…	…	…

"""