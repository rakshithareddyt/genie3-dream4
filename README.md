# GENIE3 on DREAM4 In Silico Size-100: parameter sweep and comparison with correlation

Course: Systems Biology & Network Modeling (CB 5313) · Roll no. SE23UCAB017 · Topic 6: Directed GRN inference (GENIE3)

## Question

Does GENIE3 (Huynh-Thu et al., 2010, PLoS ONE 5(9): e12776) recover known regulatory edges better than plain co-expression (\|Pearson r\|, from this course)? How sensitive is it to its parameters (K, number of trees, Random Forest vs Extra-Trees)? Can it recover edge direction?

## Data

DREAM4 In Silico Size-100 networks 1-5 (simulated expression data with gold-standard edges), obtained from the Bioconductor `DREAM4` package git repository. Raw data is not committed to this repository; `download_data.py` fetches it. I use pooled steady-state profiles (wild-type + knockouts + knockdowns, 201 samples x 100 genes).

### Why I did not use the multifactorial file

The original GENIE3 paper uses the DREAM4 multifactorial data, so I first tried those files from the Bioconductor `DREAM4` package. They did not work in my setup. On network 4, GENIE3 on the multifactorial file gave an AUPR of 0.024 and an AUROC of 0.497, which is about the same as chance (AUPR 0.021, AUROC 0.5). Pearson correlation was also at chance.

To find out whether the problem was my code or the data, I ran three checks:

1. I scored the multifactorial expression file of every network against the gold standard of every network (a 5x5 grid). All values were between about 0.46 and 0.54, with no clear diagonal, so no file matched its own gold standard better than the others.

2. I checked whether the file was transposed or its columns were mislabelled. Transposing did not help (AUROC 0.45 to 0.54), and the column means of the multifactorial file did not follow the wild-type expression values (correlation 0.01 to 0.16).

3. As a control, I used the knockout files from the same package. These matched their own gold standards clearly (AUROC on the diagonal: 0.86, 0.75, 0.78, 0.81, 0.72, compared with about 0.46 to 0.55 off the diagonal).

So the gold standards and the knockout data are correct, and the problem is specific to the multifactorial file in this copy of the data. I did not find out why. I did not test whether the file is mislabelled, corrupted or different from the original challenge file, so I make no claim about the cause.

Instead I pooled the wild-type, knockout and knockdown profiles of each network (201 samples x 100 genes) and used those. This is not the paper's setting, so my results are not a reproduction of the paper's DREAM4 results. Also, knockouts are interventions, but GENIE3 treats them as ordinary samples and ignores which gene was perturbed.

## How to run

Requires Python 3, git and an internet connection.

```bash
pip install -r requirements.txt
python download_data.py          # fetches data into ./data
python sweep.py --net 4          # K / trees / RF-vs-ET sweep, seed repeats, Pearson baseline (about 10-25 min)
python sweep.py --robust         # default setting on all 5 networks (about 5 min)
python directionality.py         # is edge direction recovered? (about 5 min)
python sweep.py --net 4 --quick  # fast smoke test (about 1 min)
```

Outputs are written to the `results/` folder.

## Code
- `genie3.py`: GENIE3 implemented with scikit-learn (one tree ensemble per target gene; raw variance-reduction importance as the edge score) and the Pearson baseline
- `evaluate.py`: AUPR and AUROC against the gold standard over all ordered gene pairs
- `sweep.py`: parameter sweep, baseline comparison, seed repeats, five-network run
- `directionality.py`: tests whether the true direction of each edge is scored higher
- `download_data.py`: downloads the data

## Results

Overall, GENIE3 was able to identify some meaningful regulatory interactions from the expression data, but the performance depended on the model and parameters used. I evaluated the predictions using AUPR and compared them with a Pearson correlation baseline. On Network 4, the Random Forest model with 100 trees gave an AUPR of 0.133 when using K ≈ 10 (the square root of 99), which increased to 0.186 when K = 99. Extra-Trees performed even better on this network, reaching an AUPR of 0.223 with K = 99, while the Pearson baseline was 0.089.

Across all five networks, GENIE3 had a higher AUPR than Pearson on 4 of 5 networks (mean AUPR 0.125 compared with 0.094). Network 3 was the exception: GENIE3 had an AUPR of 0.121, while Pearson reached 0.131. Network 4 was the most favourable case for GENIE3, so looking only at its results would overstate the overall gain.

I also checked variation between random seeds for one setting. The AUPR varied by about 0.015 across the three seeds, so small differences between settings, such as using 100 versus 400 trees, should not be treated as reliable improvements.

GENIE3 also gives a direction to predicted regulatory links. Across the five networks, it oriented 538 of the 978 true edges correctly (about 55%), with the results ranging from 51% to 63% across networks. Pearson correlation cannot determine direction because the two directions for a gene pair receive the same score, so they are effectively tied. The paper reported a lower direction error rate, about 26% at full recall (roughly 74% correct), but that result used multifactorial data, K = p − 1 and 1000 trees. My experiment used a different data setup and settings, and I did not test why the direction result was lower.

## Limitations

One limitation is that I did not use the same dataset setup as the original DREAM4 multifactorial experiment, so my results cannot be directly compared with the headline results reported in the paper. The paper's multifactorial experiment used steady-state expression profiles generated from simultaneous perturbations of all genes, whereas my experiment used the pooled dataset with different experimental conditions.

The parameter and model settings were compared on Network 4 only, so I cannot assume that the best setting there would also be best for the other four networks. I measured variation between random seeds for only one setting, so I do not know whether the same amount of variation applies to the other settings. Also, the pooled data contain knockout profiles, which are interventions, but GENIE3 treats them as ordinary samples and does not account for which gene was perturbed.

Finally, an AUPR value alone does not tell us whether every predicted interaction is biologically meaningful. GENIE3 gives a ranking of possible regulatory links, so choosing which interactions should actually be included in a final network is another challenge.

## AI-use disclosure

Claude assisted with debugging a data problem, and helping with the code in this repository. I ran all analyses. The synopsis, the annotated paper and the interpretation of the results are my own.
