# Check of the bibliographic data

Record of the check of every reference cited in DESCRIPTION, README.Rmd, the vignettes
and the help pages. Each entry was compared with the first page of the original article
and with Crossref (api.crossref.org) or PubMed. Checked on 2026-10-05, except Berger and
Boos (1994), Boschloo (1970) and Fay and Hunsberger (2021), which were checked on
2026-10-06.

| Reference | Data used in the package | Original | Crossref or PubMed | Result |
|---|---|---|---|---|
| Blackwelder (1982) | Controlled Clinical Trials 3(4), 345-353, doi:10.1016/0197-2456(82)90024-1 | PDF p. 345: Controlled Clinical Trials 3:345-353 (1982), PII 0197-2456(82)90024-1 | Crossref: volume 3, issue 4, pages 345-353, December 1982, same DOI | Agrees |
| Farrington and Manning (1990) | Statistics in Medicine 9(12), 1447-1454, doi:10.1002/sim.4780091208 | PDF p. 1447: Statistics in Medicine, Vol. 9, 1447-1454 (1990) | Crossref: volume 9, issue 12, pages 1447-1454, December 1990 | Agrees |
| Friede and Kieser (2004) | Pharmaceutical Statistics 3(4), 269-279, doi:10.1002/pst.140 | PDF p. 269: Pharmaceut. Statist. 2004; 3: 269-279, DOI 10.1002/pst.140 | Crossref: volume 3, issue 4, pages 269-279, October 2004 | Agrees |
| Friede, Mitchell and Mueller-Velten (2007) | Biometrical Journal 49(6), 903-916, doi:10.1002/bimj.200610373 | PDF p. 903: Biometrical Journal 49 (2007) 6, 903-916, DOI 10.1002/bimj.200610373 | Crossref: volume 49, issue 6, pages 903-916, print 2007 | Agrees |
| Kieser and Friede (2000) | Statistics in Medicine 19(7), 901-911 | PDF p. 901: Statist. Med. 2000; 19:901-911 (the issue is not printed) | Crossref search did not find the article. PubMed (PMID 10750058): volume 19, issue 7, pages 901-911, DOI 10.1002/(SICI)1097-0258(20000415)19:7<901::AID-SIM405>3.0.CO;2-L | Agrees |
| Kieser (2020), Chapter 21 | Springer, Cham, doi:10.1007/978-3-030-49528-2_21 | Chapter PDF p. 225: Chapter 21, (c) Springer Nature Switzerland AG 2020, https://doi.org/10.1007/978-3-030-49528-2_21 | Crossref not checked (HTTP 429, twice) | Agrees with the original |
| Mehrotra, Chan and Berger (2003) | Biometrics 59(2), 441-450, doi:10.1111/1541-0420.00051 | PDF p. 441: Biometrics 59, 441-450, June 2003 | Crossref: volume 59, issue 2, pages 441-450, print June 2003 | Agrees |
| Fay and Hunsberger (2021) | Statistics Surveys 15, 72-110, doi:10.1214/21-SS131 | PDF p. 72: Statistics Surveys, Vol. 15 (2021) 72-110, ISSN 1935-7516, https://doi.org/10.1214/21-SS131 | Crossref: volume 15, no issue, no pages, print 2021, same DOI. The DOI resolves at doi.org to the article at projecteuclid.org | Agrees |
| Berger and Boos (1994) | Journal of the American Statistical Association 89, 1012-1016, doi:10.1080/01621459.1994.10476836 | PDF p. 1012: Journal of the American Statistical Association, September 1994, Vol. 89, No. 427, Theory and Methods; the article ends on p. 1016. The cover page of the publisher gives 89:427, 1012-1016 and the same DOI | Crossref not checked (HTTP 429, twice on 2026-10-06). The DOI resolves at doi.org to the article at tandfonline.com | Agrees |
| Boschloo (1970) | Statistica Neerlandica 24(1), 1-9, doi:10.1111/j.1467-9574.1970.tb00104.x | PDF p. 1: Statistica Neerlandica 24 (1970) nr. 1; the references end on p. 9 | Crossref: volume 24, issue 1, pages 1-9, print March 1970, DOI 10.1111/j.1467-9574.1970.tb00104.x | Agrees |

Berger and Boos (1994) has been cited in DESCRIPTION, README.md and the vignette
`bbssr-statistical-methods` since version 2.0.0, and Boschloo (1970) in README.md and the
same vignette. After the check, both are also cited in the help page of `BinaryRR()`,
Boschloo (1970) in DESCRIPTION and in the text of the two vignettes
`bbssr-statistical-methods` and `bbssr-validation`, and the numbers of Boschloo (1970)
are recomputed by `inst/reproduce/reproduce-published.R`. Besides the bibliographic data,
the statements about the two articles were compared with the originals.

- Berger and Boos (1994), Sections 1 and 2: for a 1 - beta confidence set C_beta of the
  nuisance parameter under the null hypothesis, p_beta = sup over C_beta of p(theta), plus
  beta, is a valid p-value (Lemma), and beta is chosen small, "such as .001 or .0001".
  Example 2 applies the procedure to the 2 x 2 table with a .999 confidence interval for
  the common proportion. The vignette describes the same procedure and writes gamma for
  beta, after the argument `bb.gamma`.
- Boschloo (1970), Sections 3 and 6: Fisher's test is used at a raised conditional level
  gamma, chosen as high as possible such that the unconditional level does not exceed
  alpha for any value of the common probability. Rejecting when the Fisher p-value is at
  most gamma is the same as rejecting when the unconditional p-value with the Fisher
  p-value as ordering statistic is at most alpha, which is how the package defines the
  test. In the two-sided test each of the two parts of the critical region has the
  conditional level gamma / 2, and the example of Section 6 rejects because
  0.0565 < 0.114 / 2 (printed as ".565 < .144/2"). This is the `'central'` convention of
  `tsmethod`. An independent computation with Python (exact rational conditional
  p-values) gives the numbers of Sections 2, 3 and 6, including the two-sided raised
  level 0.114, and columns I and II of the table of Section 4 except one entry: the power
  of Fisher's test at alpha = .05, p1 = .6 and p2 = .1 is 0.84515004, which rounds to
  .8452 against the published .8451. Column III, the randomized test, is not part of the
  package.

The two-sided convention `tsmethod = 'blaker'` was added on 2026-10-06. Formula (2) of
Mehrotra, Chan and Berger (2003) is cited for it in the help page of `BinaryRR()`,
README.Rmd, NEWS.md and the vignettes `bbssr-statistical-methods` and `bbssr-validation`.
Fay and Hunsberger (2021) is cited for its mid-p version and its example in the help page,
NEWS.md and the two vignettes, and as a reproduced source in README.Rmd. The statements
were compared with the originals.

- Mehrotra, Chan and Berger (2003), Section 2: the two-sided p-value of Fisher's test is
  formula (2), the sum of the conditional probabilities of the tables i with
  g(i, t) <= g(x1, t), where g(x, t) is the smaller of the two one-sided tail
  probabilities, "following Blaker (2000) and Agresti and Min (2001)". The Boschloo test
  of formula (13) orders the outcomes by this p-value. Blaker (2000) and Agresti and Min
  (2001) are not cited by the package, since their originals have not been checked.
- Mehrotra, Chan and Berger (2003), Section 3.1 and Tables 1 and 3: the example values
  (CP interval (0.0080, 0.0826) and the nine p-values) are reproduced with formula (2),
  and equally with the minlike convention, which gives the same p-values there; the
  central convention does not reproduce them. In Tables 1 and 3, wherever formula (2) and
  the minlike convention give different rounded values in the configurations with groups
  of unequal size, the published values of F, B and B* are those of minlike, and the
  percentages agree with rounding to two decimals followed by rounding to one decimal.
  This is recorded in
  `inst/reproduce/reproduce-published.R` and checked by
  `tools/reference/check_mehrotra_2003.py`.
- Fay and Hunsberger (2021), Section 8 and Table 1: Blaker's two-sided ordering function
  T_B(x, beta) is the probability of the tables whose gamma, the smaller of the two tail
  probabilities, does not exceed that of x, and for 8 of 14 against 1 of 7 the p-values
  are 0.159 (Fisher-Irwin, the minlike convention), 0.087 (Blaker) and 0.157 (central).
  The values of T_B(x, 1) in Table 1 are recomputed by the reproduction script.
- Fay and Hunsberger (2021), Section 9: "the mid p-value is 0.5 times the probability of
  equality plus the probability of more extreme". The definition is stated in general
  terms and illustrated with the one-sided conditional p-value of equation (9.1), and
  Table 2 lists the mid-p adjustment as applying to any method. The package applies this
  definition to the orderings of the minlike and blaker conventions.
