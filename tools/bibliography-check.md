# Check of the bibliographic data

Record of the check of every reference cited in DESCRIPTION, README.Rmd, the vignettes
and the help pages. Each entry was compared with the first page of the original article
and with Crossref (api.crossref.org) or PubMed. Checked on 2026-10-05.

| Reference | Data used in the package | Original | Crossref or PubMed | Result |
|---|---|---|---|---|
| Blackwelder (1982) | Controlled Clinical Trials 3(4), 345-353, doi:10.1016/0197-2456(82)90024-1 | PDF p. 345: Controlled Clinical Trials 3:345-353 (1982), PII 0197-2456(82)90024-1 | Crossref: volume 3, issue 4, pages 345-353, December 1982, same DOI | Agrees |
| Farrington and Manning (1990) | Statistics in Medicine 9(12), 1447-1454, doi:10.1002/sim.4780091208 | PDF p. 1447: Statistics in Medicine, Vol. 9, 1447-1454 (1990) | Crossref: volume 9, issue 12, pages 1447-1454, December 1990 | Agrees |
| Friede and Kieser (2004) | Pharmaceutical Statistics 3(4), 269-279, doi:10.1002/pst.140 | PDF p. 269: Pharmaceut. Statist. 2004; 3: 269-279, DOI 10.1002/pst.140 | Crossref: volume 3, issue 4, pages 269-279, October 2004 | Agrees |
| Friede, Mitchell and Mueller-Velten (2007) | Biometrical Journal 49(6), 903-916, doi:10.1002/bimj.200610373 | PDF p. 903: Biometrical Journal 49 (2007) 6, 903-916, DOI 10.1002/bimj.200610373 | Crossref: volume 49, issue 6, pages 903-916, print 2007 | Agrees |
| Kieser and Friede (2000) | Statistics in Medicine 19(7), 901-911 | PDF p. 901: Statist. Med. 2000; 19:901-911 (the issue is not printed) | Crossref search did not find the article. PubMed (PMID 10750058): volume 19, issue 7, pages 901-911, DOI 10.1002/(SICI)1097-0258(20000415)19:7<901::AID-SIM405>3.0.CO;2-L | Agrees |
| Kieser (2020), Chapter 21 | Springer, Cham, doi:10.1007/978-3-030-49528-2_21 | Chapter PDF p. 225: Chapter 21, (c) Springer Nature Switzerland AG 2020, https://doi.org/10.1007/978-3-030-49528-2_21 | Crossref not checked (HTTP 429, twice) | Agrees with the original |
| Mehrotra, Chan and Berger (2003) | Biometrics 59(2), 441-450, doi:10.1111/1541-0420.00051 | PDF p. 441: Biometrics 59, 441-450, June 2003 | Crossref: volume 59, issue 2, pages 441-450, print June 2003 | Agrees |
| Berger and Boos (1994) | Journal of the American Statistical Association 89, 1012-1016, doi:10.1080/01621459.1994.10476836 | Not available. The reference list of Mehrotra, Chan and Berger (2003) gives Journal of the American Statistical Association 89, 1012-1016 | Crossref not checked (HTTP 429, twice) | Not verified with the original; the PDF is needed |
| Boschloo (1970) | Statistica Neerlandica 24, 1-9 | Not available | Not checked | Not verified with the original; the PDF is needed |

Berger and Boos (1994) has been cited in DESCRIPTION, README.md and the vignette
`bbssr-statistical-methods` since version 2.0.0, and Boschloo (1970) in README.md and the
same vignette. The two references are kept there until the originals can be checked, are
not added anywhere else, and are removed before the release if the originals cannot be
obtained.
