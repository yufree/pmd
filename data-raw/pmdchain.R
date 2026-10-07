# Build data(pmdchain): a curated database of known multi-step reaction chains,
# each an ordered sequence of paired mass distances (PMDs) for known-pathway
# screening via getchainseq(list, db = pmdchain). PMD per step is the EXACT
# monoisotopic mass of its CHNOPS element
# delta (never hand-typed), matching the sda/breadth data-raw convention. A
# per-chain `keggcount` records how many ordered, compound-connected KEGG
# reaction paths share the chain's PMD signature -- an observable degeneracy /
# specificity score (lower = rarer = more specific), computed from data(keggrall).
suppressMessages({library(pmd); library(data.table)})

# exact monoisotopic atomic masses (C,H,N,O,P,S)
am <- c(C = 12, H = 1.0078250319, N = 14.0030740052,
        O = 15.9949146221, P = 30.97376151, S = 31.97207069)
massof <- function(d) sum(d * am[c("C","H","N","O","P","S")])

# single-step transformation building blocks: signed (C,H,N,O,P,S) element deltas
blk <- list(
        methylene  = c(1, 2,0,0,0,0),   # +CH2  (elongation / methylation)
        desat      = c(0,-2,0,0,0,0),   # -2H   (desaturation)
        oxidation  = c(0, 0,0,1,0,0),   # +O    (hydroxylation / oxidation)
        hydration  = c(0, 2,0,1,0,0),   # +H2O
        dehydrate  = c(0,-2,0,-1,0,0),  # -H2O
        decarboxyl = c(-1,0,0,-2,0,0),  # -CO2
        acetyl     = c(2, 2,0,1,0,0),   # +C2H2O (acetyl / ketide)
        phosphoryl = c(0, 1,0,3,1,0),   # +HPO3
        sulfate    = c(0, 0,0,3,0,1),   # +SO3
        hexose     = c(6,10,0,5,0,0),   # +C6H10O5 (glycosylation)
        glucuronide= c(6, 8,0,6,0,0),   # +C6H8O6  (glucuronidation)
        methyl_loss= c(-1,-2,0,0,0,0),  # -CH2 (demethylation)
        amination  = c(0, 3,1,0,0,0))   # +NH3 (amination)

# curated chains: ordered lists of block names + metadata
chains <- list(
        list(id="methylene_homolog", class="natural",
             name="Methylene homolog (elongation / methylation)",
             steps=c("methylene","methylene"),
             desc="repeated +CH2; alkyl chain elongation or sequential methylation"),
        list(id="desat_elong", class="natural",
             name="Desaturation-elongation (fatty acid)",
             steps=c("desat","methylene"),
             desc="canonical fatty-acid biosynthesis motif: desaturation then chain elongation"),
        list(id="lipid_edit", class="natural",
             name="Elongation-desaturation-oxidation (lipid editing)",
             steps=c("methylene","desat","oxidation"),
             desc="heterogeneous cascade of the three most recurrent metabolic edits"),
        list(id="seq_oxidation", class="natural",
             name="Sequential oxidation (hydroxylation)",
             steps=c("oxidation","oxidation"),
             desc="successive +O additions, e.g. steroid / bile-acid hydroxylation"),
        list(id="beta_oxidation", class="natural",
             name="Beta-oxidation intermediates",
             steps=c("desat","hydration","oxidation"),
             desc="dehydrogenation, hydration, oxidation of the fatty-acyl chain"),
        list(id="diglycoside", class="natural",
             name="Diglycosylation",
             steps=c("hexose","hexose"),
             desc="successive hexose additions (di-/oligo-glycoside formation)"),
        list(id="hydroxyl_glyc", class="natural",
             name="Hydroxylation-glycosylation",
             steps=c("oxidation","hexose"),
             desc="phase-I-like oxidation followed by glycosylation (plant detox)"),
        list(id="polyketide", class="natural",
             name="Polyketide ketide extension",
             steps=c("acetyl","acetyl"),
             desc="successive C2H2O ketide units"),
        list(id="phospho_cascade", class="natural",
             name="Sequential phosphorylation",
             steps=c("phosphoryl","phosphoryl"),
             desc="successive +HPO3 (e.g. mono- to di-/tri-phosphate)"),
        list(id="amino_methyl", class="natural",
             name="Amination-methylation",
             steps=c("amination","methylene"),
             desc="amino group addition followed by N-/O-methylation"),
        list(id="phase_glucuronide", class="xenobiotic",
             name="Phase I->II: hydroxylation-glucuronidation",
             steps=c("oxidation","glucuronide"),
             desc="cytochrome-P450 oxidation followed by glucuronide conjugation"),
        list(id="phase_sulfate", class="xenobiotic",
             name="Phase I->II: hydroxylation-sulfation",
             steps=c("oxidation","sulfate"),
             desc="oxidation followed by sulfate conjugation"),
        list(id="phase_methyl", class="xenobiotic",
             name="Phase I->II: hydroxylation-methylation",
             steps=c("oxidation","methylene"),
             desc="oxidation followed by O-/N-methylation"),
        list(id="demethyl_ox", class="xenobiotic",
             name="Oxidative demethylation-oxidation",
             steps=c("methyl_loss","oxidation"),
             desc="loss of a methyl followed by oxidation"),
        list(id="dihydroxy_dehydr", class="natural",
             name="Dihydroxylation-dehydration",
             steps=c("oxidation","oxidation","dehydrate"),
             desc="two hydroxylations then water elimination"))

# ---- assemble one row per step, exact PMD per element delta
rows <- rbindlist(lapply(chains, function(ch) {
        d  <- sapply(ch$steps, function(s) blk[[s]])           # 6 x nstep
        pmd <- apply(d, 2, function(v) abs(massof(v)))
        data.table(chain_id = ch$id, name = ch$name, class = ch$class,
                   nstep = length(ch$steps), step = seq_along(ch$steps),
                   transformation = ch$steps,
                   pmd = round(pmd, 4),
                   dC = d[1,], dH = d[2,], dN = d[3,], dO = d[4,], dP = d[5,], dS = d[6,],
                   description = ch$desc)
}))

# ---- specificity: count ordered, compound-connected KEGG reaction paths whose
#      successive |PMD| match each chain (data(keggrall)), within 0.002 Da.
data(keggrall, package = "pmd")
K <- as.data.table(keggrall)[!is.na(formula1) & !is.na(formula2) &
                             formula1 != "" & formula2 != ""]
K[, p := round(pmd, 4)]
matchstep <- function(dt_from, target, tol = 0.002)
        dt_from[abs(p - target) <= tol]
keggcount <- sapply(chains, function(ch) {
        pm <- rows[chain_id == ch$id][order(step), pmd]
        # paths start as reactions matching step 1
        cur <- matchstep(K, pm[1])[, .(end = formula2, n = 1L)]
        for (s in pm[-1]) {
                nxt <- matchstep(K, s)
                cur <- merge(cur, nxt[, .(formula1, formula2)],
                             by.x = "end", by.y = "formula1", allow.cartesian = TRUE)
                if (!nrow(cur)) break
                cur <- cur[, .(end = formula2, n = 1L)]
        }
        nrow(cur)
})
names(keggcount) <- sapply(chains, `[[`, "id")
rows[, keggcount := keggcount[chain_id]]

pmdchain <- as.data.frame(rows)
save(pmdchain, file = "data/pmdchain.rda", compress = "xz")
cat("pmdchain built:", length(unique(pmdchain$chain_id)), "chains,",
    nrow(pmdchain), "steps\n")
print(unique(pmdchain[, c("chain_id","nstep","class","keggcount")]), row.names = FALSE)
