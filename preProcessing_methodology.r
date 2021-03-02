# Required packages
require(phyloseq)
require(FCBF)

setwd('~/git/MIE_Metagenomics')
source('utils.r')

## --- Create Phyloseq object --- ##
# Defining input and outdir directories
inputDir = 'data/raw/'
outDir = 'data/processed/'

# Defining biom, sample and tree files (100nt) --> trimmed sequencing
biom_file = paste0(inputDir, 'otu_table_json.biom')
tree_file = paste0(inputDir, '97_otus.tree')
sampleDir = paste0(inputDir, 'clinical_SelectFeatures_unique_patients_fecal.rds')

# Reading the three files into phyloseq compatible objects
print('Importing OTU table ... ')
otu = import_biom(biom_file)
print('Reading tree file ...')
tree = read_tree(tree_file)
print('Reading sample data ...')

sample = readRDS(sampleDir)
sample_cohort = sample_data(sample)

phyloseq = merge_phyloseq(otu, tree, sample_cohort)
saveRDS(phyloseq, file = paste0(outDir, 'phyloseq.rds'))


## --- Glomerating, outlier analysis, labeling and balancing --- ##
# Load phyloseq object
filePath = paste0(outDir, 'phyloseq.rds')
taxLevel = 'Rank6' # genera

pseq = readRDS(filePath)

print(paste0('Glomerating OTUs by ', taxLevel, ' ...'))
pseqGlom = tax_glom(pseq, taxLevel)
pats = isoForest(pseqGlom)

pseq2 = prune_samples(pats, pseqGlom)

targets = c('DIET_TYPE')
for (i in seq_along(targets)) {

  res = labeling(pseq2, targets[i])
  saveRDS(res, file = paste0(outDir, 'phyloseq_', targets[i], '.rds'))

}


## --- Split --- ##
set.seed(1993)
pseq = readRDS(paste0(outDir), 'phyloseq_COUNTRY.rds')
smp_size = floor(0.85 * nsamples(pseq))
train_ind = sample(seq_len(nsamples(pseq)), size = smp_size)
trainPats = sample_names(pseq)[train_ind]
testPats = sample_names(pseq)[-train_ind]

train = subset_samples(pseq, sample_names(pseq) %in% trainPats)
test = subset_samples(pseq, sample_names(pseq) %in% testPats)

saveRDS(test, file = paste0(outDir, 'test_COUNTRY.rds'))
print('Done split!'))
  

## --- Feature Selectino (with FCBF) --- ##
otu = as.data.frame(t(otu_table(train)))
otu = apply(otu, 2, function(x) log2(x + 1))

otu = cbind.data.frame(otu, target = get_variable(data, 'COUNTRY')
names(otu) = make.names(names(otu))

fcbf = fast.cor.FS(otu, thresh = 0.0025)

saveRDS(fcbf, file = paste0(outDir, 'fcbf_COUNTRY'))
print('Done FCBF!')
  
