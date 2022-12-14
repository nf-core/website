```bash
# Install nextflow
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

# Launch the RNAseq pipeline
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --output ./results/ \
    --genome GRCh37 \
    -profile docker

# Install nf-core tools
pip install nf-core

# List all nf-core pipelines and show available updates
nf-core list
```
