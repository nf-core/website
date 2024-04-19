---
title: meta map and ext properties
subtitle: What are the meta map and ext properties in nf-core components?
---

In nf-core DSL2 pipelines, to add sample-specific information and metadata that is carried throughout the pipeline, we use a meta variable. This avoids the need to create separate channels for each new characteristic.
The meta variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the meta variable.
The `meta map` is a [groovy map](https://www.tutorialspoint.com/groovy/groovy_maps.htm), which is like a python dictionary, as shown below:

```groovy
[id: 'test', single_end: false]
```

Thus, the information can be accessed within processes and `module.conf` files with the key i.e. `meta.id`

The meta variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the meta variable.

```groovy
input:
tuple val(meta), path(reads)
```

This pattern doesn't work out of the box with [fromFilePairs](https://www.nextflow.io/docs/edge/channel.html#fromfilepairs)

The difference between the two:

```groovy
// fromFilePairs
filepairs = [
    SRR493366,
    [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]
]

// meta map
meta_map = [
    [id: 'test', single_end: false], // meta map
    [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]
]
```

As you can see the difference, they are both [groovy lists](https://www.tutorialspoint.com/groovy/groovy_lists.htm).
However, the filepairs just has a `val` that is a string, where as the `meta_map` the first value in the list, is a [groovy map](https://www.tutorialspoint.com/groovy/groovy_maps.htm), which is like a python dictionary.
The only required value is `meta.id` for most of the modules, however, they usually contain fields like `meta.single_end` and `meta.strandedness`

### Common patterns

The `meta map` is generated with [create_fastq_channel function in the input_check subworkflow](https://github.com/nf-core/rnaseq/blob/587c61b441c5e00bd3201317d48b95a82afe6aaa/subworkflows/local/input_check.nf#L23-L45) of most nf-core pipelines. Where the meta information is easily extracted from a samplesheet that contains the input file paths.

### Generating a `meta map` from file pairs

Sometimes you want to use nf-core modules in small scripts. You don't want to make a samplesheet, or maintain a bunch of validation.
For instance, here's an example script to run fastqc

```groovy
nextflow.enable.dsl = 2

params.input = "*.fastq.gz"

include { FASTQC } from "./modules/nf-core/modules/fastqc/main"

workflow {
    ch_fastq = Channel.fromFilePairs(params.input, size: -1)
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
        }

    FASTQC ( ch_fastq )
}
```

### Sorting samples by groups

```groovy
ch_genome_bam.map {
    meta, bam ->
    fmeta = meta.findAll { it.key != 'read_group' }
    fmeta.id = fmeta.id.split('_')[0..-2].join('_')
    [ fmeta, bam ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { ch_sort_bam }
```

### Combining channel on meta subset

Sometimes it is necessary to combine multiple channels based on a subset of the meta maps.
Unfortunately this is not yet supported as the argument `by` isn't a closure in `.combine()` and `.join()` and it probably won't ([Nextflow issue #3175](https://github.com/nextflow-io/nextflow/issues/3175)).

To bypass this restriction one of the solution is to create a new map with only the necessary keys and make the junction on it. Here is an example:

```groovy
ch_input = [[["id":"Ind1","ref":"RefA"],"file1"],[["id":"Ind2","ref":"RefB"],"file2"]]
ch_ref   = [[["ref":"RefA"],"fileA"],[["ref":"RefB"],"fileB"]]

ch_join  = ch_input
            .map{metaIR, file -> [metaIR.subMap(["ref"]), metaIR, file]}
            .combine(chr_ref)
            .map{metaR, metaIR, file, ref -> [metaIR, file, ref]}
```

### Modify the meta map

There is multiple ways to modify the meta map.
Here are some examples:

```groovy
// Add to map - adding two maps makes a new Map object
ch.map { meta, files -> [ meta + [ single_end: files instanceof Path ], files ] }

// Remove certain keys (and their entries) from a map
ch.map { meta, files -> [ meta.subMap( ['id','rg'] ), files ] }
  // OR by specifying what not to include
ch.map { meta, files -> [ meta.findAll { ! it.key in ['single_end'] }, files ] }

// Split a map - use both methods of removing keys ( there is a split method for Maps, but the results are not Maps )
ch.map { meta, files -> def keyset = ['id', 'read_group']; [ meta.subMap(keyset), meta.findAll { ! it.key in keyset },  files ] }
```

### Conclusion

As you can see the `meta map` is a quite flexible way for storing meta data in channels. Feel free to add whatever other key-value pairs your pipeline may need to it. We're looking to add [Custom objects](https://github.com/nf-core/modules/issues/1338) which will lock down the usage a bit more.

## Advanced pattern

### Multimaping

It is possible with `multiMap` to split a channel in to and to call them separately afterwards.

```groovy
ch_input = reads.combine(db).multiMap{ it ->
   reads: it[0]
   db: it[1]
}
MODULE(ch_input.reads, ch_input.db)
```

### Adding additional information to the meta map

It is possible to combine a input channel with a set of parameters as follows:

```groovy
ch_input.flatMap { meta, filetype ->
    [300, 500, 1000].collect {
      def new_meta = meta.clone()
      new_meta.window_size = it
      [ new_meta, filetype]
    }
}
```

You can also combine this technique with others for more processing:

```groovy
workflow {

    input = [
        [
            [ patient: 'sample', sample: 'test', id: 'test' ],
            file ("chr21_23355001-46709983.bed")
        ],
        [
            [ patient: 'sample', sample: 'test', id: 'test' ],
            file ("chr21_2-23354000.bed")
        ],
        [
            [ patient: 'sample2', sample: 'test5', id: 'test' ],
            file ("chr21_23355001-46709983.bed")
        ],
        [
            [ patient: 'sample2', sample: 'test5', id: 'test' ],
            file ("chr21_2-23354000.bed")
        ]
    ]
    Channel.fromList ( input )
        .map { meta, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? new_meta.sample + "_" + intervals.baseName : new_meta.sample
            intervals = intervals.baseName != "no_intervals" ? intervals : []
            [new_meta, intervals]
        }.view { meta, intervals -> meta.id }
}
```
