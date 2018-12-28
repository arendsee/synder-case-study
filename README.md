[![unstable](http://badges.github.io/stability-badges/dist/unstable.svg)](http://github.com/badges/stability-badges)

# A case study for `synder`

This repo contains the code needed to generate the tables and figures used in
the synder paper and its supplementary files. And more.

See [main synder page here](https://github.com/arendsee/synder)

# Data Retrieval

All data needed to run this program is available on DataHub. To retrieve this
data, first install the data CLI tool by following the instructions available
[here](https://datahub.io/download). Then run the command:

``` sh
data get https://datahub.io/arendsee/synder-nfyc
```

This will produce a folder with the following structure:

```
arendsee/
└── synder-nfyc
    ├── archive
    │   ├── al.blast.tab
    │   ├── al.gff
    │   ├── at.gff
    │   ├── at-vs-al.tab
    │   ├── at-vs-br.tab
    │   ├── at-vs-cr.tab
    │   ├── at-vs-ep.tab
    │   ├── at-vs-es.tab
    │   ├── br.blast.tab
    │   ├── br.gff
    │   ├── cr.blast.tab
    │   ├── cr.gff
    │   ├── es.blast.tab
    │   ├── es.gff
    │   ├── nfyc.gff
    │   ├── nfyc-map.tab
    │   ├── tree.newick
    │   └── yc-vs-al.blast.tab
    ├── data
    │   └── validation_report.json
    └── datapackage.json
```
