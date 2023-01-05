declare module 'astro:content' {
    export { z } from 'astro/zod';
    export type CollectionEntry<C extends keyof typeof entryMap> = typeof entryMap[C][keyof typeof entryMap[C]] &
        Render;

    type BaseCollectionConfig<S extends import('astro/zod').ZodRawShape> = {
        schema?: S;
        slug?: (entry: {
            id: CollectionEntry<keyof typeof entryMap>['id'];
            defaultSlug: string;
            collection: string;
            body: string;
            data: import('astro/zod').infer<import('astro/zod').ZodObject<S>>;
        }) => string | Promise<string>;
    };
    export function defineCollection<S extends import('astro/zod').ZodRawShape>(
        input: BaseCollectionConfig<S>
    ): BaseCollectionConfig<S>;

    export function getEntry<C extends keyof typeof entryMap, E extends keyof typeof entryMap[C]>(
        collection: C,
        entryKey: E
    ): Promise<typeof entryMap[C][E] & Render>;
    export function getCollection<C extends keyof typeof entryMap, E extends keyof typeof entryMap[C]>(
        collection: C,
        filter?: (data: typeof entryMap[C][E]) => boolean
    ): Promise<(typeof entryMap[C][E] & Render)[]>;

    type InferEntrySchema<C extends keyof typeof entryMap> = import('astro/zod').infer<
        import('astro/zod').ZodObject<Required<ContentConfig['collections'][C]>['schema']>
    >;

    type Render = {
        render(): Promise<{
            Content: import('astro').MarkdownInstance<{}>['Content'];
            headings: import('astro').MarkdownHeading[];
            injectedFrontmatter: Record<string, any>;
        }>;
    };

    const entryMap: {
        events: {
            '2018/hackathon-scilifelab-2018.md': {
                id: '2018/hackathon-scilifelab-2018.md';
                slug: '2018/hackathon-scilifelab-2018';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2018/nextflow-camp-2018.md': {
                id: '2018/nextflow-camp-2018.md';
                slug: '2018/nextflow-camp-2018';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2019/hackathon-qbic-2019.md': {
                id: '2019/hackathon-qbic-2019.md';
                slug: '2019/hackathon-qbic-2019';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2019/hackathon-scilifelab-2019.md': {
                id: '2019/hackathon-scilifelab-2019.md';
                slug: '2019/hackathon-scilifelab-2019';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2019/nextflow-camp-2019-community-updates.md': {
                id: '2019/nextflow-camp-2019-community-updates.md';
                slug: '2019/nextflow-camp-2019-community-updates';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2019/nextflow-camp-2019-tutorial.md': {
                id: '2019/nextflow-camp-2019-tutorial.md';
                slug: '2019/nextflow-camp-2019-tutorial';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/BovReg-CRG-workshop-2020.md': {
                id: '2020/BovReg-CRG-workshop-2020.md';
                slug: '2020/bovreg-crg-workshop-2020';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/aws-webinar.md': {
                id: '2020/aws-webinar.md';
                slug: '2020/aws-webinar';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/bosc-2020-demo.md': {
                id: '2020/bosc-2020-demo.md';
                slug: '2020/bosc-2020-demo';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/bosc-2020.md': {
                id: '2020/bosc-2020.md';
                slug: '2020/bosc-2020';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/eccb2020-workshop.md': {
                id: '2020/eccb2020-workshop.md';
                slug: '2020/eccb2020-workshop';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/elixir-proteomics.md': {
                id: '2020/elixir-proteomics.md';
                slug: '2020/elixir-proteomics';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/hackathon-francis-crick-2020.md': {
                id: '2020/hackathon-francis-crick-2020.md';
                slug: '2020/hackathon-francis-crick-2020';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/hackathon-july-2020.md': {
                id: '2020/hackathon-july-2020.md';
                slug: '2020/hackathon-july-2020';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/ifc-mexico-seminar.md': {
                id: '2020/ifc-mexico-seminar.md';
                slug: '2020/ifc-mexico-seminar';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/intro-to-eager.md': {
                id: '2020/intro-to-eager.md';
                slug: '2020/intro-to-eager';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2020/jobim-2020.md': {
                id: '2020/jobim-2020.md';
                slug: '2020/jobim-2020';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-1-nf-core-into.md': {
                id: '2021/bytesize-1-nf-core-into.md';
                slug: '2021/bytesize-1-nf-core-into';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-10-institutional-profiles.md': {
                id: '2021/bytesize-10-institutional-profiles.md';
                slug: '2021/bytesize-10-institutional-profiles';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-11-dev-envs-workflows.md': {
                id: '2021/bytesize-11-dev-envs-workflows.md';
                slug: '2021/bytesize-11-dev-envs-workflows';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-12-template-sync.md': {
                id: '2021/bytesize-12-template-sync.md';
                slug: '2021/bytesize-12-template-sync';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-13-tuning-pipeline-performance.md': {
                id: '2021/bytesize-13-tuning-pipeline-performance.md';
                slug: '2021/bytesize-13-tuning-pipeline-performance';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-14-graphic-design.md': {
                id: '2021/bytesize-14-graphic-design.md';
                slug: '2021/bytesize-14-graphic-design';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-15-pipeline-first-release.md': {
                id: '2021/bytesize-15-pipeline-first-release.md';
                slug: '2021/bytesize-15-pipeline-first-release';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-16-module-test-data.md': {
                id: '2021/bytesize-16-module-test-data.md';
                slug: '2021/bytesize-16-module-test-data';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-17-pytest-workflow.md': {
                id: '2021/bytesize-17-pytest-workflow.md';
                slug: '2021/bytesize-17-pytest-workflow';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-18-dev-envs-workflows.md': {
                id: '2021/bytesize-18-dev-envs-workflows.md';
                slug: '2021/bytesize-18-dev-envs-workflows';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-19-aws-megatests.md': {
                id: '2021/bytesize-19-aws-megatests.md';
                slug: '2021/bytesize-19-aws-megatests';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-2-configs.md': {
                id: '2021/bytesize-2-configs.md';
                slug: '2021/bytesize-2-configs';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-20-nextflow-tower.md': {
                id: '2021/bytesize-20-nextflow-tower.md';
                slug: '2021/bytesize-20-nextflow-tower';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-21-nf-core-mhcquant.md': {
                id: '2021/bytesize-21-nf-core-mhcquant.md';
                slug: '2021/bytesize-21-nf-core-mhcquant';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-22-nf-core-eager.md': {
                id: '2021/bytesize-22-nf-core-eager.md';
                slug: '2021/bytesize-22-nf-core-eager';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-23-nf-core-hic.md': {
                id: '2021/bytesize-23-nf-core-hic.md';
                slug: '2021/bytesize-23-nf-core-hic';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-24-dsl2-pipeline-starter.md': {
                id: '2021/bytesize-24-dsl2-pipeline-starter.md';
                slug: '2021/bytesize-24-dsl2-pipeline-starter';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-25-nf-core-ampliseq.md': {
                id: '2021/bytesize-25-nf-core-ampliseq.md';
                slug: '2021/bytesize-25-nf-core-ampliseq';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-26-nf-core-metaboigniter.md': {
                id: '2021/bytesize-26-nf-core-metaboigniter.md';
                slug: '2021/bytesize-26-nf-core-metaboigniter';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-27-nf-core-smrna.md': {
                id: '2021/bytesize-27-nf-core-smrna.md';
                slug: '2021/bytesize-27-nf-core-smrna';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-28-nf-core-sarek.md': {
                id: '2021/bytesize-28-nf-core-sarek.md';
                slug: '2021/bytesize-28-nf-core-sarek';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-29-nf-core-coproid.md': {
                id: '2021/bytesize-29-nf-core-coproid.md';
                slug: '2021/bytesize-29-nf-core-coproid';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-3-code-structure.md': {
                id: '2021/bytesize-3-code-structure.md';
                slug: '2021/bytesize-3-code-structure';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-30-nf-core-pgdb.md': {
                id: '2021/bytesize-30-nf-core-pgdb.md';
                slug: '2021/bytesize-30-nf-core-pgdb';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-4-github-contribution-basics.md': {
                id: '2021/bytesize-4-github-contribution-basics.md';
                slug: '2021/bytesize-4-github-contribution-basics';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-5-dsl2-module-development.md': {
                id: '2021/bytesize-5-dsl2-module-development.md';
                slug: '2021/bytesize-5-dsl2-module-development';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-6-dsl2-using-modules-pipelines.md': {
                id: '2021/bytesize-6-dsl2-using-modules-pipelines.md';
                slug: '2021/bytesize-6-dsl2-using-modules-pipelines';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-7-nf-core-ci-tests.md': {
                id: '2021/bytesize-7-nf-core-ci-tests.md';
                slug: '2021/bytesize-7-nf-core-ci-tests';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-8-nf-core-offline.md': {
                id: '2021/bytesize-8-nf-core-offline.md';
                slug: '2021/bytesize-8-nf-core-offline';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/bytesize-9-json-schema.md': {
                id: '2021/bytesize-9-json-schema.md';
                slug: '2021/bytesize-9-json-schema';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/elixir-workflow-training-event.md': {
                id: '2021/elixir-workflow-training-event.md';
                slug: '2021/elixir-workflow-training-event';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/genomeweb.md': {
                id: '2021/genomeweb.md';
                slug: '2021/genomeweb';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/hackathon-march-2021.md': {
                id: '2021/hackathon-march-2021.md';
                slug: '2021/hackathon-march-2021';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/hackathon-october-2021.md': {
                id: '2021/hackathon-october-2021.md';
                slug: '2021/hackathon-october-2021';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2021/seqera-cloud-series-azure.md': {
                id: '2021/seqera-cloud-series-azure.md';
                slug: '2021/seqera-cloud-series-azure';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-31-nf-core-dualrnaseq.md': {
                id: '2022/bytesize-31-nf-core-dualrnaseq.md';
                slug: '2022/bytesize-31-nf-core-dualrnaseq';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-32-nf-core-rnaseq.md': {
                id: '2022/bytesize-32-nf-core-rnaseq.md';
                slug: '2022/bytesize-32-nf-core-rnaseq';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-33-tower-cli.md': {
                id: '2022/bytesize-33-tower-cli.md';
                slug: '2022/bytesize-33-tower-cli';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-34-dsl2-syntax.md': {
                id: '2022/bytesize-34-dsl2-syntax.md';
                slug: '2022/bytesize-34-dsl2-syntax';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-35-debugging-pipeline.md': {
                id: '2022/bytesize-35-debugging-pipeline.md';
                slug: '2022/bytesize-35-debugging-pipeline';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-36-multiqc.md': {
                id: '2022/bytesize-36-multiqc.md';
                slug: '2022/bytesize-36-multiqc';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-37-gathertown.md': {
                id: '2022/bytesize-37-gathertown.md';
                slug: '2022/bytesize-37-gathertown';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-40-software-packaging.md': {
                id: '2022/bytesize-40-software-packaging.md';
                slug: '2022/bytesize-40-software-packaging';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-41-prettier.md': {
                id: '2022/bytesize-41-prettier.md';
                slug: '2022/bytesize-41-prettier';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-Subworkflows.md': {
                id: '2022/bytesize-Subworkflows.md';
                slug: '2022/bytesize-subworkflows';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-bactopia.md': {
                id: '2022/bytesize-bactopia.md';
                slug: '2022/bytesize-bactopia';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-chipseq.md': {
                id: '2022/bytesize-chipseq.md';
                slug: '2022/bytesize-chipseq';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-custom_scripts.md': {
                id: '2022/bytesize-custom_scripts.md';
                slug: '2022/bytesize-custom_scripts';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-cutandrun.md': {
                id: '2022/bytesize-cutandrun.md';
                slug: '2022/bytesize-cutandrun';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-dsl2-coding-styles-pt1.md': {
                id: '2022/bytesize-dsl2-coding-styles-pt1.md';
                slug: '2022/bytesize-dsl2-coding-styles-pt1';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-gitpod.md': {
                id: '2022/bytesize-gitpod.md';
                slug: '2022/bytesize-gitpod';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-inkscape.md': {
                id: '2022/bytesize-inkscape.md';
                slug: '2022/bytesize-inkscape';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-nanoseq.md': {
                id: '2022/bytesize-nanoseq.md';
                slug: '2022/bytesize-nanoseq';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-nascent.md': {
                id: '2022/bytesize-nascent.md';
                slug: '2022/bytesize-nascent';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-proteinfold.md': {
                id: '2022/bytesize-proteinfold.md';
                slug: '2022/bytesize-proteinfold';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-resources-to-learn-nextflow.md': {
                id: '2022/bytesize-resources-to-learn-nextflow.md';
                slug: '2022/bytesize-resources-to-learn-nextflow';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-rnafusion.md': {
                id: '2022/bytesize-rnafusion.md';
                slug: '2022/bytesize-rnafusion';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-survey_results.md': {
                id: '2022/bytesize-survey_results.md';
                slug: '2022/bytesize-survey_results';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-viralrecon.md': {
                id: '2022/bytesize-viralrecon.md';
                slug: '2022/bytesize-viralrecon';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize-working_with_github.md': {
                id: '2022/bytesize-working_with_github.md';
                slug: '2022/bytesize-working_with_github';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize_community-update.md': {
                id: '2022/bytesize_community-update.md';
                slug: '2022/bytesize_community-update';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize_configure_lint_tests.md': {
                id: '2022/bytesize_configure_lint_tests.md';
                slug: '2022/bytesize_configure_lint_tests';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize_nf-core-airrflow.md': {
                id: '2022/bytesize_nf-core-airrflow.md';
                slug: '2022/bytesize_nf-core-airrflow';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/bytesize_nftest.md': {
                id: '2022/bytesize_nftest.md';
                slug: '2022/bytesize_nftest';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/euro-faang-workshop.md': {
                id: '2022/euro-faang-workshop.md';
                slug: '2022/euro-faang-workshop';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/hackathon-march-2022.md': {
                id: '2022/hackathon-march-2022.md';
                slug: '2022/hackathon-march-2022';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/hackathon-october-2022.md': {
                id: '2022/hackathon-october-2022.md';
                slug: '2022/hackathon-october-2022';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2022/training-october-2022.md': {
                id: '2022/training-october-2022.md';
                slug: '2022/training-october-2022';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/bytesize_funcscan.md': {
                id: '2023/bytesize_funcscan.md';
                slug: '2023/bytesize_funcscan';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/bytesize_nf-core-taxprofiler.md': {
                id: '2023/bytesize_nf-core-taxprofiler.md';
                slug: '2023/bytesize_nf-core-taxprofiler';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/denmark-unseen-bio.md': {
                id: '2023/hackathon-march-2023/denmark-unseen-bio.md';
                slug: '2023/hackathon-march-2023/denmark-unseen-bio';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/france-igdr.md': {
                id: '2023/hackathon-march-2023/france-igdr.md';
                slug: '2023/hackathon-march-2023/france-igdr';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/germany-mpi-eva.md': {
                id: '2023/hackathon-march-2023/germany-mpi-eva.md';
                slug: '2023/hackathon-march-2023/germany-mpi-eva';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/germany-qbic.md': {
                id: '2023/hackathon-march-2023/germany-qbic.md';
                slug: '2023/hackathon-march-2023/germany-qbic';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/index.md': {
                id: '2023/hackathon-march-2023/index.md';
                slug: '2023/hackathon-march-2023';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/south-africa-stellenbosch.md': {
                id: '2023/hackathon-march-2023/south-africa-stellenbosch.md';
                slug: '2023/hackathon-march-2023/south-africa-stellenbosch';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/spain-crg.md': {
                id: '2023/hackathon-march-2023/spain-crg.md';
                slug: '2023/hackathon-march-2023/spain-crg';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/stockholm.md': {
                id: '2023/hackathon-march-2023/stockholm.md';
                slug: '2023/hackathon-march-2023/stockholm';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/sweden-scilifelab.md': {
                id: '2023/hackathon-march-2023/sweden-scilifelab.md';
                slug: '2023/hackathon-march-2023/sweden-scilifelab';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/uk-igc-edinburgh.md': {
                id: '2023/hackathon-march-2023/uk-igc-edinburgh.md';
                slug: '2023/hackathon-march-2023/uk-igc-edinburgh';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/uk-wellcome-campus.md': {
                id: '2023/hackathon-march-2023/uk-wellcome-campus.md';
                slug: '2023/hackathon-march-2023/uk-wellcome-campus';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/usa-san-jose.md': {
                id: '2023/hackathon-march-2023/usa-san-jose.md';
                slug: '2023/hackathon-march-2023/usa-san-jose';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/usa-university-texas.md': {
                id: '2023/hackathon-march-2023/usa-university-texas.md';
                slug: '2023/hackathon-march-2023/usa-university-texas';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/hackathon-march-2023/usa-university-wyoming.md': {
                id: '2023/hackathon-march-2023/usa-university-wyoming.md';
                slug: '2023/hackathon-march-2023/usa-university-wyoming';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
            '2023/training-march-2023.md': {
                id: '2023/training-march-2023.md';
                slug: '2023/training-march-2023';
                body: string;
                collection: 'events';
                data: InferEntrySchema<'events'>;
            };
        };
    };

    type ContentConfig = typeof import('./config');
}
