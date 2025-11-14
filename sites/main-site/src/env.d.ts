/// <reference path="../.astro/types.d.ts" />
/// <reference types="astro/client" />

declare module 'astro:content' {
  interface DataEntryMap {
    'advisories': {
      id: string;
      data: {
        title: string;
        subtitle: string;
        category: ('pipelines' | 'modules' | 'subworkflows' | 'configuration') | ('pipelines' | 'modules' | 'subworkflows' | 'configuration')[];
        type: ('known_regression' | 'incompatibility' | 'security' | 'performance' | 'data_corruption' | 'scientific_advice' | 'other') | ('known_regression' | 'incompatibility' | 'security' | 'performance' | 'data_corruption' | 'scientific_advice' | 'other')[];
        severity: 'low' | 'medium' | 'high' | 'critical';
        publishedDate: Date;
        reporter?: (string | Record<string, string>)[];
        reviewer?: (string | Record<string, string>)[];
        pipelines?: string[] | Array<{name: string; versions: string[]}>;
        modules?: string[];
        subworkflows?: string[];
        configuration?: string[];
        nextflowVersions?: string[];
        nextflowExecutors?: (
          | 'AWS Batch'
          | 'Azure Batch'
          | 'Bridge'
          | 'Flux Executor'
          | 'Google Cloud Batch'
          | 'Google Life Sciences'
          | 'HTCondor'
          | 'HyperQueue'
          | 'Kubernetes'
          | 'Local'
          | 'LSF'
          | 'Moab'
          | 'NQSII'
          | 'OAR'
          | 'PBS/Torque'
          | 'PBS Pro'
          | 'SGE'
          | 'SLURM'
        )[];
        softwareDependencies?: (
          | 'Apptainer'
          | 'Charliecloud'
          | 'Docker'
          | 'Podman'
          | 'Sarus'
          | 'Shifter'
          | 'Singularity'
          | 'Conda'
          | 'Spack'
          | 'Wave'
        )[] | Array<{
          name: 'Apptainer' | 'Charliecloud' | 'Docker' | 'Podman' | 'Sarus' | 'Shifter' | 'Singularity' | 'Conda' | 'Spack' | 'Wave';
          versions: string[];
        }>;
        references?: Array<{
          title: string;
          description: string;
          url: string;
        }>;
      };
    };
  }
}
