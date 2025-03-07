/// <reference path="../.astro/types.d.ts" />
/// <reference types="astro/client" />

// Define the content collection types
declare module 'astro:content' {
  interface DataEntryMap {
    'pipeline-results': {
      id: string;
      data: {
        pipeline: string;
        version: string;
        results_path: string;
        content: {
          Key: string;
          LastModified: Date;
          ETag: string;
          Size: number;
          StorageClass: string;
        }[];
        commonPrefixes: {
          Prefix: string;
        }[];
      };
    };
  }
}
