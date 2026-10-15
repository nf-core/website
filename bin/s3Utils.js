import { S3Client, ListObjectsV2Command, GetObjectCommand } from "@aws-sdk/client-s3";

/**
 * Creates an anonymous S3 client for public bucket access
 */
export function createS3Client() {
  return new S3Client({
    region: "eu-west-1",
    signer: { sign: async (request) => request },
    credentials: {
      accessKeyId: "",
      secretAccessKey: "",
    },
  });
}

/**
 * Fetches all objects with a given prefix from S3
 * @param {string} bucketName - S3 bucket name
 * @param {string} prefix - Object key prefix
 * @param {string} [delimiter] - Optional delimiter for folder-like structure
 * @returns {Promise<{bucketContents: Array, commonPrefixes: Array}>} Object containing Contents and CommonPrefixes arrays
 */
export async function fetchS3Objects(bucketName, prefix, delimiter) {
  const client = createS3Client();

  let bucketContents = [];
  let commonPrefixes = [];
  let isTruncated = true;
  let continuationToken;

  try {
    while (isTruncated) {
      const response = await client.send(
        new ListObjectsV2Command({
          Bucket: bucketName,
          Prefix: prefix,
          Delimiter: delimiter,
          ContinuationToken: continuationToken,
        }),
      );

      if (response.KeyCount === 0) {
        break;
      }

      if (response.Contents) {
        bucketContents.push(...response.Contents);
      }

      if (response.CommonPrefixes) {
        commonPrefixes.push(...response.CommonPrefixes);
      }

      isTruncated = response.IsTruncated ?? false;
      continuationToken = response.NextContinuationToken;
    }

    return { bucketContents, commonPrefixes };
  } catch (error) {
    console.warn(`Failed to fetch S3 objects for prefix ${prefix}:`, error);
    return { bucketContents: [], commonPrefixes: [] };
  }
}

/**
 * Fetches and parses CO2 footprint summary file from S3 for a given pipeline and commit SHA
 * @param {string} pipelineName - Name of the pipeline (e.g., 'rnaseq')
 * @param {string} commitSha - Git commit SHA for the release
 * @returns {Promise<{filename: string, co2e_emissions: number|null, energy_consumption: number|null}|null>} CO2 footprint data object with filename, co2e_emissions, and energy_consumption, or null if not found
 */
export async function fetchCO2FootprintFiles(pipelineName, commitSha) {
  console.log(`Fetching CO2 footprint files for ${pipelineName} ${commitSha}`);
  const bucketName = "nf-core-awsmegatests";
  const prefix = `${pipelineName}/results-${commitSha}/`;

  const { bucketContents } = await fetchS3Objects(bucketName, prefix);

  // Filter for CO2 footprint summary file matching the pattern
  const co2FootprintFile = bucketContents.find((file) => {
    const fileName = file.Key?.split("/").pop() || "";
    return /co2footprint_summary_.*\.txt$/.test(fileName);
  });
  console.log(`CO2 footprint file found: ${co2FootprintFile?.Key}`);

  if (!co2FootprintFile) {
    return null;
  }

  const client = createS3Client();

  try {
    const response = await client.send(
      new GetObjectCommand({
        Bucket: bucketName,
        Key: co2FootprintFile.Key,
      }),
    );

    // Read the stream and convert to string
    const bodyContents = await response.Body?.transformToString();

    if (!bodyContents) {
      return null;
    }

    // Parse CO2e emissions and energy consumption
    const co2Match = bodyContents.match(/CO2e emissions:\s*([\d.]+)\s*g/);
    const energyMatch = bodyContents.match(/Energy consumption:\s*([\d.]+)\s*Wh/);

    const fileName = co2FootprintFile.Key?.split("/").pop() || "";

    return {
      filename: fileName,
      co2e_emissions: co2Match ? parseFloat(co2Match[1]) : null,
      energy_consumption: energyMatch ? parseFloat(energyMatch[1]) : null,
    };
  } catch (error) {
    console.warn(`Failed to fetch CO2 footprint file ${co2FootprintFile.Key}:`, error);
    return null;
  }
}

// Runs of one megatest launch (e.g. a test matrix) start within seconds of each other, retries hours or days later
const MEGATEST_LAUNCH_WINDOW_MS = 60 * 60 * 1000;

/**
 * Gets the start time of a Nextflow run from its execution report file name
 * (e.g. `execution_report_2026-09-23_13-00-35.html`), falling back to the upload time
 * @param {{Key: string, LastModified: Date}} report - S3 object of the execution report
 * @returns {number} Start time in milliseconds since epoch
 */
function reportStartTime(report) {
  const match = report.Key.match(/execution_report_(\d{4}-\d{2}-\d{2})_(\d{2})-(\d{2})-(\d{2})\.html$/);
  return match
    ? Date.parse(`${match[1]}T${match[2]}:${match[3]}:${match[4]}Z`)
    : new Date(report.LastModified).getTime();
}

/**
 * Parses the run outcome and the first line of the error message from a Nextflow execution report
 * @param {string} html - (Beginning of the) execution report
 * @returns {{status: "success"|"failed"|"unknown", error?: string}} Run outcome
 */
function parseExecutionReport(html) {
  if (html.includes("Workflow execution completed successfully")) {
    return { status: "success" };
  }
  if (!html.includes("Workflow execution completed unsuccessfully")) {
    return { status: "unknown" };
  }
  const message = html.match(/full error message was:<\/p>\s*<pre[^>]*>([\s\S]*?)<\/pre>/)?.[1] ?? "";
  const error = message
    .replace(/<[^>]*>/g, "")
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g, "'")
    .replace(/&amp;/g, "&")
    // remove terminal colour codes, e.g. `\x1b[0;31m` or a bare `[0;31m`
    .replace(/\x1b?\[[0-9;]*m/g, "")
    .split("\n")
    .map((line) => line.trim())
    .find(Boolean);
  return { status: "failed", error: error?.slice(0, 200) ?? "" };
}

/**
 * Fetches the AWS megatest status for a given pipeline release from the Nextflow execution reports on S3.
 * Reports are looked up in `pipeline_info/` of the results folder, or of its direct subfolders if it has none
 * (pipelines like rnaseq run each test configuration with its own outdir). Only the runs of the latest launch
 * are kept, and reports already parsed in a previous build (same key and ETag) are not downloaded again.
 * @param {string} pipelineName - Name of the pipeline (e.g., 'rnaseq')
 * @param {string} tagName - Release tag the results belong to
 * @param {string} commitSha - Git commit SHA for the release
 * @param {object} [previous] - Megatest status from the previous build, used as cache
 * @returns {Promise<{tag_name: string, has_results: boolean, runs: Array<{report: string, etag: string|null, started: string, status: string, error?: string}>}>} Megatest status object
 */
export async function fetchMegatestStatus(pipelineName, tagName, commitSha, previous) {
  const bucketName = "nf-core-awsmegatests";
  const prefix = `${pipelineName}/results-${commitSha}/`;
  const status = { tag_name: tagName, has_results: false, runs: [] };

  const { bucketContents, commonPrefixes } = await fetchS3Objects(bucketName, prefix, "/");
  if (bucketContents.length === 0 && commonPrefixes.length === 0) {
    return status;
  }
  status.has_results = true;

  const subDirs = commonPrefixes.map((p) => p.Prefix);
  const runDirs = subDirs.includes(`${prefix}pipeline_info/`) ? [prefix] : subDirs;
  const reports = (
    await Promise.all(runDirs.map((dir) => fetchS3Objects(bucketName, `${dir}pipeline_info/execution_report`)))
  ).flatMap(({ bucketContents: reportFiles }) => reportFiles);
  if (reports.length === 0) {
    return status;
  }

  const latestStart = Math.max(...reports.map(reportStartTime));
  const latestLaunch = reports.filter((report) => latestStart - reportStartTime(report) < MEGATEST_LAUNCH_WINDOW_MS);

  const client = createS3Client();
  status.runs = await Promise.all(
    latestLaunch.map(async (report) => {
      const run = { report: report.Key, etag: report.ETag, started: new Date(reportStartTime(report)).toISOString() };
      const cached = previous?.runs?.find((r) => r.report === report.Key && r.etag === report.ETag);
      if (cached) {
        return cached;
      }
      try {
        // The status banner is near the top of the report, so only fetch the first 512 kB
        const response = await client.send(
          new GetObjectCommand({ Bucket: bucketName, Key: report.Key, Range: "bytes=0-524287" }),
        );
        return { ...run, ...parseExecutionReport((await response.Body?.transformToString()) ?? "") };
      } catch (error) {
        console.warn(`Failed to fetch execution report ${report.Key}:`, error);
        // Without ETag, so the report is fetched again in the next build
        return { ...run, etag: null, status: "unknown" };
      }
    }),
  );

  return status;
}
